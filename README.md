# TEfastR

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24238942.svg)](https://doi.org/10.5281/zenodo.24238942)

Fast, vectorized gene and transposable element (TE) read counting in R.

TEfastR re-implements the read-counting stage of
[TEtranscripts](https://github.com/mhammell-laboratory/TEtranscripts)
(Jin et al. 2015, [doi:10.1093/bioinformatics/btv422](https://doi.org/10.1093/bioinformatics/btv422))
using Bioconductor. It counts reads in a BAM file against gene and TE
annotations and uses an EM algorithm to split multi-mapped reads across TE
instances. The output is the same gene/TE count table that TEtranscripts
passes to DESeq2.

TEtranscripts handles one alignment at a time in Python. TEfastR works on the
whole BAM at once:

- **One bulk BAM read** with `GenomicAlignments::readGAlignments()`
- **Bulk overlap queries** with `GenomicRanges::findOverlaps()`, run for all
  reads together
- **Sparse-matrix SQUAREM EM** (`Matrix`) to redistribute multi-mapped reads

It follows TEtranscripts' assignment rules: unique vs. multi-mapped
classification, `--stranded` matching, proportional ambiguity resolution, and
the same SQUAREM-accelerated EM.

## Performance and validation

On a real paired-end BAM with 44.5 million alignments (human, GENCODE v50 +
RepeatMasker), TEfastR loaded the GTFs, counted the reads, and ran 100 EM
iterations in **7.3 minutes on a single core**. Peak memory was about 54 GB.

On a chromosome 14 subset (3,894 genes and TE elements, 6.4M alignments),
per-feature counts correlate with Python `TEcount` at **r = 0.99999**, and 91%
of features agree to within 1%.

The test suite compares TEfastR with TEtranscripts' own counting core on
synthetic single-end and paired-end data. It covers unique assignment,
ambiguity, EM, pairing, all three strandedness modes, and fragment-length
estimation.

> As with any reimplementation of a published method, spot-check TEfastR
> against TEtranscripts on a subset of your data before publishing results.

## Installation

TEfastR depends on Bioconductor packages, so install those first:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("GenomicAlignments", "GenomicRanges", "IRanges",
                       "Rsamtools", "S4Vectors"))

remotes::install_github("sulevk/TEfastR", build_vignettes = TRUE)
```

Requires R >= 4.1.

## Quick start

Build the annotation indexes once per project. They don't depend on the BAM,
so you can reuse them for every sample.

```r
library(TEfastR)

te_gtf   <- system.file("extdata", "te.gtf",   package = "TEfastR")
gene_gtf <- system.file("extdata", "gene.gtf", package = "TEfastR")
bam      <- system.file("extdata", "test.bam", package = "TEfastR")

te_index   <- build_te_annotation(te_gtf)
gene_index <- build_gene_annotation(gene_gtf)

result <- count_te(bam, gene_index, te_index,
                   stranded = "no", mode = "multi", iteration = 100)
result
head(as.data.frame(result))

write_count_table(result, "sample1.cntTable")   # TEtranscripts .cntTable format
```

`count_te()` returns a `TEfastR_result` with:

| Element              | Contents                                                               |
|----------------------|------------------------------------------------------------------------|
| `gene_counts`        | Counts per gene                                                        |
| `te_instance_counts` | Counts per individual TE instance                                      |
| `te_element_counts`  | TE instances summed by `gene_id:family_id:class_id` (TEtranscripts' output level) |
| `stats`              | Unique, non-unique, and unannotated read counts, and average read length |

### Main arguments

| Argument     | Values                            | Meaning                                                    |
|--------------|-----------------------------------|------------------------------------------------------------|
| `stranded`   | `"no"`, `"forward"`, `"reverse"`  | Library strandedness (TEtranscripts `--stranded`)          |
| `mode`       | `"multi"` (default), `"uniq"`     | Split multi-mappers with EM, or count unique reads only    |
| `iteration`  | integer, default `100`            | EM iterations. `0` keeps the naive even split              |
| `fragLength` | integer                           | Fragment length for single-end data. Estimated automatically for paired-end |
| `maxL`       | integer, default `500`            | Paired-end pairs with a larger insert are left out of the fragment-length estimate |

## Multiple samples and DESeq2

`count_te_batch()` counts several BAMs and returns a feature × sample matrix
that you can pass straight to DESeq2:

```r
bams <- c(treatment1 = "treatment_rep1.bam", treatment2 = "treatment_rep2.bam",
          control1   = "control_rep1.bam",   control2   = "control_rep2.bam")

counts <- count_te_batch(bams, gene_index, te_index, stranded = "no", mode = "multi")

library(DESeq2)
coldata <- data.frame(group = factor(c("treatment", "treatment", "control", "control"),
                                     levels = c("control", "treatment")),
                      row.names = colnames(counts))
dds <- DESeqDataSetFromMatrix(round(counts), coldata, design = ~ group)
res <- results(DESeq(dds))
```

The counts need `round()` because EM makes TE counts fractional.
TEtranscripts' generated R script rounds them the same way.

## Per-locus TE counts with coordinates

TEfastR counts every TE copy in the GTF before summing copies per element, so
the per-copy counts come with their genomic locations:

```r
loci <- te_locus_counts(result, te_index)      # loci with count > 0
head(loci)
#>                 locus_id chrom start end strand ... subfamily family class length count
write_locus_table(result, te_index, "sample.locus.tsv")

# feature x sample matrix with one row per TE copy (TElocal-style), for DESeq2
locus_counts <- count_te_batch(bams, gene_index, te_index, level = "locus")
coords <- te_index$loci   # rows match rownames(locus_counts) after the genes
```

Locus IDs use TElocal's format, `transcript_id:gene_id:family_id:class_id`
(for example `L1M2_dup229:L1M2:L1:LINE`). Coordinates are 1-based and
inclusive, as in the GTF. Treat counts for young, near-identical copies as
EM estimates. Counts for older, more divergent copies, which get mostly
uniquely mapped reads, are much more reliable.

## Command line

A command-line wrapper that mirrors TEtranscripts' `TEcount` is installed
with the package:

```sh
SCRIPT=$(Rscript -e 'cat(system.file("scripts", "run_te_count.R", package = "TEfastR"))')

Rscript "$SCRIPT" --bam sample.bam --gtf genes.gtf --te te.gtf \
  --stranded no --mode multi --iteration 100 --out sample.cntTable
```

Other options: `--fragLength`, `--feature` (default `exon`),
`--id_attribute` (default `gene_id`), and `--locus_out sample.locus.tsv`
(also write per-locus TE counts with coordinates).

## Input formats

- **BAM**: coordinate- or name-sorted, indexed or not. Unlike TEtranscripts,
  TEfastR has no `--sortByPos` step. Reads are grouped by name with a hash
  join, so the BAM doesn't need to be name-sorted first. A BAM must be all
  single-end or all paired-end.
- **Gene GTF**: a standard GTF, for example GENCODE. Uses `exon` features and
  the `gene_id` attribute by default.
- **TE GTF**: a TEtranscripts-style TE GTF, where each record has `gene_id`,
  `transcript_id`, `family_id`, and `class_id` attributes.

## Known differences from TEtranscripts

- CIGAR `D` (deletion) keeps a read in one aligned block. TEtranscripts
  splits the read at every CIGAR operation. Only `N` (splice) splits a block
  in both tools. This matters only for reads with a small indel right next to
  a feature boundary.
- Single-end or paired-end is decided once from the first alignment, not per
  read.

The vignette explains these differences in more detail. It also describes
the bugs that real-data validation found and how they were fixed:

```r
vignette("TEfastR")
```

## How to cite

If you use TEfastR, please cite it with its Zenodo DOI:

> Koks, S. TEfastR: fast, vectorized gene and transposable element read
> counting in R. Zenodo. https://doi.org/10.5281/zenodo.24238942

| DOI | Use |
|---|---|
| [10.5281/zenodo.24238942](https://doi.org/10.5281/zenodo.24238942) | **All versions.** Always resolves to the latest release; use this to cite the package in general |
| [10.5281/zenodo.24238943](https://doi.org/10.5281/zenodo.24238943) | **v0.2.0 only.** Use this when a methods section must name the exact version used |

In R, `citation("TEfastR")` prints the reference (and the TEtranscripts paper
below). GitHub's "Cite this repository" button gives it in APA and BibTeX.
Every later release gets its own version DOI on
[Zenodo](https://doi.org/10.5281/zenodo.24238942).

Please also cite the original TEtranscripts paper:

> Jin Y, Tam OH, Paniagua E, Hammell M. TEtranscripts: a package for including
> transposable elements in differential expression analysis of RNA-seq
> datasets. *Bioinformatics* 31(22):3593–3599 (2015).
> doi:10.1093/bioinformatics/btv422

## License

Artistic-2.0
