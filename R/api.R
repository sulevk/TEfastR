#' Count gene and TE reads in one BAM file
#'
#' The main entry point. Counts a single BAM file against a gene index (from
#' [build_gene_annotation()]) and a TE index (from [build_te_annotation()]),
#' reproducing TEtranscripts' `TEcount`/`TEtranscripts` counting semantics:
#' unique-read gene/TE assignment with proportional ambiguity resolution, and
#' a SQUAREM EM step that redistributes multi-mapped reads across TE
#' instances.
#'
#' Unlike TEtranscripts, no `sortByPos`/name-sort step is needed: reads are
#' grouped by name via a hash join, so a coordinate-sorted, indexed BAM is
#' read directly, in whatever order it's stored in.
#'
#' @param bam Path to a BAM file (coordinate- or name-sorted; indexed or not).
#' @param gene_index A `TEfastR_gene_index`, from [build_gene_annotation()].
#' @param te_index A `TEfastR_te_index`, from [build_te_annotation()].
#' @param stranded One of `"no"`, `"forward"`, `"reverse"` -- library
#'   strandedness, as in TEtranscripts' `--stranded`.
#' @param mode One of `"multi"` (distribute multi-mapped reads via EM,
#'   TEtranscripts' default) or `"uniq"` (count unique mappers only).
#' @param iteration Number of SQUAREM EM iterations. `0` disables EM (TE
#'   multi-mappers get their naive, evenly-split prior instead).
#' @param fragLength For single-end data, the average fragment length to use
#'   for TE effective-length correction. Ignored (and estimated from the
#'   data) for paired-end BAMs unless there aren't enough proper pairs to
#'   estimate it, in which case this value is required.
#' @param maxL For paired-end data, proper pairs with an insert size over
#'   this are excluded from the fragment-length estimate used above (as
#'   likely discordant/mismapped pairs). Matches TEtranscripts' `--maxL`
#'   (default `500`).
#' @return A `TEfastR_result` object (a list) with elements `gene_counts`
#'   (named numeric vector), `te_instance_counts` (named numeric vector, one
#'   per TE instance), `te_element_counts` (named numeric vector, TE
#'   instances summed to elements), and `stats` (a list with `uniq_reads`,
#'   `nonunique`, `empty`, `avg_read_length`). Pass it to
#'   [write_count_table()] or [as.data.frame.TEfastR_result()].
#' @examples
#' te_gtf <- system.file("extdata", "te.gtf", package = "TEfastR")
#' gene_gtf <- system.file("extdata", "gene.gtf", package = "TEfastR")
#' bam <- system.file("extdata", "test.bam", package = "TEfastR")
#' te_index <- build_te_annotation(te_gtf)
#' gene_index <- build_gene_annotation(gene_gtf)
#' result <- count_te(bam, gene_index, te_index)
#' head(result$gene_counts)
#' @export
count_te <- function(bam, gene_index, te_index, stranded = c("no", "forward", "reverse"),
                      mode = c("multi", "uniq"), iteration = 100L, fragLength = 0L, maxL = 500L) {
  stranded <- match.arg(stranded)
  mode <- match.arg(mode)
  stopifnot(inherits(gene_index, "TEfastR_gene_index"), inherits(te_index, "TEfastR_te_index"))

  res <- .count_transcript_abundance(bam, gene_index, te_index, stranded, mode,
                                      as.integer(iteration), as.integer(fragLength), as.integer(maxL))

  ele_counts <- as.numeric(tapply(res$te_counts, te_index$ele_idx, sum))
  ele_counts <- ele_counts[match(seq_along(te_index$elements), sort(unique(te_index$ele_idx)))]

  structure(list(
    gene_counts = stats::setNames(res$gene_counts, gene_index$genes),
    te_instance_counts = res$te_counts,
    te_element_counts = stats::setNames(ele_counts, te_index$elements),
    stats = list(uniq_reads = res$uniq_reads, nonunique = res$nonunique,
                 empty = res$empty, avg_read_length = res$avg_read_length),
    bam = bam
  ), class = "TEfastR_result")
}

#' @export
print.TEfastR_result <- function(x, ...) {
  cat("<TEfastR_result>", basename(x$bam), "\n")
  cat(sprintf("  %d genes, %d TE elements\n", length(x$gene_counts), length(x$te_element_counts)))
  cat(sprintf("  unique reads: %d  non-unique: %d  unannotated: %d\n",
              x$stats$uniq_reads, x$stats$nonunique, x$stats$empty))
  invisible(x)
}

#' Convert a `count_te()` result to a data frame
#'
#' @param x A `TEfastR_result`, from [count_te()].
#' @param ... Ignored.
#' @return A two-column data frame (`id`, `count`) with one row per gene and
#'   one row per TE element, in TEtranscripts' `.cntTable` row order (genes,
#'   then TE elements).
#' @export
as.data.frame.TEfastR_result <- function(x, ...) {
  data.frame(id = c(names(x$gene_counts), names(x$te_element_counts)),
             count = c(unname(x$gene_counts), unname(x$te_element_counts)),
             stringsAsFactors = FALSE)
}

#' Write a count table in TEtranscripts' `.cntTable` format
#'
#' @param result A `TEfastR_result`, from [count_te()].
#' @param path Output file path.
#' @return `path`, invisibly.
#' @export
write_count_table <- function(result, path) {
  out <- as.data.frame(result)
  data.table::fwrite(out, path, sep = "\t", col.names = FALSE)
  invisible(path)
}

#' Per-locus TE counts with genomic coordinates
#'
#' Joins a [count_te()] result's per-instance TE counts back to each
#' instance's genomic location from the TE index -- the TElocal-style view of
#' the same counts that [write_count_table()] reports summed per element.
#'
#' Per-locus counts for young, near-identical TE copies depend heavily on the
#' EM's split of multi-mapped reads and should be read as estimates;
#' uniquely-mappable (older, more divergent) loci are far more reliable.
#'
#' @param result A `TEfastR_result`, from [count_te()].
#' @param te_index The `TEfastR_te_index` used to produce `result`.
#' @param min_count Keep only loci with a count strictly above this value.
#'   Default `0` drops the (usually millions of) zero-count loci; use `-Inf`
#'   to keep every locus.
#' @return A data frame: `locus_id`, `chrom`, `start`, `end` (1-based,
#'   inclusive), `strand`, `transcript_id`, `subfamily`, `family`, `class`,
#'   `length`, `count`, ordered as in the TE GTF.
#' @examples
#' te_gtf <- system.file("extdata", "te.gtf", package = "TEfastR")
#' gene_gtf <- system.file("extdata", "gene.gtf", package = "TEfastR")
#' bam <- system.file("extdata", "test.bam", package = "TEfastR")
#' te_index <- build_te_annotation(te_gtf)
#' result <- count_te(bam, build_gene_annotation(gene_gtf), te_index)
#' te_locus_counts(result, te_index)
#' @export
te_locus_counts <- function(result, te_index, min_count = 0) {
  stopifnot(inherits(result, "TEfastR_result"), inherits(te_index, "TEfastR_te_index"))
  .check_locus_index(te_index, length(result$te_instance_counts))
  out <- te_index$loci
  out$count <- result$te_instance_counts
  out <- out[out$count > min_count, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Write per-locus TE counts with coordinates to a TSV file
#'
#' Writes [te_locus_counts()] as a tab-separated table with a header row.
#'
#' @inheritParams te_locus_counts
#' @param path Output file path.
#' @return `path`, invisibly.
#' @export
write_locus_table <- function(result, te_index, path, min_count = 0) {
  data.table::fwrite(te_locus_counts(result, te_index, min_count), path, sep = "\t",
                     quote = FALSE)
  invisible(path)
}

.check_locus_index <- function(te_index, n_counts) {
  if (is.null(te_index$loci))
    stop("This TE index has no locus table; rebuild it with build_te_annotation().")
  if (nrow(te_index$loci) != n_counts)
    stop("TE index does not match this result (", nrow(te_index$loci), " loci vs ",
         n_counts, " instance counts); pass the te_index used for counting.")
}

#' Count several BAM files and assemble a feature x sample count matrix
#'
#' Convenience wrapper around [count_te()] for multi-sample projects: builds
#' the annotation indexes once (implicitly, via the `gene_index`/`te_index`
#' you pass in) and reuses them across every BAM, then assembles a single
#' feature x sample matrix -- the same shape `DESeq2::DESeqDataSetFromMatrix`
#' expects.
#'
#' @param bams A named or unnamed character vector of BAM paths. If unnamed,
#'   sample names are derived from the file names (extension stripped).
#' @param gene_index A `TEfastR_gene_index`, from [build_gene_annotation()].
#' @param te_index A `TEfastR_te_index`, from [build_te_annotation()].
#' @param level `"element"` (default: TE instances summed per element, as in
#'   TEtranscripts) or `"locus"` (one row per TE instance, as in TElocal;
#'   row names are `te_index$loci$locus_id`, so coordinates can be joined
#'   back from `te_index$loci`).
#' @param ... Passed on to [count_te()] (`stranded`, `mode`, `iteration`,
#'   `fragLength`).
#' @return A numeric matrix, genes then TE elements (or TE loci) in rows,
#'   samples in columns.
#' @examples
#' te_gtf <- system.file("extdata", "te.gtf", package = "TEfastR")
#' gene_gtf <- system.file("extdata", "gene.gtf", package = "TEfastR")
#' bam <- system.file("extdata", "test.bam", package = "TEfastR")
#' te_index <- build_te_annotation(te_gtf)
#' gene_index <- build_gene_annotation(gene_gtf)
#' mat <- count_te_batch(c(sample1 = bam), gene_index, te_index)
#' dim(mat)
#' @export
count_te_batch <- function(bams, gene_index, te_index, ..., level = c("element", "locus")) {
  level <- match.arg(level)
  if (is.null(names(bams)) || any(names(bams) == "")) {
    names(bams) <- tools::file_path_sans_ext(basename(bams))
  }
  if (level == "locus") .check_locus_index(te_index, te_index$n)
  results <- lapply(bams, count_te, gene_index = gene_index, te_index = te_index, ...)
  te_ids <- if (level == "locus") te_index$loci$locus_id else names(results[[1]]$te_element_counts)
  te_of <- if (level == "locus") function(r) r$te_instance_counts else function(r) r$te_element_counts
  ids <- c(names(results[[1]]$gene_counts), te_ids)
  mat <- vapply(results, function(r) c(r$gene_counts, te_of(r)), numeric(length(ids)))
  rownames(mat) <- ids
  colnames(mat) <- names(bams)
  mat
}
