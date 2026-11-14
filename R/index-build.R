#' Build a TE annotation index
#'
#' Parses a TE GTF and builds the [GenomicRanges::GRanges()] index used by
#' [count_te()]. Build this once per project and reuse it across every BAM
#' file -- it does not depend on the BAM.
#'
#' @param te_gtf Path to the TE GTF file.
#' @param feature_type GTF feature (column 3) to keep. Default `"exon"`.
#' @return A `TEfastR_te_index` object (a list); pass it to [count_te()].
#'   Its `loci` element is a data frame with one row per TE instance (genomic
#'   locus), in the same order as `count_te()`'s `te_instance_counts`:
#'   `locus_id` (TElocal-style `transcript_id:gene_id:family_id:class_id`),
#'   `chrom`, `start`, `end` (1-based, inclusive, as in the GTF), `strand`,
#'   `transcript_id`, `subfamily` (`gene_id`), `family`, `class`, `length`.
#' @examples
#' te_gtf <- system.file("extdata", "te.gtf", package = "TEfastR")
#' te_index <- build_te_annotation(te_gtf)
#' te_index$n
#' @export
build_te_annotation <- function(te_gtf, feature_type = "exon") {
  te_dt <- read_te_gtf(te_gtf, feature_type)
  gr <- GenomicRanges::GRanges(te_dt$chrom, IRanges::IRanges(te_dt$start, te_dt$end),
                                strand = te_dt$strand)
  elements <- sort(unique(te_dt$ele_name))
  len <- te_dt$end - te_dt$start + 1L
  # TElocal's locus ID format; made unique in case a GTF reuses transcript_ids
  locus_id <- make.unique(paste(te_dt$transcript_id, te_dt$ele_name, sep = ":"), sep = "_")
  loci <- data.frame(locus_id = locus_id, chrom = te_dt$chrom, start = te_dt$start,
                     end = te_dt$end, strand = te_dt$strand,
                     transcript_id = te_dt$transcript_id, subfamily = te_dt$gene_id,
                     family = te_dt$family_id, class = te_dt$class_id, length = len,
                     stringsAsFactors = FALSE)
  structure(list(
    gr = gr,
    strand = te_dt$strand,
    ele_idx = match(te_dt$ele_name, elements),
    elements = elements,
    length = len,
    n = nrow(te_dt),
    loci = loci
  ), class = "TEfastR_te_index")
}

#' Build a gene annotation index
#'
#' Parses a gene GTF and builds the three strand-bucketed
#' [GenomicRanges::GRanges()] indexes used by [count_te()] (genes on `"+"`,
#' genes on `"-"`, and all genes combined for unstranded reads -- mirroring
#' TEtranscripts' own strand-matching rule). Build this once per project and
#' reuse it across every BAM file.
#'
#' @param gene_gtf Path to the gene GTF file.
#' @param feature_type GTF feature (column 3) to keep. Default `"exon"`.
#' @param id_attribute Attribute key used as the gene identifier. Default
#'   `"gene_id"`.
#' @return A `TEfastR_gene_index` object (a list); pass it to [count_te()].
#' @examples
#' gene_gtf <- system.file("extdata", "gene.gtf", package = "TEfastR")
#' gene_index <- build_gene_annotation(gene_gtf)
#' length(gene_index$genes)
#' @export
build_gene_annotation <- function(gene_gtf, feature_type = "exon", id_attribute = "gene_id") {
  gene_dt <- read_gene_gtf(gene_gtf, feature_type, id_attribute)
  genes <- sort(unique(gene_dt$gene))
  gidx <- match(gene_dt$gene, genes)
  gr_all <- GenomicRanges::GRanges(gene_dt$chrom, IRanges::IRanges(gene_dt$start, gene_dt$end),
                                    strand = "*")
  structure(list(
    genes = genes,
    all   = gr_all,
    plus  = gr_all[gene_dt$strand == "+"],
    minus = gr_all[gene_dt$strand == "-"],
    gene_idx_all   = gidx,
    gene_idx_plus  = gidx[gene_dt$strand == "+"],
    gene_idx_minus = gidx[gene_dt$strand == "-"]
  ), class = "TEfastR_gene_index")
}
