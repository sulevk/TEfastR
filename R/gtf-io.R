.gtf_cols <- c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attr")

#' @importFrom data.table fread
#' @importFrom stringi stri_match_first_regex
.extract_attr <- function(attr_col, key) {
  m <- stringi::stri_match_first_regex(attr_col, paste0(key, ' "([^"]*)"'))
  m[, 2]
}

#' Read a TE annotation GTF
#'
#' Parses a TEtranscripts-format TE GTF (one row per TE instance, with
#' `gene_id`, `transcript_id`, `family_id`, and `class_id` attributes -- the
#' same format used by `TEtranscripts --TE`). Mirrors
#' `TEToolkit/TEindex.py`'s `TEfeatures.build()`.
#'
#' @param path Path to the TE GTF file.
#' @param feature_type GTF feature (column 3) to keep. Default `"exon"`.
#' @return A [data.table::data.table()] with one row per TE instance.
#' @export
read_te_gtf <- function(path, feature_type = "exon") {
  dt <- data.table::fread(path, sep = "\t", header = FALSE, quote = "",
                           col.names = .gtf_cols, showProgress = FALSE)
  dt <- dt[feature == feature_type]
  if (nrow(dt) == 0) stop("No '", feature_type, "' features found in TE GTF.")
  dt[, `:=`(
    transcript_id = .extract_attr(attr, "transcript_id"),
    gene_id       = .extract_attr(attr, "gene_id"),
    family_id     = .extract_attr(attr, "family_id"),
    class_id      = .extract_attr(attr, "class_id")
  )]
  bad <- is.na(dt$transcript_id) | is.na(dt$gene_id) | is.na(dt$family_id) | is.na(dt$class_id)
  if (any(bad))
    stop(sprintf("TE GTF format error: %d line(s) missing gene_id/transcript_id/family_id/class_id.",
                  sum(bad)))
  dt[, ele_name := paste(gene_id, family_id, class_id, sep = ":")]
  dt[, te_idx := .I]
  dt[, attr := NULL]
  dt[]
}

#' Read a gene annotation GTF
#'
#' Parses a standard gene-model GTF (e.g. GENCODE), keeping one row per
#' feature of type `feature_type` and extracting `id_attribute` as the
#' gene identifier. Mirrors `TEToolkit/GeneFeatures.py`.
#'
#' @param path Path to the gene GTF file.
#' @param feature_type GTF feature (column 3) to keep. Default `"exon"`.
#' @param id_attribute Attribute key used as the gene identifier. Default
#'   `"gene_id"`.
#' @return A [data.table::data.table()] with one row per feature.
#' @export
read_gene_gtf <- function(path, feature_type = "exon", id_attribute = "gene_id") {
  dt <- data.table::fread(path, sep = "\t", header = FALSE, quote = "",
                           col.names = .gtf_cols, showProgress = FALSE)
  dt <- dt[feature == feature_type]
  if (nrow(dt) == 0) stop("No '", feature_type, "' features found in gene GTF.")
  dt[, gene := .extract_attr(attr, id_attribute)]
  dt <- dt[!is.na(gene)]
  dt[, attr := NULL]
  dt[]
}
