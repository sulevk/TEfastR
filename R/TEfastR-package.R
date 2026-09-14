#' TEfastR: fast, vectorized gene/TE read counting
#'
#' @description
#' TEfastR re-implements the counting stage of TEtranscripts
#' (\url{https://github.com/mhammell-laboratory/TEtranscripts}) using
#' Bioconductor's C-backed bulk interval operations instead of a
#' per-alignment Python loop. See `vignette("TEfastR")` for a full guide.
#'
#' @section Key functions:
#' * [build_gene_annotation()] / [build_te_annotation()] -- parse GTFs once,
#'   reuse across samples.
#' * [count_te()] -- count one BAM file against a pair of annotation indexes.
#' * [count_te_batch()] -- count several BAM files and assemble a
#'   feature x sample count matrix, ready for `DESeq2`.
#' * [write_count_table()] -- write a `count_te()` result in TEtranscripts'
#'   `.cntTable` format.
#'
#' @keywords internal
"_PACKAGE"

## data.table's NSE columns trigger R CMD check "no visible binding" notes;
## declare them once here rather than sprinkling globalVariables() calls.
utils::globalVariables(c(
  "feature", "attr", "transcript_id", "gene_id", "family_id", "class_id",
  "ele_name", "te_idx", "gene", "gene_idx", "rkey", "is_read1", "is_read2",
  "group_row", "n1", "n2", "is_unique_group", "nloci", "N", "locus_id",
  "row1", "row2", "idx", "w", "n_cand", "n_te_loci", "read_id", "mult",
  ".I", ".N", ".", ":="
))
