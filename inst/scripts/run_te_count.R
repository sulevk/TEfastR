#!/usr/bin/env Rscript
# Command-line wrapper around TEfastR::count_te(), for use outside R/RStudio
# (matching TEtranscripts' own TEcount CLI). Find this script's installed
# path with:
#   system.file("scripts", "run_te_count.R", package = "TEfastR")
#
# Usage:
#   Rscript run_te_count.R --bam sample.bam --gtf genes.gtf --te te.gtf \
#     [--stranded no|forward|reverse] [--mode multi|uniq] [--iteration 100] \
#     [--fragLength 0] [--feature exon] [--id_attribute gene_id] [--out sample.cntTable]

suppressPackageStartupMessages(library(TEfastR))

.parse_args <- function(argv) {
  args <- list(bam = NULL, gtf = NULL, te = NULL, stranded = "no", mode = "multi",
               iteration = "100", fragLength = "0", feature = "exon",
               id_attribute = "gene_id", out = NULL)
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--") || !(substring(key, 3) %in% names(args)))
      stop("Unknown or malformed argument: ", key)
    args[[substring(key, 3)]] <- argv[[i + 1]]
    i <- i + 2
  }
  if (is.null(args$bam) || is.null(args$gtf) || is.null(args$te))
    stop("Usage: --bam <file> --gtf <gene GTF> --te <TE GTF> ",
         "[--stranded no|forward|reverse] [--mode uniq|multi] [--iteration N] ",
         "[--fragLength N] [--feature exon] [--id_attribute gene_id] [--out path]")
  if (is.null(args$out))
    args$out <- paste0(tools::file_path_sans_ext(basename(args$bam)), ".cntTable")
  args
}

args <- .parse_args(commandArgs(trailingOnly = TRUE))

message("Loading TE annotation...")
te_index <- build_te_annotation(args$te, feature_type = args$feature)
message(sprintf("  %d TE instances, %d elements", te_index$n, length(te_index$elements)))

message("Loading gene annotation...")
gene_index <- build_gene_annotation(args$gtf, feature_type = args$feature,
                                     id_attribute = args$id_attribute)
message(sprintf("  %d genes", length(gene_index$genes)))

message("Reading alignments and counting...")
result <- count_te(args$bam, gene_index, te_index, stranded = args$stranded, mode = args$mode,
                    iteration = as.integer(args$iteration), fragLength = as.integer(args$fragLength))
print(result)

write_count_table(result, args$out)
message("Wrote ", args$out)
