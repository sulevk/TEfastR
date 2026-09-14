# Overlaps every CIGAR block of every candidate alignment against the gene
# index in 3 bulk calls (unstranded/"+"/"-" reads), instead of one interval
# query per read. Returns a data.table(locus_id, gene_idx) with one row per
# distinct (locus, gene) hit -- the R analogue of Gene_annotation()'s
# per-locus set() union in TEtranscripts.
.overlap_genes <- function(blk_gr, blk_locus, blk_dir, gene_idx) {
  hits_list <- list()

  do_one <- function(sel, target_gr, target_idx) {
    if (!any(sel) || length(target_gr) == 0) return(NULL)
    h <- GenomicRanges::findOverlaps(blk_gr[sel], target_gr, ignore.strand = TRUE)
    if (length(h) == 0) return(NULL)
    data.table::data.table(locus_id = blk_locus[sel][S4Vectors::queryHits(h)],
                            gene_idx = target_idx[S4Vectors::subjectHits(h)])
  }

  hits_list[[1]] <- do_one(blk_dir == 0L, gene_idx$all,   gene_idx$gene_idx_all)
  hits_list[[2]] <- do_one(blk_dir == 1L, gene_idx$plus,  gene_idx$gene_idx_plus)
  hits_list[[3]] <- do_one(blk_dir == -1L, gene_idx$minus, gene_idx$gene_idx_minus)

  out <- data.table::rbindlist(hits_list)
  if (nrow(out) == 0) return(data.table::data.table(locus_id = integer(0), gene_idx = integer(0)))
  unique(out)
}

# Overlaps every CIGAR block against the (strand-unsplit) TE index in one
# bulk call, then keeps a hit only when the read is unstranded or its
# direction matches the TE's own strand -- mirrors TEfeatures.TE_annotation().
.overlap_tes <- function(blk_gr, blk_locus, blk_dir, te_idx) {
  if (length(te_idx$gr) == 0)
    return(data.table::data.table(locus_id = integer(0), te_idx = integer(0)))
  h <- GenomicRanges::findOverlaps(blk_gr, te_idx$gr, ignore.strand = TRUE)
  if (length(h) == 0)
    return(data.table::data.table(locus_id = integer(0), te_idx = integer(0)))
  dir <- blk_dir[S4Vectors::queryHits(h)]
  te_strand_char <- te_idx$strand[S4Vectors::subjectHits(h)]
  read_strand_char <- c("-", ".", "+")[dir + 2L]
  keep <- read_strand_char == "." | read_strand_char == te_strand_char
  out <- data.table::data.table(locus_id = blk_locus[S4Vectors::queryHits(h)][keep],
                                 te_idx = S4Vectors::subjectHits(h)[keep])
  unique(out)
}
