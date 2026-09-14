# One row per (read, locus). A "locus" is one r1/r2 alignment-pair candidate
# for that read; row1/row2 index into `aln`. Mirrors the pairing rules in
# TEToolkit/Counting.py's main while-loop, including its handling of an
# unequal, both-nonzero count of read1/read2 records for one multi-mapped
# read name: such a read contributes no loci at all (silently counted as
# unannotated), preserved here for parity with TEtranscripts.
.build_loci <- function(at, is_paired, mode) {
  if (!is_paired) {
    keep <- if (mode == "uniq") at[is_unique_group == TRUE] else at
    return(data.table::data.table(rkey = keep$rkey, row1 = keep$row, row2 = NA_integer_))
  }

  uniq_grp <- at[is_unique_group == TRUE]
  # exactly one representative locus per unique-group key
  u1 <- uniq_grp[is_read1 == TRUE, .(rkey, row1 = row)]
  u2 <- uniq_grp[is_read2 == TRUE, .(rkey, row2 = row)]
  uniq_keys <- unique(uniq_grp$rkey)
  uniq_loci <- merge(data.table::data.table(rkey = uniq_keys), u1, by = "rkey", all.x = TRUE, sort = FALSE)
  uniq_loci <- merge(uniq_loci, u2, by = "rkey", all.x = TRUE, sort = FALSE)

  if (mode == "uniq") return(uniq_loci[, .(rkey, row1, row2)])

  multi_grp <- at[is_unique_group == FALSE]
  if (nrow(multi_grp) == 0) return(uniq_loci[, .(rkey, row1, row2)])

  g1 <- multi_grp[is_read1 == TRUE][order(group_row)]
  g2 <- multi_grp[is_read2 == TRUE][order(group_row)]
  g1[, idx := seq_len(.N), by = rkey]
  g2[, idx := seq_len(.N), by = rkey]
  n1 <- g1[, .N, by = rkey]; data.table::setnames(n1, "N", "n1")
  n2 <- g2[, .N, by = rkey]; data.table::setnames(n2, "N", "n2")
  cnt <- merge(n1, n2, by = "rkey", all = TRUE)
  cnt[is.na(n1), n1 := 0L]; cnt[is.na(n2), n2 := 0L]

  only1 <- cnt[n2 == 0L, rkey]
  only2 <- cnt[n1 == 0L, rkey]
  matched <- cnt[n1 == n2 & n1 > 0L, rkey]
  # n1 != n2 and both > 0: silently dropped, as in Counting.py

  multi_loci <- data.table::rbindlist(list(
    if (length(only1) > 0) g1[rkey %in% only1, .(rkey, row1 = row, row2 = NA_integer_)] else NULL,
    if (length(only2) > 0) g2[rkey %in% only2, .(rkey, row1 = NA_integer_, row2 = row)] else NULL,
    if (length(matched) > 0) {
      a <- g1[rkey %in% matched]; b <- g2[rkey %in% matched]
      m <- merge(a[, .(rkey, idx, row1 = row)], b[, .(rkey, idx, row2 = row)], by = c("rkey", "idx"))
      m[, .(rkey, row1, row2)]
    } else NULL
  ), use.names = TRUE)

  data.table::rbindlist(list(uniq_loci[, .(rkey, row1, row2)], multi_loci), use.names = TRUE)
}

# CIGAR-derived aligned blocks (via grglist) for every locus, unioning r1's
# and r2's blocks, tagged with this locus's overlap direction. Mirrors
# fetch_exon() + the direction logic at the top of ovp_annotation() in
# TEToolkit/Counting.py.
#
# Known difference from TEtranscripts: CIGAR "D" (deletion) stays within one
# contiguous aligned block here (via GenomicAlignments::grglist()), whereas
# fetch_exon() starts a new block at every CIGAR op, splitting the read into
# extra tiny fragments around small deletions. Only "N" (splice/gap) breaks a
# block in both implementations.
.build_locus_blocks <- function(loci, aln, stranded) {
  rows_needed <- stats::na.omit(unique(c(loci$row1, loci$row2)))
  grl <- GenomicAlignments::grglist(aln[rows_needed])
  flat <- unlist(grl, use.names = FALSE)
  blk_row <- rep(rows_needed, S4Vectors::elementNROWS(grl))

  strand1 <- as.character(GenomicAlignments::strand(aln))  # per-row own strand
  dir_r1 <- rep(0L, nrow(loci))
  has1 <- !is.na(loci$row1)
  dir_r1[has1] <- ifelse(strand1[loci$row1[has1]] == "-", -1L, 1L)
  dir_final <- ifelse(has1, dir_r1, 1L)
  has2 <- !is.na(loci$row2)
  # r2 present and r2 is forward ("+") => direction = -1 (mirrors ovp_annotation)
  idx2 <- which(has2)
  dir_final[idx2[strand1[loci$row2[idx2]] == "+"]] <- -1L

  if (stranded == "no") dir_final[] <- 0L
  if (stranded == "reverse") dir_final <- dir_final * -1L

  dir_by_row1 <- rep(NA_integer_, length(rows_needed)); names(dir_by_row1) <- as.character(rows_needed)
  locus_by_row <- rep(NA_integer_, length(rows_needed)); names(locus_by_row) <- as.character(rows_needed)
  if (any(has1)) { locus_by_row[as.character(loci$row1[has1])] <- loci$locus_id[has1]
                    dir_by_row1[as.character(loci$row1[has1])] <- dir_final[has1] }
  if (any(has2)) { locus_by_row[as.character(loci$row2[has2])] <- loci$locus_id[has2]
                    dir_by_row1[as.character(loci$row2[has2])] <- dir_final[has2] }

  blk_key <- as.character(blk_row)
  data.table::data.table(
    seqnames = as.character(GenomicAlignments::seqnames(flat)),
    start = GenomicAlignments::start(flat), end = GenomicAlignments::end(flat),
    locus_id = locus_by_row[blk_key],
    dir = dir_by_row1[blk_key]
  )
}

# `uniq_loci` must be the UNIQUE-read subset of loci (one locus per read) --
# Counting.py only accumulates this estimate from reads it has already
# classified as uniquely mapped, inside the `len(multi_read1) <= 1 and
# len(multi_read2) <= 1` branch. It also floor-divides the final average
# (`avgReadLength = avgReadLength // tmp_cnt`, Python integer division) and,
# for paired-end, discards pairs with an insert size over `max_l` before
# counting them toward its first-10000-reads cap. Both details matter: this
# average feeds TE effective-length = TE length - avg_read_length + 1, and
# for short/truncated TE instances that difference can sit close enough to
# zero that a fractional-vs-floored estimate flips it from included to
# excluded (or a skewed estimate from outlier insert sizes does the same).
.estimate_read_length <- function(aln, uniq_loci, is_paired, frag_length, max_l = 500L) {
  if (!is_paired) {
    if (frag_length > 0) return(frag_length)
    return(floor(mean(GenomicAlignments::qwidth(aln)[seq_len(min(10000L, length(aln)))])))
  }
  both <- uniq_loci[!is.na(row1) & !is.na(row2)]
  if (nrow(both) > 0) {
    pos1 <- GenomicAlignments::start(aln)[both$row1]
    pos2 <- GenomicAlignments::start(aln)[both$row2]
    within_maxl <- abs(pos1 - pos2) <= max_l
    both <- both[within_maxl][seq_len(min(10000L, sum(within_maxl)))]
  }
  if (nrow(both) == 0) {
    if (frag_length > 0) return(frag_length)
    stop("Not enough paired reads to estimate fragment length; pass fragLength explicitly.")
  }
  pos1 <- GenomicAlignments::start(aln)[both$row1]
  pos2 <- GenomicAlignments::start(aln)[both$row2]
  qw2 <- GenomicAlignments::qwidth(aln)[both$row2]
  floor(mean(abs(pos1 - pos2) + qw2))
}
