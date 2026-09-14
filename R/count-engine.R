# Internal counting engine. Faithful port of
# TEToolkit/Counting.py::count_transcript_abundance(), but vectorized: one
# bulk BAM read + bulk overlap calls for ALL reads, instead of a per-read
# Python loop; and a sparse-matrix EM instead of a pure-Python one. See
# build_loci()/build_locus_blocks() for the per-locus pairing/CIGAR logic and
# em.R for the EM step.
.count_transcript_abundance <- function(bam_path, gene_idx, te_idx, stranded, mode,
                                         num_itr, frag_length, max_l = 500L) {
  aln <- .load_alignments(bam_path)
  n <- length(aln)
  qname <- S4Vectors::mcols(aln)$qname
  flag  <- S4Vectors::mcols(aln)$flag
  is_paired <- bitwAnd(flag[1], .FLAG_PAIRED) != 0

  read_key <- if (is_paired) .strip_mate_suffix(qname) else qname
  is_read1 <- bitwAnd(flag, .FLAG_READ1) != 0
  is_read2 <- bitwAnd(flag, .FLAG_READ2) != 0
  aln_strand <- as.character(GenomicAlignments::strand(aln))

  at <- data.table::data.table(row = seq_len(n), rkey = read_key, is_read1 = is_read1,
                                is_read2 = is_read2, aln_strand = aln_strand, orig = seq_len(n))
  at[, group_row := seq_len(.N), by = rkey]  # stable within-group order (file order)

  if (is_paired) {
    at[, `:=`(n1 = sum(is_read1), n2 = sum(is_read2)), by = rkey]
    at[, is_unique_group := n1 <= 1L & n2 <= 1L]
  } else {
    at[, n1 := .N, by = rkey]
    at[, is_unique_group := n1 == 1L]
  }

  loci <- .build_loci(at, is_paired, mode)

  n_gene <- length(gene_idx$genes)
  n_te   <- te_idx$n
  gene_counts <- numeric(n_gene)
  te_uniq_counts <- numeric(n_te)
  te_multi_counts <- numeric(n_te)  # naive (pre-EM) prior, as in TEtranscripts

  if (nrow(loci) == 0) {
    return(list(gene_counts = gene_counts, te_counts = te_uniq_counts,
                uniq_reads = 0L, nonunique = 0L, empty = 0L, avg_read_length = 0))
  }

  loci[, locus_id := .I]

  blocks <- .build_locus_blocks(loci, aln, stranded)
  blk_gr <- GenomicRanges::GRanges(blocks$seqnames, IRanges::IRanges(blocks$start, blocks$end))
  gene_hits <- .overlap_genes(blk_gr, blocks$locus_id, blocks$dir, gene_idx)
  te_hits   <- .overlap_tes(blk_gr, blocks$locus_id, blocks$dir, te_idx)

  # per-read (not per-locus) alignment count, to split unique vs multi reads,
  # exactly mirroring Counting.py's len(alignments_per_read) check.
  read_nloci <- loci[, .N, by = rkey]
  data.table::setnames(read_nloci, "N", "nloci")
  loci <- merge(loci, read_nloci, by = "rkey", all.x = TRUE, sort = FALSE)

  uniq_read_loci  <- loci[nloci == 1L]
  multi_read_loci <- loci[nloci > 1L]

  ## ---- unique reads: gene first, TE second (Counting.py's uniq-read path) ----
  ug <- gene_hits[locus_id %in% uniq_read_loci$locus_id]
  ut <- te_hits[locus_id %in% uniq_read_loci$locus_id]

  gene_n_per_locus <- ug[, .N, by = locus_id]
  gene_ambig_loci  <- gene_n_per_locus[N > 1, locus_id]
  gene_simple_loci <- gene_n_per_locus[N == 1, locus_id]
  if (length(gene_simple_loci) > 0) {
    simple_gene <- ug[locus_id %in% gene_simple_loci]
    tab <- simple_gene[, .N, by = gene_idx]
    gene_counts[tab$gene_idx] <- gene_counts[tab$gene_idx] + tab$N
  }
  # A unique read's single locus overlapping multiple genes (e.g. two
  # transcript isoforms' exons) is split evenly and IMMEDIATELY -- mirrors
  # Counting.py::parse_annotations_gene's `len(annot_gene) == 1` branch. This
  # is *not* deferred to resolve_annotation_ambiguity the way TE's equivalent
  # case is; getting this wrong makes a confidently-assigned neighboring gene
  # silently absorb every ambiguous read (see the TEfastR vignette/tests for
  # the real-data case that exposed this).
  if (length(gene_ambig_loci) > 0) {
    ambig_gene <- ug[locus_id %in% gene_ambig_loci]
    ambig_gene[, w := 1 / .N, by = locus_id]
    tab <- ambig_gene[, .(w = sum(w)), by = gene_idx]
    gene_counts[tab$gene_idx] <- gene_counts[tab$gene_idx] + tab$w
  }

  no_gene_loci <- setdiff(uniq_read_loci$locus_id, ug$locus_id)  # unique reads with 0 gene hits
  ut_for_no_gene <- ut[locus_id %in% no_gene_loci]
  te_n_per_locus <- ut_for_no_gene[, .N, by = locus_id]
  te_simple_loci <- te_n_per_locus[N == 1, locus_id]
  te_ambig_loci  <- te_n_per_locus[N > 1, locus_id]
  if (length(te_simple_loci) > 0) {
    simple_te <- ut_for_no_gene[locus_id %in% te_simple_loci]
    tab <- simple_te[, .N, by = te_idx]
    te_uniq_counts[tab$te_idx] <- te_uniq_counts[tab$te_idx] + tab$N
  }
  te_leftover_ids <- list(); te_leftover_w <- numeric(0)
  if (length(te_ambig_loci) > 0) {
    grp <- split(ut_for_no_gene[locus_id %in% te_ambig_loci]$te_idx,
                 ut_for_no_gene[locus_id %in% te_ambig_loci]$locus_id)
    te_leftover_ids <- c(te_leftover_ids, grp)
    te_leftover_w <- c(te_leftover_w, rep(1, length(grp)))
  }
  empty_uniq <- length(setdiff(no_gene_loci, ut_for_no_gene$locus_id))

  gene_leftover_ids <- list()
  gene_leftover_w <- numeric(0)

  ## ---- multi-locus reads: TE first, gene second ----
  mt <- te_hits[locus_id %in% multi_read_loci$locus_id]
  mg <- gene_hits[locus_id %in% multi_read_loci$locus_id]

  # Build per-read candidate TE lists across ALL of that read's TE-hit loci
  # (flat, with repeats across loci) -- this is the multi_algn list from
  # parse_annotations_TE() in TEtranscripts, needed both for the naive prior
  # and for EM's incidence matrix.
  #
  # Important: TEtranscripts' ovp_annotation() only appends a locus to
  # annot_TE when that locus overlaps >=1 TE (`if len(TEs) > 0: annot_TE.
  # append(TEs)`); loci with zero TE hits are skipped entirely, not
  # represented by an empty entry. So `len(annot_TE)` -- the divisor in the
  # naive prior formula -- is the count of *TE-hit* loci for a read, not its
  # total locus count (`nloci`, which also includes loci with no TE hit at
  # all, e.g. a genomic region with no repeat annotation). Using `nloci`
  # there systematically underweights every candidate for reads that have
  # any non-TE-hit loci.
  read_of_locus <- unique(multi_read_loci[, .(locus_id, rkey)])
  mt2 <- merge(mt, read_of_locus, by = "locus_id", all.x = TRUE, sort = FALSE)
  n_cand_per_locus <- mt2[, .N, by = locus_id]
  data.table::setnames(n_cand_per_locus, "N", "n_cand")
  n_te_loci_per_read <- unique(mt2[, .(locus_id, rkey)])[, .N, by = rkey]
  data.table::setnames(n_te_loci_per_read, "N", "n_te_loci")
  mt2 <- merge(mt2, n_cand_per_locus, by = "locus_id", all.x = TRUE, sort = FALSE)
  mt2 <- merge(mt2, n_te_loci_per_read, by = "rkey", all.x = TRUE, sort = FALSE)

  reads_with_any_te <- unique(mt2$rkey)
  n_reads_multi <- length(reads_with_any_te)

  if (n_reads_multi > 0) {
    read_id_map <- stats::setNames(seq_along(reads_with_any_te), reads_with_any_te)
    mt2[, read_id := read_id_map[rkey]]

    # naive prior: 1 / (n_te_loci_for_read * n_cand_at_that_locus), summed per TE
    mt2[, w := 1 / (n_te_loci * n_cand)]
    prior_tab <- mt2[, .(w = sum(w)), by = te_idx]
    te_multi_counts[prior_tab$te_idx] <- te_multi_counts[prior_tab$te_idx] + prior_tab$w

    # EM incidence matrix: entries counted with multiplicity across loci
    inc_tab <- mt2[, .(mult = .N), by = .(read_id, te_idx)]
    M <- Matrix::sparseMatrix(i = inc_tab$read_id, j = inc_tab$te_idx, x = inc_tab$mult,
                               dims = c(n_reads_multi, n_te))
  } else {
    M <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(0, n_te))
  }

  # Reads with genuinely 0 TE hits across all loci fall back to gene
  # assignment. Mirrors Counting.py's parse_annotations_gene `len(annot_gene)
  # > 1` branch: ALWAYS deferred to resolve_annotation_ambiguity, as ONE
  # leftover entry per READ (not per locus) that pools the distinct genes
  # hit across every one of that read's loci (loci with 0 or exactly 1 gene
  # hit are pooled in too, not counted immediately). The weight divisor is
  # the count of *gene-hit* loci for the read (mirroring ovp_annotation's
  # `if len(genes) > 0: annot_gene.append(...)` filtering -- see the same
  # note above the TE naive-prior calculation), not its total locus count.
  reads_no_te <- setdiff(unique(multi_read_loci$rkey), reads_with_any_te)
  if (length(reads_no_te) > 0) {
    loci_no_te <- multi_read_loci[rkey %in% reads_no_te]
    mg_no_te <- mg[locus_id %in% loci_no_te$locus_id]
    read_of_locus_no_te <- unique(loci_no_te[, .(locus_id, rkey)])
    mg_no_te <- merge(mg_no_te, read_of_locus_no_te, by = "locus_id", all.x = TRUE, sort = FALSE)
    per_read_genes <- unique(mg_no_te[, .(rkey, gene_idx)])
    gene_grp <- split(per_read_genes$gene_idx, per_read_genes$rkey)

    n_gene_loci_per_read <- unique(mg_no_te[, .(locus_id, rkey)])[, .N, by = rkey]
    w_per_read <- stats::setNames(ifelse(n_gene_loci_per_read$N > 1, 1 / n_gene_loci_per_read$N, 1),
                                   n_gene_loci_per_read$rkey)

    gene_leftover_ids <- c(gene_leftover_ids, gene_grp)
    gene_leftover_w <- c(gene_leftover_w, unname(w_per_read[names(gene_grp)]))
    empty_multi <- length(reads_no_te) - length(gene_grp)
  } else {
    empty_multi <- 0L
  }

  ## ---- resolve ambiguity, then EM ----
  if (length(gene_leftover_ids) > 0)
    gene_counts <- .resolve_ambiguity(gene_counts, gene_leftover_ids, gene_leftover_w)
  if (length(te_leftover_ids) > 0)
    te_uniq_counts <- .resolve_ambiguity(te_uniq_counts, te_leftover_ids, te_leftover_w)

  avg_read_length <- .estimate_read_length(aln, uniq_read_loci, is_paired, frag_length, max_l)

  if (num_itr > 0 && n_reads_multi > 0) {
    # EMAlgorithm.py::EMestimate() is called with an all-zero uniq_counts
    # placeholder (`te_tmp_counts`), not the real per-TE unique counts --
    # the EM's internal means/effective-length calculation never sees the
    # real unique counts, which are only added to its output afterward
    # (Counting.py's `te_counts = te_counts + new_te_multi_counts`, where
    # `te_counts` already holds the real unique counts at that point).
    # Passing the real counts in here instead skews the EM's starting point
    # so TEs with many confident unique reads absorb ambiguous multi-mapped
    # reads that belong elsewhere.
    new_multi <- .em_estimate(numeric(n_te), te_multi_counts, M, te_idx$length,
                               num_itr, avg_read_length)
    te_counts <- te_uniq_counts + new_multi
  } else {
    te_counts <- te_uniq_counts + te_multi_counts
  }

  uniq_reads <- nrow(unique(uniq_read_loci[, .(rkey)]))
  nonunique <- length(unique(multi_read_loci$rkey))
  empty <- empty_uniq + empty_multi

  list(gene_counts = gene_counts, te_counts = te_counts,
       uniq_reads = uniq_reads, nonunique = nonunique, empty = empty,
       avg_read_length = avg_read_length)
}
