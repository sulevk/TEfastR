# Regression tests cross-checked, during development, against TEtranscripts'
# own Python counting core (TEToolkit/Counting.py) on identical synthetic
# inputs. Covers: unique gene/TE assignment, single-locus ambiguity (gene and
# TE), a genuine multi-mapped read (exercises the EM step), paired-end
# pairing/singletons, and stranded (no/forward/reverse) matching.

make_bam <- function(sam_body, tmp, name) {
  sam_header <- "@HD\tVN:1.4\n@SQ\tSN:chr1\tLN:50000\n"
  sam_path <- file.path(tmp, paste0(name, ".sam"))
  writeLines(paste0(sam_header, sam_body), sam_path, sep = "")
  Rsamtools::asBam(sam_path, file.path(tmp, name), overwrite = TRUE, indexDestination = TRUE)
}

write_gtf <- function(lines, tmp, name) {
  path <- file.path(tmp, name)
  writeLines(lines, path)
  path
}

expect_counts <- function(result, expected) {
  df <- as.data.frame(result)
  got <- stats::setNames(df$count, df$id)
  for (nm in names(expected)) {
    expect_true(nm %in% names(got), info = paste("missing feature", nm))
    expect_equal(got[[nm]], expected[[nm]], tolerance = 1e-6, label = nm)
  }
}

test_that("single-end: unique gene/TE, ambiguity, multi-mapped (EM)", {
  te_gtf <- system.file("extdata", "te.gtf", package = "TEfastR")
  gene_gtf <- system.file("extdata", "gene.gtf", package = "TEfastR")
  bam <- system.file("extdata", "test.bam", package = "TEfastR")

  te_index <- build_te_annotation(te_gtf)
  gene_index <- build_gene_annotation(gene_gtf)
  result <- count_te(bam, gene_index, te_index, stranded = "no", mode = "multi", iteration = 100L)

  expect_counts(result, list(GeneA = 1, GeneC = 0.5, GeneD = 0.5,
                             "TE1:FAM1:CLASS1" = 1, "TE2:FAM2:CLASS2" = 0.5,
                             "TE3:FAM2:CLASS2" = 0.5, "TE4:FAM3:CLASS3" = 0.5,
                             "TE5:FAM3:CLASS3" = 0.5))
  expect_equal(result$stats$uniq_reads, 4L)
  expect_equal(result$stats$nonunique, 1L)
  expect_equal(result$stats$empty, 0L)
})

test_that("paired-end: unique pair, singleton, multi-mapped pair (EM)", {
  tmp <- withr::local_tempdir()
  gene_gtf <- write_gtf('chr1\tt\texon\t100\t300\t.\t+\t.\tgene_id "GeneA";', tmp, "gene2.gtf")
  te_gtf <- write_gtf(c(
    'chr1\tt\texon\t1000\t1500\t.\t+\t.\tgene_id "TEA"; transcript_id "TEA_i1"; family_id "F1"; class_id "C1";',
    'chr1\tt\texon\t2000\t2500\t.\t+\t.\tgene_id "TEB"; transcript_id "TEB_i1"; family_id "F1"; class_id "C1";',
    'chr1\tt\texon\t3000\t3500\t.\t+\t.\tgene_id "TEC"; transcript_id "TEC_i1"; family_id "F2"; class_id "C2";'
  ), tmp, "te2.gtf")
  sam <- paste0(
    "readP_unique\t99\tchr1\t120\t60\t30M\t=\t220\t130\t", strrep("A", 30), "\t", strrep("I", 30), "\n",
    "readP_unique\t147\tchr1\t220\t60\t30M\t=\t120\t-130\t", strrep("A", 30), "\t", strrep("I", 30), "\n",
    "readP_multi\t99\tchr1\t1020\t60\t30M\t=\t1050\t60\t", strrep("A", 30), "\t", strrep("I", 30), "\n",
    "readP_multi\t147\tchr1\t1050\t60\t30M\t=\t1020\t-60\t", strrep("A", 30), "\t", strrep("I", 30), "\n",
    "readP_multi\t355\tchr1\t2020\t0\t30M\t=\t2050\t60\t", strrep("A", 30), "\t", strrep("I", 30), "\n",
    "readP_multi\t403\tchr1\t2050\t0\t30M\t=\t2020\t-60\t", strrep("A", 30), "\t", strrep("I", 30), "\n",
    "readP_singleton\t65\tchr1\t3020\t60\t30M\t*\t0\t0\t", strrep("A", 30), "\t", strrep("I", 30), "\n"
  )
  bam <- make_bam(sam, tmp, "test2")

  te_index <- build_te_annotation(te_gtf)
  gene_index <- build_gene_annotation(gene_gtf)
  result <- count_te(bam, gene_index, te_index, stranded = "no", mode = "multi", iteration = 100L)

  expect_counts(result, list(GeneA = 1, "TEA:F1:C1" = 0.5, "TEB:F1:C1" = 0.5, "TEC:F2:C2" = 1))
})

test_that("stranded no/forward/reverse match TE strand as expected", {
  tmp <- withr::local_tempdir()
  gene_gtf <- write_gtf('chr1\tt\texon\t1\t2\t.\t+\t.\tgene_id "Dummy";', tmp, "gene3.gtf")
  te_gtf <- write_gtf(
    'chr1\tt\texon\t4000\t4400\t.\t-\t.\tgene_id "TEminus"; transcript_id "TEminus_i1"; family_id "F"; class_id "C";',
    tmp, "te3.gtf")
  sam <- paste0("read_fwd\t0\tchr1\t4050\t60\t50M\t*\t0\t0\t", strrep("A", 50), "\t", strrep("I", 50), "\n")
  bam <- make_bam(sam, tmp, "test3")

  te_index <- build_te_annotation(te_gtf)
  gene_index <- build_gene_annotation(gene_gtf)

  expect_counts(count_te(bam, gene_index, te_index, stranded = "no", iteration = 0L),
                list("TEminus:F:C" = 1))
  expect_counts(count_te(bam, gene_index, te_index, stranded = "forward", iteration = 0L),
                list("TEminus:F:C" = 0))
  expect_counts(count_te(bam, gene_index, te_index, stranded = "reverse", iteration = 0L),
                list("TEminus:F:C" = 1))
})

test_that("single-locus multi-gene ambiguity is an immediate uniform split, not deferred", {
  # Regression test for a real bug found on real GENCODE data: a gene with an
  # exclusive exon (giving it a nonzero "confident" baseline count) plus an
  # exon shared with a second, otherwise-exclusive-free gene (e.g. a small
  # ncRNA nested inside a multi-isoform host gene's shared terminal exon).
  # TEtranscripts splits every such ambiguous read evenly and IMMEDIATELY
  # (Counting.py::parse_annotations_gene's `len(annot_gene) == 1` branch) --
  # it does NOT defer to resolve_annotation_ambiguity the way TE's equivalent
  # case does. Deferring here (this package's original bug) makes the
  # confidently-assigned gene silently absorb every ambiguous read instead of
  # splitting 50/50, because the other gene never has a nonzero baseline to
  # split proportionally against.
  tmp <- withr::local_tempdir()
  gene_gtf <- write_gtf(c(
    'chr1\tt\texon\t100\t500\t.\t+\t.\tgene_id "GeneA";',   # GeneA's own exclusive exon
    'chr1\tt\texon\t1000\t1100\t.\t+\t.\tgene_id "GeneA";', # GeneA's shared terminal exon
    'chr1\tt\texon\t1000\t1100\t.\t+\t.\tgene_id "GeneB";'  # GeneB: only this shared exon
  ), tmp, "gene4.gtf")
  te_gtf <- write_gtf(
    'chr1\tt\texon\t9000\t9100\t.\t+\t.\tgene_id "TEx"; transcript_id "TEx_i1"; family_id "F"; class_id "C";',
    tmp, "te4.gtf")
  sam <- paste0(
    paste0("excl", 1:4, "\t0\tchr1\t", c(150, 200, 250, 300), "\t60\t50M\t*\t0\t0\t",
           strrep("A", 50), "\t", strrep("I", 50), collapse = "\n"), "\n",
    paste0("shared", 1:4, "\t0\tchr1\t", c(1020, 1030, 1040, 1050), "\t60\t50M\t*\t0\t0\t",
           strrep("A", 50), "\t", strrep("I", 50), collapse = "\n"), "\n"
  )
  bam <- make_bam(sam, tmp, "test4")

  te_index <- build_te_annotation(te_gtf)
  gene_index <- build_gene_annotation(gene_gtf)
  result <- count_te(bam, gene_index, te_index, stranded = "no", iteration = 0L)

  # 4 exclusive reads + 4 ambiguous reads split 50/50 => GeneA = 4 + 2 = 6, GeneB = 2
  expect_counts(result, list(GeneA = 6, GeneB = 2))
})

test_that("fragment-length estimate matches TEtranscripts: floored, maxL-filtered, unique-only", {
  # Regression test for a real bug found on real data: Counting.py estimates
  # the average fragment length from UNIQUE reads only, drops proper pairs
  # with an insert size over `maxL` (default 500) before counting them
  # toward its first-10000 cap, and floor-divides the final average
  # (`avgReadLength // tmp_cnt`, Python integer division). That estimate
  # feeds TE effective length = TE length - fragment length + 1, gating
  # whether a TE is considered at all (effLen <= 0 => excluded). Getting any
  # of the three wrong can flip a short TE from included to excluded (or
  # vice versa) and hand its multi-mapped reads to the wrong element
  # entirely -- exactly what a single un-filtered outlier insert size does
  # here: it would inflate the fragment-length estimate enough to exclude
  # ShortTE and give the ambiguous read entirely to LongTE instead.
  tmp <- withr::local_tempdir()
  gene_gtf <- write_gtf('chr1\tt\texon\t1\t2\t.\t+\t.\tgene_id "Dummy";', tmp, "gene5.gtf")
  te_gtf <- write_gtf(c(
    'chr1\tt\texon\t20000\t20150\t.\t+\t.\tgene_id "ShortTE"; transcript_id "ShortTE_i1"; family_id "F"; class_id "C";',
    'chr1\tt\texon\t30000\t32000\t.\t+\t.\tgene_id "LongTE"; transcript_id "LongTE_i1"; family_id "F"; class_id "C";'
  ), tmp, "te5.gtf")

  # 9 unique pairs with a normal ~100bp fragment, establishing the estimate
  p1 <- 1000 + (0:8) * 300; p2 <- p1 + 70
  norm <- paste0(
    "norm", 0:8, "\t99\tchr1\t", p1 + 1, "\t60\t30M\t=\t", p2 + 1, "\t100\t",
    strrep("A", 30), "\t", strrep("I", 30), "\n",
    "norm", 0:8, "\t147\tchr1\t", p2 + 1, "\t60\t30M\t=\t", p1 + 1, "\t-100\t",
    strrep("A", 30), "\t", strrep("I", 30), collapse = "\n"
  )
  # 1 unique outlier pair, insert size 9000 -- must be excluded by maxL
  op1 <- 1000 + 9 * 300; op2 <- op1 + 9000
  outlier <- paste0(
    "outlier\t99\tchr1\t", op1 + 1, "\t60\t30M\t=\t", op2 + 1, "\t9030\t",
    strrep("A", 30), "\t", strrep("I", 30), "\n",
    "outlier\t147\tchr1\t", op2 + 1, "\t60\t30M\t=\t", op1 + 1, "\t-9030\t",
    strrep("A", 30), "\t", strrep("I", 30)
  )
  # 1 multi-mapped pair: locus1 -> ShortTE, locus2 -> LongTE
  multi <- paste0(
    "multi\t99\tchr1\t20031\t60\t30M\t=\t20061\t60\t", strrep("A", 30), "\t", strrep("I", 30), "\n",
    "multi\t147\tchr1\t20061\t60\t30M\t=\t20031\t-60\t", strrep("A", 30), "\t", strrep("I", 30), "\n",
    "multi\t355\tchr1\t30031\t0\t30M\t=\t30061\t60\t", strrep("A", 30), "\t", strrep("I", 30), "\n",
    "multi\t403\tchr1\t30061\t0\t30M\t=\t30031\t-60\t", strrep("A", 30), "\t", strrep("I", 30)
  )
  sam <- paste0(norm, "\n", outlier, "\n", multi, "\n")
  bam <- make_bam(sam, tmp, "test5")

  te_index <- build_te_annotation(te_gtf)
  gene_index <- build_gene_annotation(gene_gtf)
  result <- count_te(bam, gene_index, te_index, stranded = "no", iteration = 100L)

  expect_equal(result$stats$avg_read_length, 100)
  # matches TEtranscripts' Python counting core on this exact fixture
  expect_equal(unname(result$te_element_counts["ShortTE:F:C"]), 0.9999795652363824, tolerance = 1e-6)
  expect_equal(unname(result$te_element_counts["LongTE:F:C"]), 2.043476e-05, tolerance = 1e-6)
})

test_that("count_te_batch assembles a feature x sample matrix", {
  te_gtf <- system.file("extdata", "te.gtf", package = "TEfastR")
  gene_gtf <- system.file("extdata", "gene.gtf", package = "TEfastR")
  bam <- system.file("extdata", "test.bam", package = "TEfastR")
  te_index <- build_te_annotation(te_gtf)
  gene_index <- build_gene_annotation(gene_gtf)

  mat <- count_te_batch(c(s1 = bam, s2 = bam), gene_index, te_index, iteration = 100L)
  expect_equal(dim(mat), c(length(gene_index$genes) + length(te_index$elements), 2L))
  expect_equal(colnames(mat), c("s1", "s2"))
  expect_equal(mat["GeneA", "s1"], mat["GeneA", "s2"])
})
