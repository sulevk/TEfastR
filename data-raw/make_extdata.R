# Generates inst/extdata/{gene.gtf, te.gtf, test.bam(.bai)}: the small
# synthetic single-end dataset used by roxygen @examples, the vignette, and
# tests/testthat. Covers unique gene/TE assignment, single-locus ambiguity
# (gene and TE), and a genuine multi-mapped read (exercises the EM step).
# Re-run with: Rscript data-raw/make_extdata.R
suppressPackageStartupMessages(library(Rsamtools))

out_dir <- file.path("inst", "extdata")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

writeLines(c(
  'chr1\tt\texon\t100\t200\t.\t+\t.\tgene_id "GeneA";',
  'chr1\tt\texon\t800\t900\t.\t+\t.\tgene_id "GeneC";',
  'chr1\tt\texon\t850\t950\t.\t+\t.\tgene_id "GeneD";'
), file.path(out_dir, "gene.gtf"))

writeLines(c(
  'chr1\tt\texon\t500\t600\t.\t+\t.\tgene_id "TE1"; transcript_id "TE1_i1"; family_id "FAM1"; class_id "CLASS1";',
  'chr1\tt\texon\t1000\t1100\t.\t+\t.\tgene_id "TE2"; transcript_id "TE2_i1"; family_id "FAM2"; class_id "CLASS2";',
  'chr1\tt\texon\t1050\t1150\t.\t+\t.\tgene_id "TE3"; transcript_id "TE3_i1"; family_id "FAM2"; class_id "CLASS2";',
  'chr1\tt\texon\t1200\t1300\t.\t+\t.\tgene_id "TE4"; transcript_id "TE4_i1"; family_id "FAM3"; class_id "CLASS3";',
  'chr1\tt\texon\t1400\t1500\t.\t+\t.\tgene_id "TE5"; transcript_id "TE5_i1"; family_id "FAM3"; class_id "CLASS3";'
), file.path(out_dir, "te.gtf"))

sam_header <- "@HD\tVN:1.4\n@SQ\tSN:chr1\tLN:10000\n"
sam_body <- paste0(
  "read_gene_unique\t0\tchr1\t150\t60\t50M\t*\t0\t0\t", strrep("A", 50), "\t", strrep("I", 50), "\n",
  "read_te_unique\t0\tchr1\t530\t60\t50M\t*\t0\t0\t", strrep("A", 50), "\t", strrep("I", 50), "\n",
  "read_gene_ambig\t0\tchr1\t870\t60\t20M\t*\t0\t0\t", strrep("A", 20), "\t", strrep("I", 20), "\n",
  "read_te_ambig\t0\tchr1\t1070\t60\t20M\t*\t0\t0\t", strrep("A", 20), "\t", strrep("I", 20), "\n",
  "read_multi\t0\tchr1\t1230\t60\t20M\t*\t0\t0\t", strrep("A", 20), "\t", strrep("I", 20), "\n",
  "read_multi\t256\tchr1\t1430\t0\t20M\t*\t0\t0\t", strrep("A", 20), "\t", strrep("I", 20), "\n"
)
sam_path <- file.path(out_dir, "test.sam")
writeLines(paste0(sam_header, sam_body), sam_path, sep = "")
asBam(sam_path, file.path(out_dir, "test"), overwrite = TRUE, indexDestination = TRUE)
file.remove(sam_path)

message("Wrote ", out_dir)
