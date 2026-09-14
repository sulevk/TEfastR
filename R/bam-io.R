# SAM flag bits (see the SAM spec): paired=0x1, read1=0x40, read2=0x80.
.FLAG_PAIRED <- 0x1L
.FLAG_READ1  <- 0x40L
.FLAG_READ2  <- 0x80L

#' @importFrom GenomicAlignments readGAlignments
.load_alignments <- function(bam_path) {
  flag <- Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE,
                                  isNotPassingQualityControls = FALSE,
                                  isSecondaryAlignment = NA)
  param <- Rsamtools::ScanBamParam(flag = flag, what = c("qname", "flag"))
  aln <- GenomicAlignments::readGAlignments(bam_path, param = param, use.names = FALSE)
  if (length(aln) == 0) stop("No alignments passed filtering in ", bam_path)
  aln
}

# TEtranscripts strips a trailing "/1"/"/2" (or, failing that, ".1"/".2")
# mate suffix from paired-end read names before grouping by name. Vectorized
# over the whole qname column at once.
.strip_mate_suffix <- function(qn) {
  slash <- regexpr("/", qn, fixed = TRUE)
  out <- qn
  has_slash <- slash > 0
  out[has_slash] <- substr(qn[has_slash], 1L, slash[has_slash] - 1L)
  dotted <- !has_slash & grepl("\\.[12]$", qn)
  out[dotted] <- substr(qn[dotted], 1L, nchar(qn[dotted]) - 2L)
  out
}
