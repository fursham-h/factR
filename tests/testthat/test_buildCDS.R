context("Test buildCDS function")

library(BSgenome.Mmusculus.UCSC.mm10)

test_that("Test buildCDS robustness", {
  expect_error(buildCDS(query_gtf, ref_cds, Mmusculus))
  expect_error(buildCDS(query_gtf, query_gtf, Mmusculus))
  expect_error(buildCDS(query_gtf, ref_gtf, Mmusculus))
})

test_that("Test .runbuildCDS", {
  out <- .runbuildCDS(matched_query_gtf, ref_gtf, Mmusculus, "test")
  expect_equal(is(out, "GRanges"), T)
  expect_equal(length(out), 49)
  expect_equal(BiocGenerics::start(out), BiocGenerics::start(unlist(query_cds)))
})

# 
# test_that("Test .getCDSranges", {
#   out <- .getCDSranges(query_exons[[1]], 900, 2000, "test")
#   expect_equal(out$phase, c(0, 1))
# })
# 
# test_that("Test .getCDSstart", {
#   out <- .getCDSstart(query_exons[[1]], ref_cds[[1]], Mmusculus)
#   expect_equal(out$ORF_start, "Annotated ATG")
#   expect_equal(out$fiveUTRlength, 287)
# 
#   out <- .getCDSstart(query_exons[[1]][3:14], ref_cds[[1]], Mmusculus)
#   expect_equal(out$ORF_start, "Internal ATG")
#   expect_equal(out$fiveUTRlength, 51)
# 
#   out <- .getCDSstart(query_exons[[1]][6], ref_cds[[1]], Mmusculus)
#   expect_equal(out$ORF_start, "Inferred frame")
#   expect_equal(out$fiveUTRlength, 0)
# 
#   out <- .getCDSstart(query_exons[[1]][8], ref_cds[[1]], Mmusculus)
#   expect_equal(out$ORF_start, "Not found")
#   expect_equal(out$fiveUTRlength, 0)
# })
# 
# test_that("Test .getCDSstop", {
#   out <- .getCDSstop(query_exons[[1]], Mmusculus, fiveUTRlength = 287)
#   expect_equal(out$ORF_found, T)
#   expect_equal(out$threeUTRlength, 1155)
# 
#   out <- .getCDSstop(query_exons[[1]], Mmusculus, fiveUTRlength = 3100)
#   expect_equal(out$ORF_found, F)
#   expect_equal(out$threeUTRlength, 0)
# })
# 
# q2r <- .prepq2r(query_exons, ref_exons, ref_cds, query_exons[0])
# test_that("Test .prepq2r", {
#   expect_equal(q2r$ref_transcript_id, c(2, 1, 1, 1))
#   expect_equal(q2r$coverage, c(1, 1, 3, 3))
# })


