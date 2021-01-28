context("Test buildCDS function")

data(query_gtf)
data(query_cds)
data(matched_query_gtf)
data(new_query_gtf)
data(ref_gtf)
data(ref_cds)

library(BSgenome.Mmusculus.UCSC.mm10)

test_that("Test buildCDS robustness", {
  expect_error(buildCDS(query_gtf, ref_cds, Mmusculus))
  expect_error(buildCDS(query_gtf, query_gtf, Mmusculus))
  expect_error(buildCDS(query_gtf, ref_gtf, Mmusculus))
})

test_that("Test .runbuildCDS", {
  out <- .runbuildCDS(matched_query_gtf, ref_gtf, Mmusculus, "test")
  expect_equal(is(out, "GRanges"), TRUE)
  expect_equal(length(out), 49)
  expect_equal(BiocGenerics::start(out), BiocGenerics::start(unlist(query_cds)))
})

