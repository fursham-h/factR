context("Test buildCDS function")

library(BSgenome.Mmusculus.UCSC.mm10)
out1 <- buildCDS(query_exons, ref_cds, Mmusculus, q2rcovs)
q2rcovserr <- q2rcovs
q2rcovserr$transcript_id <- c('1', '2', '3', '4')

test_that("Test proper input arguments", {
  expect_error(buildCDS(query_exons[[1]], ref_cds[[1]], Mmusculus, q2rcovs))
  expect_error(buildCDS(query_exons, ref_cds, Mmusculus, q2rcovs, ids = c(4,5)))
  expect_error(buildCDS(query_exons, ref_cds, Mmusculus, q2rcovs, ids = c(1,1)))
  expect_error(buildCDS(query_exons, ref_cds, Mmusculus, q2rcovserr))
})

test_that("Test structure of outputs", {
  expect_equal(is(out1, "GRangesList"), TRUE)
  expect_equal(length(out1), 4)
})

test_that("Test value of outputs", {
  expect_identical(ranges(out1), ranges(query_cds))
})
