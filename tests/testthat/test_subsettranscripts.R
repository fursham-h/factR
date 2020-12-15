context("Test subset function")

test_that("Test subset robustness", {
  expect_error(subsetNewTranscripts(query_gtf, ref_cds))
  expect_error(buildCDS(query_gtf))
})

test_that("Test .subsetTranscripts", {
  out <- .subsetTranscripts(matched_query_gtf, ref_gtf, F, F)
  expect_equal(is(out, "GRanges"), T)
  expect_equal(length(out), 27)
  expect_equal(unique(out$transcript_id), c("transcript3", "transcript4"))
})

test_that("Test .subsetTranscripts by.name", {
  out <- .subsetTranscripts(matched_query_gtf, matched_query_gtf, T, F)
  expect_equal(length(out), 0)
})

test_that("Test .subsetTranscripts by.CDS", {
  out <- .subsetTranscripts(new_query_gtf, ref_gtf, F, T)
  expect_equal(length(out), 49)
  expect_equal(unique(out$transcript_id), c("transcript3", "transcript4"))
})