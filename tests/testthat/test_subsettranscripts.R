context("Test subset function")

test_that("Test subset robustness", {
  expect_error(subsetNewTranscripts(query_gtf, ref_cds))
  expect_error(buildCDS(query_gtf))
})

test_that("Test .subsetTranscripts", {
  out <- .subsetTranscripts(matched_query_gtf, ref_gtf)
  expect_equal(is(out, "GRanges"), T)
  expect_equal(length(out), 27)
  expect_equal(unique(out$transcript_id), c("transcript3", "transcript4"))
})
