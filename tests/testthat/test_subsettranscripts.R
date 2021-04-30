context("Test subset function")

data(query_gtf)
data(matched_query_gtf)
data(new_query_gtf)
data(ref_gtf)
data(ref_cds)

test_that("Test subset robustness", {
  expect_error(subsetNewTranscripts(query_gtf, ref_cds))
  expect_error(buildCDS(query_gtf))
})

test_that("Test .subsetTranscripts", {
  out <- .subsetTranscripts(matched_query_gtf, ref_gtf, "exon")
  expect_equal(is(out, "GRanges"), TRUE)
  expect_equal(length(out), 27)
  expect_equal(unique(out$transcript_id), c("transcript3", "transcript4"))
})

test_that("Test .subsetTranscripts by.name", {
  out <- .subsetTranscripts(matched_query_gtf, ref_gtf, "intron")
  expect_equal(length(out), 27)
  expect_equal(unique(out$transcript_id), c("transcript3", "transcript4"))
})

test_that("Test .subsetTranscripts by.CDS", {
  out <- .subsetTranscripts(new_query_gtf, ref_gtf, "cds")
  expect_equal(length(out), 49)
  expect_equal(unique(out$transcript_id), c("transcript3", "transcript4"))
})
