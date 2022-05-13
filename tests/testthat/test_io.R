context("Test io function")

# setup GRanges
gtf <- system.file("extdata", "sc_merged_sample.gtf.gz", package = "factR")
gtf <- factR::importGTF(gtf)


test_that("Test import functionality", {
  expect_true(is_gtf(gtf))
  expect_equal(length(gtf), 8117)
  expect_error(factR::importGTF("fake.gtf"))
})
