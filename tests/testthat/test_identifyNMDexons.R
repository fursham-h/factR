context("Test identifyNMDexons function")


test_that(".grouptxbycdsstart functionality", {
  # output from new_query_gtf as input
  out <- .grouptxbycdsstart(new_query_gtf)
  expect_equal(nrow(out), 4)
  expect_equal(out[1,1], "chr10_79854714_+")
  
  # test reverse strand
  test.gtf <- new_query_gtf
  strand(test.gtf) <-  "-"
  out <- .grouptxbycdsstart(test.gtf)
  expect_equal(nrow(out), 4)
  expect_equal(unique(out$coord), c("chr10_79863274_-", "chr10_79862472_-"))
})
