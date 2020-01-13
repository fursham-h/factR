context("Test predictNMD functions")

test_that("Test .testNMD", {
  out <- .testNMD(query_exons[[1]], query_cds[[1]])
  expect_equal(out$is_NMD, F)
  expect_equal(out$dist_to_lastEJ, -133)
  expect_equal(out$threeUTRlength, 1155)
  
  out <- .testNMD(query_exons[[3]], query_cds[[3]])
  expect_equal(out$is_NMD, T)
  expect_equal(out$dist_to_lastEJ, 361)
  expect_equal(out$threeUTRlength, 641)
  expect_equal(out$dist_to_downEJs, c(66,283,361))
})

test_that("Test proper input arguments", {
  expect_error(predictNMD(query_exons[[1]], query_cds))
  expect_error(predictNMD(query_exons[[1]], q2rcovs))
  expect_error(predictNMD(q2r, q2rcovs))
})

test_that("Test the structure of outputs", {
  
  out1 <- predictNMD(query_exons[[1]], query_cds[[1]])
  out2 <- predictNMD(query_exons, query_cds, return = "all")
  
  
  expect_equal(names(out1), c(
    "is_NMD", "dist_to_lastEJ",
    "num_of_down_EJs", "dist_to_downEJs",
    "threeUTRlength"
  ))
  expect_equal(names(out2), c(
    "transcript", "is_NMD", "dist_to_lastEJ",
    "num_of_down_EJs", "dist_to_downEJs",
    "threeUTRlength"
  ))
  expect_equal(nrow(out2), 4)
})
