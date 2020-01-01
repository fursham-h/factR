context("Test predictNMD functions")

out1 <- predictNMD(query_exons[[1]], query_cds[[1]])
out2 <- predictNMD(query_exons, query_cds, return = "all")

test_that("Test the structure of outputs", {
  expect_equal(names(out1), c(
    "is_NMD", "dist_to_lastEJ",
    "num_of_down_EJs", "dist_to_downEJs",
    "threeUTRlength"
  ))
  expect_equal(names(out2), c(
    "exons", "is_NMD", "dist_to_lastEJ",
    "num_of_down_EJs", "dist_to_downEJs",
    "threeUTRlength"
  ))
  expect_equal(nrow(out2), 4)
})

test_that("Test the value of outputs", {
  expect_equal(out2$dist_to_lastEJ, c(-133, -133, 361, -133))
  expect_equal(out2$num_of_down_EJs, c(0, 0, 3, 0))
  expect_equal(out2$dist_to_downEJs, c("", "", "66,283,361", ""))
  expect_equal(out2$threeUTRlength, c(1155, 1494, 641, 1082))
})
