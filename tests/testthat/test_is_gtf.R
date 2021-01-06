context("Test is_gtf function")

test_that("Test the value of outputs", {
  expect_equal(is_gtf(query_gtf), TRUE)
  expect_equal(is_gtf(query_exons[1]), FALSE)
  expect_equal(is_gtf(query_exons[[1]]), FALSE)
  expect_equal(is_gtf(data.frame("name" = "Me")), FALSE)
})
