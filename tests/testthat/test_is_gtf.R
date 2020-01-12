context("Test is_gtf function")

test_that("Test the value of outputs", {
  expect_equal(is_gtf(query_gtf), T)
  expect_equal(is_gtf(query_exons[1]), F)
  expect_equal(is_gtf(query_exons[[1]]), F)
  expect_equal(is_gtf(q2rcovs), F)
})
