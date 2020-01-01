context("Test getCoverages function")

out1 <- calcCovs(query_exons, ref_exons, q2r)
out2 <- calcCovs(query_exons, ref_exons, q2r, return = "all")

test_that("Test the structure of outputs", {
  expect_equal(is(out1, "data.frame"), TRUE)
  expect_equal(is(out2, "data.frame"), TRUE)
  expect_equal(nrow(out1), 4)
  expect_equal(nrow(out2), 8)
})

test_that("Test the value of outputs", {
  expect_equal(out1, q2rcovs)
  expect_equal(round(out2$coverage[1:4], 3), c(1, 0.901, 0.901, 1))
})
