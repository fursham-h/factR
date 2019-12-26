context("Test buildCDS function")

out1 = buildCDS(query_exons, ref_cds, Mmusculus, q2rcovs)

test_that("Test the structure of outputs", {
  
  expect_equal(is(out1,'GRangesList'), TRUE)
  expect_equal(length(out1), 4)
  
})

test_that("Test the value of outputs", {
  
  expect_equal(ranges(out1), ranges(query_cds)) 
 
})
