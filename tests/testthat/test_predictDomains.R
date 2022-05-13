context("Test io function")

library(BSgenome.Mmusculus.UCSC.mm10)
data(new_query_gtf)


test_that("Test error catches", {
  expect_error(predictDomains(new_query_gtf))
  
})
