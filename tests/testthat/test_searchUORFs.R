context("Test searchuORFs functions")

out <- searchuORFs(query_exons, query_cds, Mmusculus)

test_that("Test the structure of outputs", {
  expect_equal(names(out), c("exons", "cds"))
  expect_equal(names(out$exons), c("uORF_1_transcript1"))
  expect_equal(sum(lengths(out$exons)), 14)
  expect_equal(sum(lengths(out$cds)), 1)
})

test_that("Test the value of outputs", {
  expect_equal(BiocGenerics::start(out$exons)[[1]][1:3], c(79854427, 79856504, 79858752))
  expect_equal(BiocGenerics::end(out$exons)[[1]][1:3], c(79854721, 79856534, 79858824))

  expect_equal(BiocGenerics::start(out$cds)[[1]], 79854489)
  expect_equal(BiocGenerics::end(out$cds)[[1]], 79854518)

  expect_equal(out$cds[[1]]$frame, 0)
  expect_equal(out$cds[[1]]$phase, 0)
})
