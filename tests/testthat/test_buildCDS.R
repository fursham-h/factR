context("Test buildCDS function")

library(BSgenome.Mmusculus.UCSC.mm10)

test_that("Test buildCDS robustness", {
  expect_error(buildCDS(query_gtf, ref_cds, Mmusculus))
  expect_error(buildCDS(query_gtf, query_gtf, Mmusculus))
  expect_error(buildCDS(query_gtf, ref_gtf, Mmusculus))
  
  query_gtf_seqmatched <- matchSeqLevels(query_gtf, Mmusculus)
  expect_error(buildCDS(query_gtf_seqmatched, ref_gtf, Mmusculus))
  
  
  expect_error(buildCDS(query_exons, ref_cds, Mmusculus, q2rcovs, ids = c(4,5)))
  expect_error(buildCDS(query_exons, ref_cds, Mmusculus, q2rcovs, ids = c(1,1)))
  expect_error(buildCDS(query_exons, ref_cds, Mmusculus, q2rcovserr))
})

test_that("Test .getCDSranges", {
  out <- .getCDSranges(query_exons[[1]], 900, 2000, 'test')
  expect_equal(out$phase, c(0,1))
})

test_that("Test .getCDSstart", {
  out <- .getCDSstart(query_exons[[1]], ref_cds[[1]], Mmusculus)
  expect_equal(out$ORF_start, "Annotated ATG")
  expect_equal(out$fiveUTRlength, 287)
  
  out <- .getCDSstart(query_exons[[1]][3:14], ref_cds[[1]], Mmusculus)
  expect_equal(out$ORF_start, "Internal ATG")
  expect_equal(out$fiveUTRlength, 51)
  
  out <- .getCDSstart(query_exons[[1]][6], ref_cds[[1]], Mmusculus)
  expect_equal(out$ORF_start, "Inferred frame")
  expect_equal(out$fiveUTRlength, 0)
  
  out <- .getCDSstart(query_exons[[1]][8], ref_cds[[1]], Mmusculus)
  expect_equal(out$ORF_start, "Not found")
  expect_equal(out$fiveUTRlength, 0)
})

test_that("Test .getCDSstop", {
  out <- .getCDSstop(query_exons[[1]], Mmusculus, fiveUTRlength = 287)
  expect_equal(out$ORF_found, T)
  expect_equal(out$threeUTRlength, 1155)
  
  out <- .getCDSstop(query_exons[[1]], Mmusculus, fiveUTRlength = 3100)
  expect_equal(out$ORF_found, F)
  expect_equal(out$threeUTRlength, 0)
})

test_that("Test .getCDSgr", {
  out <- .getCDSgr(query_exons, ref_cds, Mmusculus, q2rcovs)
  expect_equal(is(out, "GRanges"), T)
  expect_equal(length(out), 49)
  expect_equal(BiocGenerics::start(out), BiocGenerics::start(unlist(query_cds)))
})

test_that("Test .calcCoverage", {
  expect_equal(signif(.calcCoverage(query_exons[[1]], query_exons[[2]], 'mean'),3), 0.901)
  expect_equal(signif(.calcCoverage(query_exons[[1]], query_exons[[2]], 'query'),3), 0.91)
  expect_equal(signif(.calcCoverage(query_exons[[1]], query_exons[[2]], 'ref'),3), 0.893)
})
