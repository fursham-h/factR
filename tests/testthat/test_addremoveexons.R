context("Test add/remove exon function")

# setup GRanges
gr1 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 2),
  ranges = IRanges::IRanges(
    start = c(1, 501),
    end = c(100, 600)
  )
)

gr2 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 3),
  ranges = IRanges::IRanges(
    start = c(1, 201, 501),
    end = c(100, 300, 600)
  )
)

gr3 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 2),
  ranges = IRanges::IRanges(
    start = c(1,  401),
    end = c(100, 600)
  )
)

gr4 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 1),
  ranges = IRanges::IRanges(
    start = c(1),
    end = c(600)
  )
)

exon1 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 1),
  ranges = IRanges::IRanges(
    start = c(201),
    end = c(300)
  )
)

exon2 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 1),
  ranges = IRanges::IRanges(
    start = c(401),
    end = c(500)
  )
)

exon3 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 1),
  ranges = IRanges::IRanges(
    start = c(101),
    end = c(500)
  )
)

exonset1 <- c(exon1,exon2,exon3)
exonset2 <- exonset1 %>% as.data.frame() %>% 
  dplyr::mutate(seqnames = "2") %>% 
  makeGRangesFromDataFrame()
exonset3 <- exonset1 %>% as.data.frame() %>% 
  dplyr::mutate(strand = "-") %>% 
  makeGRangesFromDataFrame()

grl1 <- GenomicRanges::GRangesList(gr1,gr1,gr1)
names(grl1) <- c("test1", "test2", "tes3")

grl2 <- GenomicRanges::GRangesList(gr2,gr3,gr4)
names(grl2) <- c("test4", "test5", "test6")


test_that("Test addExonstoTx function", {
 out <- addExonstoTx(grl1, exonset1)
 expect_equal(out[[1]], gr2)
 expect_equal(out[[2]], gr3)
 expect_equal(out[[3]], gr4)
 
 out <- addExonstoTx(grl1, exonset2)

})

test_that("Test removeExonsfromTx function", {
  out <- removeExonsfromTx(grl2, exonset1)
  expect_equal(out[[1]], gr1)
  expect_equal(out[[2]], gr1)
  expect_equal(out[[3]], gr1)
})

