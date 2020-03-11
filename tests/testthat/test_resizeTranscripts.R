context("Test trimTranscripts function")

# setup GRanges
gr1 <- GenomicRanges::GRanges(
  seqnames = "chr1", strand = c("+", "+", "+"),
  ranges = IRanges(
    start = c(1, 500, 1000),
    end = c(100, 600, 1100)
  )
)
gr2 <- GenomicRanges::GRanges(
  seqnames = "chr1", strand = c("-", "-", "-"),
  ranges = IRanges(
    start = c(1, 500, 1000),
    end = c(100, 600, 1100)
  )
)
grl <- GenomicRanges::GRangesList(gr1=gr1,gr2=gr2)

out1 <- trimTranscripts(gr1, 20, 80)
out2 <- trimTranscripts(gr2, 20, 80)
out3 <- trimTranscripts(gr1, 110, 150)
out4 <- trimTranscripts(gr2, 110, 150)

out5 <- trimTranscripts(grl, 20, 80)
out6 <- trimTranscripts(grl, c(20,110))


test_that("Test sample output", {
  expect_equal(BiocGenerics::start(out1), c(21, 500, 1000))
  expect_equal(BiocGenerics::end(out1), c(100, 600, 1020))

  expect_equal(BiocGenerics::start(out2), c(1000, 500, 81))
  expect_equal(BiocGenerics::end(out2), c(1080, 600, 100))

  expect_equal(BiocGenerics::start(out3), 510)
  expect_equal(BiocGenerics::end(out3), 551)

  expect_equal(BiocGenerics::start(out4), 550)
  expect_equal(BiocGenerics::end(out4), 591)
  
  expect_equal(as.numeric(BiocGenerics::unlist(BiocGenerics::start(out5))), 
               c(21, 500, 1000, 1000, 500, 81))
  expect_equal(as.numeric(BiocGenerics::unlist(BiocGenerics::end(out5))), 
               c(100, 600, 1020, 1080, 600, 100))
  
  expect_equal(as.numeric(BiocGenerics::unlist(BiocGenerics::start(out6))), 
               c(21, 500, 1000, 500, 1))
  expect_equal(as.numeric(BiocGenerics::unlist(BiocGenerics::end(out6))), 
               c(100, 600, 1100, 591, 100))
})
