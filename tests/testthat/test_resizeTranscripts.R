context("Test resizeTranscript function")

# setup GRanges
gr1 <- GenomicRanges::GRanges(seqnames = "chr1", strand = c("+", "+", "+"),
              ranges = IRanges(start = c(1,500,1000), 
                               end = c(100,600,1100)))
gr2 <- GenomicRanges::GRanges(seqnames = "chr1", strand = c("-", "-", "-"),
               ranges = IRanges(start = c(1,500,1000), 
                                end = c(100,600,1100)))

out1 = resizeGRangesTranscripts(gr1, 20,80)
out2 = resizeGRangesTranscripts(gr2, 20,80)
out3 = resizeGRangesTranscripts(gr1, 110,150)
out4 = resizeGRangesTranscripts(gr2, 110,150)
test_that("Test sample output", {

  expect_equal(BiocGenerics::start(out1), c(21,500,1000))
  expect_equal(BiocGenerics::end(out1), c(100,600,1020))
  
  expect_equal(BiocGenerics::start(out2), c(1000,500,81))
  expect_equal(BiocGenerics::end(out2), c(1080,600,100))
  
  expect_equal(BiocGenerics::start(out3), 510)
  expect_equal(BiocGenerics::end(out3), 551)
  
  expect_equal(BiocGenerics::start(out4), 550)
  expect_equal(BiocGenerics::end(out4), 591)
})


