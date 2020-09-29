context("Test _each function")

# setup GRanges
gr1 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 5),
  ranges = IRanges(
    start = c(1, 401, 801, 1001, 1301),
    end = c(100, 500, 900, 1100, 1400)
  ),
  score = c(250, 100, 500, 750, 325)
)
gr2 <- GenomicRanges::GRanges(
  seqnames = "2", strand = rep("-", 6),
  ranges = IRanges(
    start = c(201, 401, 601, 801, 981, 1201),
    end = c(300, 500, 700, 920, 1100, 1250)
  ),
  score = c(1, 1000, 500, 50, 725, 300)
)

grl <- GenomicRanges::GRangesList("gr1" = gr1, "gr2" = gr2)

test_that("Test sorteach functionality", {
  out <- sorteach(grl, exonorder)
  expect_equal(GenomicRanges::start(out$gr1), c(1, 401, 801, 1001, 1301))
  expect_equal(GenomicRanges::start(out$gr2), c(1201, 981, 801, 601, 401, 201))

  out <- sorteach(grl, dplyr::desc(start))
  expect_equal(GenomicRanges::start(out$gr1), c(1301, 1001, 801, 401, 1))
  expect_equal(GenomicRanges::start(out$gr2), c(1201, 981, 801, 601, 401, 201))

  out <- sorteach(grl, score)
  expect_equal(GenomicRanges::start(out$gr1), c(401, 1, 1301, 801, 1001))
  expect_equal(GenomicRanges::start(out$gr2), c(201, 801, 1201, 601, 981, 401))
})

test_that("Test filtereach functionality", {
  out <- filtereach(grl, score > 500)
  expect_equal(as.numeric(lengths(out)), c(1, 2))

  out <- filtereach(grl, seqnames == "2")
  expect_equal(as.numeric(lengths(out)), 6)
})

test_that("Test mutateeach functionality", {
  out <- mutateeach(grl, mid = ceiling((end - start) / 2))
  expect_equal(unlist(out)$mid, c(50, 50, 50, 50, 50, 50, 50, 50, 60, 60, 25))

  out <- mutateeach(grl, coord = paste0(seqnames, ":", start, "-", end))
  expect_equal(unlist(out)[1, 2]$coord, "1:1-100")
})
