context("Test alternative splicing function")

# setup GRanges
gr1 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 5),
  ranges = IRanges::IRanges(
    start = c(1, 401, 801, 1001, 1301),
    end = c(100, 500, 900, 1100, 1400)
  ),
  transcript_id = "transcript1",
  gene_id = "geneA",
  gene_name = "gene001"
)
gr2 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 6),
  ranges = IRanges(
    start = c(201, 401, 601, 801, 981, 1201),
    end = c(300, 500, 700, 920, 1100, 1250)
  ),
  transcript_id = "transcript2",
  gene_id = "geneA",
  gene_name = "gene001"
)
gr3 <- GenomicRanges::GRanges(
  seqnames = "1", strand = rep("+", 2),
  ranges = IRanges(
    start = c(350, 801),
    end = c(500, 1150)
  ),
  transcript_id = "transcript3",
  gene_id = "geneA", 
  gene_name = "gene001"
)
grcomb1 <- c(gr1, gr2, gr3)


test_that("Test .runAS functionality", {
  out <- .runAS(grcomb1)
  expect_equal(out$AStype, c("FE", "TS", "CE", "SD", "RI", "TE", "SA", "LE"))
  expect_equal(GenomicRanges::start(out), c(201, 350, 601, 901, 901, 1101, 981, 1201))
  expect_equal(GenomicRanges::width(out), c(100, 51, 100, 20, 100, 50, 20, 50))

  GenomicRanges::strand(grcomb1) <- "-"
  out <- .runAS(grcomb1)
  expect_equal(out$AStype, c("LE", "TE", "CE", "SA", "RI", "TS", "SD", "FE"))
})

test_that("Test annotateAS robustness", {
  expect_error(labelSplicedSegment(data.frame("id" = "test")))
  expect_error(labelSplicedSegment(query_exons))
})

test_that("Test compareAS robustness", {
  expect_error(compareSplicedSegment(data.frame("id" = "test")))
  expect_warning(compareSplicedSegment(query_exons[[1]], query_exons[[3]], data.frame("id" = "test")))
  expect_warning(compareSplicedSegment(query_exons, data.frame("id" = "test")))
})
