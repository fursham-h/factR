context("Test matching functions")

data(query_gtf)
data(ref_gtf)

matched1 <- matchChromosomes(query_gtf, ref_gtf)
matched2 <- matchGeneInfo(matched1, ref_gtf)

test_that("Test the structure of outputs", {
  expect_equal(is(matched1, "GRanges"), TRUE)
  expect_equal(is(matched2, "GRanges"), TRUE)
  expect_equal(length(matched1), 56)
  expect_equal(length(matched2), 56)
  expect_equal(names(S4Vectors::mcols(matched2)), c(
    "type", "transcript_id", "gene_id", "gene_name",
    "old_gene_id", "match_level"
  ))
})

test_that("Test the value of outputs", {
  expect_equal(as.character(GenomeInfoDb::seqnames(matched1))[1], "chr10")
  expect_equal(unique(mcols(matched2)$gene_id), "ENSMUSG00000006498.14")
  expect_equal(unique(mcols(matched2)$match_level), 4)
})
