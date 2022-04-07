context("Test identifying NMD exon function")

query_cds_gtf <- suppressMessages(buildCDS(matched_query_gtf, ref_gtf, Mmusculus))

test_that("Test identiyNMDexons checks", {
    expect_error(identifyNMDexons(query_gtf))
    expect_error(identifyNMDexons(query_exons, Mmusculus))
    expect_error(identifyNMDexons(query_gtf, Mmusculus))
})

test_that("Test identiyNMDexons output", {
    out <- suppressMessages(identifyNMDexons(query_cds_gtf, Mmusculus))
    expect_equal(length(out), 1)
    expect_equal(as.data.frame(out)$start, 79862014)
    expect_equal(as.data.frame(out)$end, 79862047)
    expect_equal(out$Exontype, "NMD-Repressive")
    #expect_equal(out$phastCons60way.UCSC.mm10, 1)
})
