
.calcCoverage <- function(tx1, tx2, over) {
  chrom <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(tx1)))
  cov <- suppressWarnings(GenomicRanges::coverage(c(tx1, tx2)))
  index <- which(names(cov) == chrom)
  cov <- cov[[index]]
  cov_val <- S4Vectors::runValue(cov)
  cov_len <- S4Vectors::runLength(cov)
  
  if (over == "mean") {
    denom <- sum(BiocGenerics::width(tx1), BiocGenerics::width(tx2)) / 2
  } else if (over == "query") {
    denom <- sum(BiocGenerics::width(tx1))
  } else if (over == "ref") {
    denom <- sum(BiocGenerics::width(tx2))
  }
  return(sum(cov_len[cov_val == 2]) / denom)
}


# 
# # prepare a dict of stop codons for pattern matching
# list_stopcodons <- Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
# pdict_stopcodons <- Biostrings::PDict(list_stopcodons)

# # testingggg
# out <- BiocParallel::bpmapply(function(x, y, z) {
#   
#   strand = as.character(strand(y)[1])
#   y <- sort(y, decreasing = strand == '-')
#   z <- sort(z, decreasing = strand == '-')
#   z <- resize(z, width = width(z) - z$phase, fix = 'end')
#   a <- subsetByOverlaps(z, y, type = 'within')
#   
#   if (length(a) == 0) {
#     return(NULL)
#   }
#   
#   starttoend <- dplyr::bind_cols(as.data.frame(range(y)), as.data.frame(range(a[1]))) %>%
#     dplyr::mutate(newstart = ifelse(strand == '-', start ,start1)) %>%
#     mutate(newend = ifelse(strand == '-', end1 ,end)) %>% 
#     dplyr::select(seqnames, start = newstart, end = newend, strand) %>% 
#     GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
#   
#   starttoend <- GenomicRanges::intersect(y, starttoend)
#   queryseq <- unlist(BSgenome::getSeq(fasta, starttoend))
#   
#   
#   # search for in-frame stop codons
#   stopcodons <- Biostrings::matchPDict(pdict_stopcodons, queryseq) %>%
#     unlist() %>%
#     as.data.frame() %>%
#     dplyr::filter(end %% 3 == 0) %>%
#     dplyr::arrange(start)
#   
#   # return if no in-frame stop codons are found
#   if (nrow(stopcodons) == 0) {
#     return(NULL)
#   } 
#   threeUTRlength <- length(queryseq) - stopcodons[1,2] + 3
#   cds <- resizeTranscript(starttoend, end = threeUTRlength)
#   
#   # return if no in-frame stop codons are found
#   if (length(cds) == 0) {
#     return(NULL)
#   } 
#   cds$type <- 'CDS'
#   cds$transcript_id <- x
#   cds$phase <- rev(cumsum(rev(width(cds)) %% 3) %% 3) 
#   return(as.data.frame(cds))
#   
# }, query2ref$transcript_id, 
# query_exons[query2ref$transcript_index], 
# ref_cds[query2ref$ref_transcript_id],
# BPPARAM = BiocParallel::MulticoreParam(), SIMPLIFY = F
# ) %>% dplyr::bind_rows()
