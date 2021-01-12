identifyNMDexons <- function(x, NMD.result = NULL) {
    
    # checks etc, will add later
    ## check for x format, check for CDS
    ## check for NMD.result, if yes check if all tx have NMD annotation
    
    # retain transcripts with CDS only
    x <- x[x$transcript_id %in% x[x$type=="CDS"]$transcript_id]
    
    # get a list of all alternative exons
    AS.exons <- labelSplicedSegment(x)
    
    # get a list of reduced introns (borrowed from splicinFunctions.R)
    exonsbytx <- S4Vectors::split(x[x$type == "exon"], ~transcript_id)
    intronsbytx <- GenomicRanges::psetdiff(BiocGenerics::unlist(range(exonsbytx)), exonsbytx)
    intronsreduced <- GenomicRanges::reduce(unlist(intronsbytx))
    
    # group AS exons by reduced intron
    AS.exons$group <-  GenomicRanges::findOverlaps(AS.exons, intronsreduced, select = "first")
}