identifyNMDexons_old2 <- function(x, NMD.result = NULL) {
    
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
    
    # retrieve transcripts containing reduced intron
    tx.intronsreduced <- GenomicRanges::findOverlaps(intronsreduced, range(exonsbytx))
    tx.intronsreduced <- as.data.frame(tx.intronsreduced)
    tx.intronsreduced$transcript_id <- names(exonsbytx)[tx.intronsreduced$subjectHits]
    
    # generate non-redundant table of exons and its groupings
    exons.list <- AS.exons %>% 
        as.data.frame() %>% 
        tidyr::unite("ASexon", c(seqnames, start, end, strand, AStype)) %>% 
        dplyr::select(ASexon, transcript_id, group)
    
    # for each AS exon, determine which transcripts splices or skips that exon
    full.exons.list <- exons.list %>%
        dplyr::select(ASexon, group) %>% 
        dplyr::distinct() %>% 
        dplyr::left_join(tx.intronsreduced[,c(1,3)], 
                         by = c("group"="queryHits")) %>% 
        dplyr::left_join(exons.list %>% dplyr::mutate(spliced="spliced"),
                         by = c("ASexon", "group", "transcript_id")) %>% 
        tidyr::replace_na(list(spliced="skipped"))
    
    # fill NMD information
    full.exons.NMD.list <- full.exons.list %>% 
       dplyr::left_join(NMD.result[c("transcript","is_NMD")],
                         by = c("transcript_id"="transcript"))
    
    # testinggg
    NMD.exons.annotate <- full.exons.NMD.list %>% 
        dplyr::select(-transcript_id) %>% 
        dplyr::group_by(ASexon, spliced) %>% 
        dplyr::mutate(data = list(is_NMD)) %>% 
        dplyr::select(-is_NMD) %>% 
        dplyr::distinct() %>% 
        tidyr::spread(spliced, data) %>% 
        dplyr::mutate(NMDclass = "no") %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(NMDclass = ifelse(all(spliced), "Poison", NMDclass)) %>% 
        dplyr::mutate(NMDclass = ifelse(all(skipped), "ORF-maintain", NMDclass)) 
    
    return(NMD.exons.annotate)
}