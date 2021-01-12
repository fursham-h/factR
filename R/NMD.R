identifyNMDexons_old <- function(x, NMD.result = NULL) {
    
    # checks etc, will add later
    ## check for x format, check for CDS
    ## check for NMD.result, if yes check if all tx have NMD annotation
    
    # group transcripts by same CDS start
    x.grouped <- .grouptxbycdsstart(x)
    
    # remove groups with only 1 transcript and nest
    x.grouped <- x.grouped %>% 
        dplyr::group_by(coord) %>% 
        dplyr::filter(dplyr::n() > 1) %>% 
        dplyr::ungroup()
    
    # append groupings to GTF and prepare for AS labelling
    x.appended <- x %>% 
        as.data.frame() %>% 
        dplyr::left_join(x.grouped, by = "transcript_id") %>% 
        dplyr::filter(!is.na(coord)) %>% 
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # get alternative segments for each group of transcripts
    AS.out <- labelSplicedSegment(x.appended, group.by = "coord")
    
    # generate exhaustive list of AS exons and transcripts that contain/missing
    exons.list <- .getexonslist(AS.out, x.grouped)
    
    # append NMD information
    exons.list <- exons.list %>% 
        dplyr::left_join(NMD.result %>% dplyr::select(transcript, is_NMD),
                         by = c("transcript_id"="transcript"))
    
    # MAJOR TEST: see if theory works
    ## Poixon exon: all true splice must return true NMD
    ## ORF-maintiang: all TRUE splice must return false and one true in FALSE splice
    exons.list %>% 
        dplyr::select(-transcript_id) %>% 
        dplyr::group_by(ASexon, spliced) %>% 
        dplyr::mutate(data = list(is_NMD)) %>% 
        dplyr::select(-is_NMD) %>% 
        dplyr::distinct() %>% 
        tidyr::spread(spliced, data) %>% 
        dplyr::mutate(NMDclass = "no") %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(NMDclass = ifelse(all(`TRUE`), "Poison", NMDclass)) %>% 
        dplyr::mutate(NMDclass = ifelse(sum(`TRUE`)==0 & any(`FALSE`), 
                                        "ORF-maintain", NMDclass)) 
    
    
    return(exons.list)
}

.grouptxbycdsstart <- function(x) {
    
    # create grangeslist of CDS and sort by exon order
    cds.grl <- S4Vectors::split(x[x$type == "CDS"], ~transcript_id)
    cds.grl <- sorteach(cds.grl, exonorder)
    
    # return translation start coordinate for each transcript
    cds.atg.grl <- GenomicRanges::resize(cds.grl, width = 1) %>% 
        filtereach(dplyr::row_number() == 1)
    
    # group transcripts by translation start
    groupsbystart <- cds.atg.grl %>% 
        as.data.frame() %>% 
        dplyr::mutate(coord = paste0(seqnames, "_", start, "_", strand)) %>% 
        dplyr::select(coord, transcript_id) %>% 
        dplyr::distinct()
        
    # return df of grouped transcripts
    return(groupsbystart)
}

.getexonslist <- function(AS.out, x.grouped) {
    
    # get a list of AS exons with its parent CDS groupings
    exons.list <- AS.out %>% 
        as.data.frame() %>% 
        tidyr::unite("ASexon", c(seqnames, start, end, strand, AStype)) %>% 
        dplyr::select(ASexon, transcript_id) 
    
    # merge with x.grouped to annotate which transcripts contain or missing of
    # the exons
    full.exons.list <- exons.list %>% 
        dplyr::full_join(x.grouped, by = "transcript_id") %>% 
        dplyr::select(ASexon, coord) %>% 
        dplyr::distinct() %>% 
        dplyr::left_join(x.grouped, by = "coord") %>% 
        dplyr::select(ASexon, transcript_id) %>% 
        dplyr::left_join(exons.list %>% mutate(spliced = TRUE),
                         by = c("ASexon", "transcript_id")) %>% 
        tidyr::replace_na(list(spliced=FALSE))
    
    return(full.exons.list)
    
}