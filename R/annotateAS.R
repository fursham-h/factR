annotateAS <- function(x, as.data.frame = FALSE, append = FALSE) {
  # catch missing args
  mandargs <- c("x")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    rlang::abort(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  
  # retrieve input object names and carry out checks
  argnames <- as.character(match.call())[-1]
  .ASchecks(x, argnames)
  
  # run Alternative splicing annotation function
  annotatedAS <- .annotate(x)
  
  # return
  if (as.data.frame) {
    return(annotatedAS %>% as.data.frame() %>%
      dplyr::mutate(coord = paste0(seqnames, ':', start, '-', end)) %>%
      dplyr::select(gene_id, transcript_id, coord, AStype))
  } else if (append) {
    return(c(x, annotatedAS))
  } else {
    return(annotatedAS)
  }
}

.ASchecks <- function(x, names) {
  # check inputs are gtf
  if (!is_gtf(x)) {
    rlang::abort(sprintf(
      "`%s` is not a gtf GRanges object",
      names[1]
    ))
  }
}


.annotate <- function(x) {
  
  x <- x %>% as.data.frame() %>%
    dplyr::filter(type == 'exon') %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(pos = dplyr::row_number()) %>% 
    dplyr::mutate(pos = ifelse(pos == 1, 'First', pos)) %>%
    dplyr::mutate(pos = ifelse(pos == dplyr::n(), 'Last', pos)) %>%
    dplyr::select(seqnames, start, end, strand, type, gene_id, transcript_id, pos) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  # get reduced intron boundaries
  exonsbytx <- S4Vectors::split(x, ~transcript_id)
  intronsbytx <- GenomicRanges::psetdiff(unlist(range(exonsbytx)), exonsbytx)
  intronsreduced <- GenomicRanges::reduce(unlist(GenomicRanges::psetdiff(unlist(range(exonsbytx)), exonsbytx)))
  
  # get exons that overlap with reduced introns
  altexons <- IRanges::findOverlapPairs(x, intronsreduced)
  
  # annotate AS exons
  altannotate <- altexons %>% as.data.frame() %>%
    dplyr::mutate(AStype = 'CE') %>%
    dplyr::mutate(AStype = ifelse(first.X.start < second.start & first.X.end < second.end & !first.pos %in% c('First','Last'), 'SD', AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start > second.start & first.X.end > second.end & !first.pos %in% c('First','Last'), 'SA', AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start < second.start & first.X.end < second.end & first.pos %in% c('First','Last'), 'LE', AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start > second.start & first.X.end > second.end & first.pos %in% c('First','Last'), 'FE', AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start < second.start & first.X.end > second.end, 'RI', AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start > second.start & first.X.end < second.end & first.pos == 'First', 'FE', AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start > second.start & first.X.end < second.end & first.pos == 'Last', 'LE', AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.strand == '-', chartr("DAFL", "ADLF", AStype), AStype))
  
  # get segments that fall within intron and annotate that segment
  altexons <- GenomicRanges::pintersect(altexons)
  if (length(altexons) == 0) {
    altexons$pos <- altexons$hit <- NULL
    return(altexons)
  } else {
    altexons$type <- 'AS'
    altexons$AStype <- altannotate$AStype
    altexons$pos <- altexons$hit <- NULL
    
    return(altexons)
  }
}
