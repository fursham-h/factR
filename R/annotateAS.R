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
  annotatedAS <- .runAS(x[x$type == 'exon'])
  
  # return
  if (as.data.frame) {
    return(annotatedAS %>% as.data.frame() %>%
      dplyr::mutate(coord = paste0(seqnames, ':', start, '-', end)) %>%
      dplyr::select(transcript_id, coord, AStype))
  } else if (append) {
    return(c(x, annotatedAS))
  } else {
    return(annotatedAS)
  }
}

#' Compare and classify alternative spliced segments
#'
#' @param exons
#' In pair-wise mode: GRanges object containing exons for a particular transcript.
#'
#' In intra-list mode: GRangesList bject containing exons for each transcripts.
#' Transcripts will be paired and compared based on its gene family or groups.
#' In order to do so, object has to contain gene_id attribute. Alternatively, user
#' may provide a dataframe with a list of transcripts and its groupings as a `groupings`
#' argument. See `groupings`.
#' @param ...
#' In pair-wise mode, argument is a GRanges object containing exons for a
#' particular transcript as a comparison to `exons`
#' @param groupings
#' Dataframe describing the groupings of the transcripts in `exons`. Ideally, Transcripts should be
#' grouped by gene families. Therefore, first column in the dataframe is a list of gene_id or
#' gene_names and the second column is the names of transcripts in `exons` which fall into
#' the gene groupings. This argument is optional if gene_id metadata is present in `exons`
#'
#'
#' @return
#' In pairwise mode: a GRangesList object with included and skipped segments in `exons`
#' as compared to ...
#'
#' In intralist mode: a dataframe with coordinates and information of alternative segments
#' between every transcripts in its groupings.
#' @export
#'
compareAS <- function(exons, ...) {
  
  # catch missing args
  mandargs <- c("exons")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    rlang::abort(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  # define global variables
  tx.id <- index.y <- index.x <- tx.id.x <- tx.id.y <- NULL
  Gene <- seqnames <- compare.to <- coord <- strand <- NULL
  . <- AS.type <- AS.direction <- NULL
  
  argnames <- as.character(match.call())[-1]
  
  if (is(exons, "GRanges")) {
    if (!'transcript_id' %in% names(S4Vectors::mcols(exons))) {
      exons$transcript_id <- 'transcript0'
    } 
    exons <- S4Vectors::split(exons, ~transcript_id)
  }
  if (!is(exons, "GRangesList")) {
    rlang::abort()
  }
  
  # if a second GRanges object is given, carry out pairwise comparison
  if (length(list(...)) > 0) {
    # multiple comparisons not supported. can be a feature in the future
    # if (length(list(...)) > 1) {
    #   # idea for multiple comparison. convert grangeslist to df for comparison
    #   rlang::warn("Multiple comparisons is not yet supported in pair-wise mode. First item in ... used")
    # }
    dots <- list(...)
    newdots <- dots[unlist(lapply(dots, function(x){is(x,"GRanges")}))]
    
    if (length(newdots) == 0) {
      # return warning and proceed to intra-list mode
    } else {
      if (length(newdots) < length(dots)) {
        # return warning for elements which are not GRanges
      }
      newdots <- as(newdots, "GRangesList")
      
      newdotsmeta <- newdots %>% as.data.frame()
      if (!'transcript_id' %in% names(newdotsmeta)) {
        names(newdots) <- paste0('transcript', as.character(c(1:length(newdots))))
      } else {
        newdots <- newdots %>%
          as.data.frame() %>%
          dplyr::mutate(transcript_id = ifelse(is.na(transcript_id), paste0('transcript',group), transcript_id)) %>%
          dplyr::mutate(group_name = transcript_id) %>%
          dplyr::select(-group) %>%
          GenomicRanges::makeGRangesListFromDataFrame(split.field = 'group_name', keep.extra.columns = T)
      }
    }
    exons <- c(exons, newdots)
  }

  return(S4Vectors::split(.runAS(exons), ~transcript_id))
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


.runAS <- function(x) {
  
  x <- x %>% as.data.frame() %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(pos = dplyr::row_number()) %>% 
    dplyr::mutate(pos = ifelse(pos == 1, 'First', pos)) %>%
    dplyr::mutate(pos = ifelse(pos == dplyr::n(), 'Last', pos)) %>%
    dplyr::select(seqnames, start, end, strand, transcript_id, pos) %>%
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
