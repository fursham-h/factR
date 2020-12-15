#' Subset GTF GRanges object for de novo transcripts
#' 
#' @description
#' subsetNewTranscripts will compare query and reference GTF GRanges and by
#' default, return query transcripts with different exon structures from 
#' transcripts in reference. In addition, query GTF can be further filtered
#' by removing transcripts of the same name (by.names = TRUE) and of same CDS info
#' (by.CDS = TRUE) as reference transcripts
#' 
#' 
#' @param query 
#' GRanges object containing query GTF data.
#' @param ref 
#' GRanges object containing reference GTF data.
#' @param by.names
#' Whether to remove query transcripts with same name as reference transcripts. 
#' Default: FALSE
#' @param by.CDS
#' Whether to remove query transcripts with CDS structure as reference transcripts.
#' Default: FALSE
#' 
#' 
#' 
#' @return
#' Filtered GRanges GTF object
#' 
#' @export
#'
#' @examples
#' subsetNewTranscripts(matched_query_gtf, ref_gtf)
subsetNewTranscripts <- function(query, ref, by.names = FALSE, by.CDS = FALSE) {
  
  # catch missing args
  mandargs <- c("query", "ref")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    rlang::abort(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  
  # retrieve input object names
  argnames <- as.character(match.call())[-1]
  
  # carry out input checks
  .checkinputs(query, ref, argnames)
  
  # run subsetting and return new transcripts
  return(.subsetTranscripts(query, ref, by.names, by.CDS))
}

.checkinputs <- function(query, ref, argnames) {
  # check inputs are gtf
  if (any(!is_gtf(query, ref))) {
    rlang::abort(sprintf(
      "%s is/are not gtf GRanges",
      paste(argnames[1:2][!is_gtf(query, ref)], collapse = ",")
    ))
  }
  
  
  # catch unmatched seqlevels
  if (suppressWarnings(!has_consistentSeqlevels(query, ref))) {
    rlang::abort(sprintf(
      "`%s` and `%s` has unmatched seqlevel styles. 
Try running: %s <- matchChromosomes(%s, %s)",
      argnames[1], argnames[2], argnames[1], argnames[1], argnames[2]
    ))
  } else {
    outmsg <- has_consistentSeqlevels(query,ref)
  }
}

.subsetTranscripts <- function(query, ref, by.names, by.CDS){
  
  # Default filtering, by exon coordinates
  # create exons by transcripts
  query_exons <- S4Vectors::split(query[query$type == "exon"], ~transcript_id)
  ref_exons <- S4Vectors::split(ref[ref$type == "exon"], ~transcript_id)
  
  # search for exact query and ref matches
  fulloverlap <- GenomicRanges::findOverlaps(query_exons, ref_exons,
                                             type = "equal", select = "first"
  )
  
  # return transcripts that are not found in reference
  query <- query[!query$transcript_id %in% names(query_exons)[!is.na(fulloverlap)]]
  
  # Optional fitering, by name and CDS coordinates
  if (by.names) {
    query <- query[!query$transcript_id %in% ref$transcript_id]
  }
  if (by.CDS) {
    # create exons by transcripts
    query_CDS <- S4Vectors::split(query[query$type == "CDS"], ~transcript_id)
    ref_CDS <- S4Vectors::split(ref[ref$type == "CDS"], ~transcript_id)
    
    # search for exact query and ref matches
    fullCDSoverlap <- GenomicRanges::findOverlaps(query_CDS, ref_CDS,
                                               type = "equal", select = "first")
    
    # return transcripts that are not found in reference
    query <- query[!query$transcript_id %in% names(query_CDS)[!is.na(fullCDSoverlap)]]
  }
  
  return(query)
}