#' Subset GTF GRanges object for de novo transcripts
#' 
#' @description
#' subsetNewTranscripts will compare query and reference GTF GRanges and return
#' transcripts which are only found in query (de novo).
#' 
#' 
#' @param query 
#' GRanges object containing query GTF data.
#' @param ref 
#' GRanges object containing reference GTF data.
#'
#' @return
#' GRanges object containing exon coordinates from transcripts which are only
#' found in query.
#' 
#' @export
#'
#' @examples
#' subsetNewTranscripts(matched_query_gtf, ref_gtf)
subsetNewTranscripts <- function(query, ref) {
  
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
  return(.subsetTranscripts(query, ref))
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
  }
}

.subsetTranscripts <- function(query, ref){
  
  # create exons by transcripts
  query_exons <- S4Vectors::split(query[query$type == "exon"], ~transcript_id)
  ref_exons <- S4Vectors::split(ref[ref$type == "exon"], ~transcript_id)
  
  # search for exact query and ref matches
  fulloverlap <- GenomicRanges::findOverlaps(query_exons, ref_exons,
                                             type = "equal", select = "first"
  )
  
  # return transcripts that are not found in reference
  return(query[!query$transcript_id %in% names(query_exons)[fulloverlap]])
}