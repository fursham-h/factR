#' Obtain percent coverage between query and reference transcripts
#'
#' @param query
#' GRangesList object containing exons for each query transcript. Transcripts
#' have to be listed in query2ref dataframe, else coverage values will not be
#' calculated
#' @param ref
#' GRangesList object containing CDS for each reference transcript.
#' @param query2ref
#' Dataframe with at least 2 columns: query transcript_id and its
#' reference transcript_id. Query and ref transcript_ids have to match transcript
#' names in query and refCDS objects. Transcripts with missing corrresponding
#' GRanges object will return error
#' @param ids
#' Numeric vector stating which columns of query2ref dataframe contain the
#' query and reference transcript_ids respectively.
#' @param return
#' If query2ref contain query transcripts with multiple comparisons, function
#' will return comparisons with 'best' coverage value by default. Alternatively,
#' a full report can be returned by setting argument to 'all'
#'
#' @return
#' Dataframe from query2ref with coverage values appended as a new column
#' @export
#' @author Fursham Hamid
#'
#' @examples
#' getCoverages(query_exons, ref_exons, q2r)
getCoverages <- function(query, ref, query2ref, ids = c(1, 2),
                         return = c("best", "all")) {

  # Plans: set 'over' arg to allow user to choose the denominator

  # catch missing args
  mandargs <- c("query", "ref", "query2ref")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  
  # define global variables
  unnanotatedq <- unnanotatedr <- coverage <- NULL

  # catch unmatched seqlevels
  if (GenomeInfoDb::seqlevelsStyle(query) != GenomeInfoDb::seqlevelsStyle(ref)) {
    stop("query and ref has unmatched seqlevel styles. try matching matchSeqLevels function")
  }

  # sanity check if all tx in q2r have GRanges object
  if (!all(query2ref[[ids[1]]] %in% names(query))) {
    missing <- sum(!query2ref[[ids[1]]] %in% names(query))
    stop(sprintf(
      "%s query transcripts have missing GRanges object",
      missing
    ))
  }
  if (!all(query2ref[[ids[2]]] %in% names(ref))) {
    missing <- sum(!query2ref[[ids[2]]] %in% names(ref))
    stop(sprintf(
      "%s reference CDSs have missing GRanges object",
      missing
    ))
  }

  # sanity check if query and ref names are in q2f df
  if (all(!names(query) %in% query2ref[[ids[1]]])) {
    unannotatedq <- sum((!names(query) %in% query2ref[ids[1]]))
    rlang::warn(sprintf(
      "%s query transcript ids were missing from query2ref df",
      unnanotatedq
    ))
  }
  if (all(!names(ref) %in% query2ref[[ids[2]]])) {
    unannotatedr <- sum((!names(ref) %in% query2ref[ids[2]]))
    rlang::warn(sprintf(
      "%s reference CDS ids were missing from query2ref df",
      unnanotatedr
    ))
  }

  # extract colnames and prepare outputCDS
  txname <- names(query2ref)[ids[1]]
  refname <- names(query2ref)[ids[2]]

  # get Coverage values for all comparisons
  out <- BiocParallel::bpmapply(function(x, y) {
    covrep <- countCoverage_(query[[x]], ref[[y]])
    return(covrep)
  }, query2ref[[txname]], query2ref[[refname]],
  BPPARAM = BiocParallel::MulticoreParam()
  )
  query2ref$coverage <- out

  # return best coverage for each tx by default
  if (return[1] == "best") {
    query2ref <- query2ref %>%
      dplyr::arrange(!!as.symbol(txname), dplyr::desc(coverage)) %>%
      dplyr::distinct(!!as.symbol(txname), .keep_all = T)
  }
  return(query2ref)
}

countCoverage_ <- function(tx1, tx2) {
  chrom <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(tx1)))
  cov <- GenomicRanges::coverage(c(tx1, tx2))
  index <- which(names(cov) == chrom)
  cov <- cov[[index]]
  cov_val <- S4Vectors::runValue(cov)
  cov_len <- S4Vectors::runLength(cov)
  return(sum(cov_len[cov_val == 2]) / ((sum(cov_len[cov_val == 1]) + (sum(cov_len[cov_val == 2]) * 2)) / 2))
}
