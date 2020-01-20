#' Predict sensitivity of mRNA transcripts to NMD
#'
#' @param x
#' Can be a GRanges object containing exon and CDS transcript features in GTF
#' format.
#' 
#' Can be a GRangesList object containing exon features for a list of transcripts.
#' If so, `cds` argument have to be provided.
#' 
#' Can be a GRanges object containing exon features for a transcript. If so, `cds` 
#' argument have to be provided.
#' 
#' @param ...
#' Logical conditions to pass to dplyr::filter to subset transcripts for analysis.
#' Variables are metadata information found in `x` and multiple conditions can be 
#' provided delimited by comma. Example: transcript_id == "transcript1"
#' 
#' @param cds
#' If `x` is a GRangesList object, `cds` has to be a GRangesList containing CDS features
#' for the list of transcripts in `x`. List names in `x` and `cds` have to match.
#' 
#' If `x` is a GRanges object, `cds` has to be a GRanges containing CDS features
#' for the transcript in `x`.
#' 
#' @param NMD_threshold
#' Minimum distance of stop_codon to last exon junction (EJ) which triggers NMD.
#' Default = 50bp

#'
#' @return
#' Dataframe with prediction of NMD sensitivity and NMD features:
#'
#' is_NMD: logical value in prediciting transcript sensitivity to NMD
#'
#' dist_to_lastEJ: Integer value of the number of bases between the first
#' base of the stop_codon to the last base of EJ. A positive value indicates that
#' the last EJ is downstream of the stop_codon.
#'
#' num_of_down_EJs: Number of EJs downstream of the stop_codon.
#'
#' dist_to_downEJs: Concatenated integer values of the number of bases between 
#' the first base of the stop_codon to the last base of each downstream EJs.
#' @export
#' @author Fursham Hamid
#'
#' @examples
#'
#' ### Examples with GRanges objects
#' predictNMD(query_exons$transcript1, query_cds$transcript1) # NMD-insensitive
#' predictNMD(query_exons$transcript3, query_cds$transcript3) # NMD-sensitive
#'
#' ### Examples with GRangesList object
#' predictNMD(query_exons, query_cds)
#' predictNMD(query_exons, query_cds, which = c("transcript1", "transcript3"))
predictNMD <- function(x, ..., cds = NULL, NMD_threshold = 50) {

  # catch missing args
  mandargs <- c("x")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    rlang::abort(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  # define global variable
  #is_NMD <- NULL
  
  # retrieve input object names
  argnames <- as.character(match.call())[-1]

  # check if exons is a gtf or both exons and cds are GR or GRlist
  if (is_gtf(x)) {
    exons <- S4Vectors::split(x[x$type == 'exon'], ~transcript_id)
    cds <- S4Vectors::split(x[x$type == 'CDS'], ~transcript_id)
  } else if (all(is(x) %in% is(cds))) {
    if (is(x, "GRanges")) {
      exons <- GRangesList('transcript' = x)
      cds <- GRangesList('transcript' = cds)
    } else if (is(x, "GRangesList")) {
      exons <- x
    } else {
      errorobj <- paste(unique(c(is(x)[1], is(cds)[2])), collapse = ', ')
      rlang::abort(sprintf(
        "input object types %s not supported", errorobj))
    }
  } else {
    txtype <- is(x)[1]
    cdstype <- is(cds)[1]
    rlang::abort(sprintf(
      "`%s` is type %s but `%s` is type %s",
      argnames[2],cdstype,argnames[1], txtype
    ))
  }
  
  # catch unmatched seqlevels
  if (GenomeInfoDb::seqlevelsStyle(exons)[1] != GenomeInfoDb::seqlevelsStyle(cds)[1]) {
    rlang::abort(sprintf(
      "`%s` and `%s` has unmatched seqlevel styles. try running: 
      \t\t%s <- matchSeqLevels(%s, %s)",
      argnames[1], argnames[2], argnames[1], argnames[1], argnames[2])
  )}

  # subset list for `which` and check if all exons have a cds entry
  totest <- .prepTotest(exons, cds, argnames, ...)
  
  # run NMD analysis and return result
  return(.testNMD(exons[totest], cds[totest], NMD_threshold))
}



.testNMD <- function(x, y, threshold) {
  out <- tibble::tibble('transcript' = as.character(), 
                 'stop_to_lastEJ' = as.double(),
                 'num_of_downEJs' = as.integer(),
                 'stop_to_downEJs' = as.character(),
                 'threeUTRlength' = as.double(),
                 'is_NMD' = as.logical())
  
  # sort all exons by strand first
  x <- sorteach(x, exonorder)

  toStopRange <- dplyr::bind_cols(as.data.frame(range(x)), as.data.frame(range(y))) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(newstart = ifelse(strand == '-', start1 ,start)) %>%
    dplyr::mutate(newend = ifelse(strand == '-', end ,end1)) %>% 
    dplyr::select(group:seqnames, start = newstart, end = newend, strand) %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  toStopWidth <- sum(width(GenomicRanges::pintersect(x, toStopRange)))
  EJtoStop <- cumsum(width(x)) - toStopWidth
  
  out <- dplyr::bind_rows(out, lapply(EJtoStop, function(x){
    id <- ifelse(!is.null(names(x)), names(x)[1], 'transcript')
    x <- sort(x, decreasing = T)
    threeUTR <- x[1]
    dist_to_last <- x[2]
    is_NMD <- ifelse(dist_to_last > threshold, T, F)
    dist_to_eachEJ <- rev(x[-1][x[-1] > 0])
    
    
    
    return(tibble::tibble('transcript' = id, 
                          'stop_to_lastEJ' = dist_to_last,
                          'num_of_downEJs' = length(dist_to_eachEJ),
                          'stop_to_downEJs' = paste(dist_to_eachEJ, collapse = ','),
                          'threeUTRlength' = threeUTR,
                          'is_NMD' = is_NMD))
  })) 
  return(out)
}


.prepTotest <- function(x, y, arg, ...) {
  totest <- names(x) # prepare vector with names for testing
  if (!missing(...)) {
    
    totest <- tryCatch(
      names(filtereach(x, ...)),
      error = function(e){
        rlang::abort(sprintf(
        "Metadata in ... not found in `%s`",
                                   arg[1]))}
    )
    if (length(totest) == 0) {
      return(NULL)
    } 
  }
  txwithcds <- intersect(totest, names(y)) # subset transcripts with cds info
  if (length(txwithcds) == 0) {
    rlang::abort("all transcripts have missing cds info")
  }
  if (length(txwithcds) < length(totest)) {
    skiptest <- length(totest) - length(txwithcds)
    rlang::warn(sprintf(
      "%s transcript(s) have missing cds info and was not analyzed",
      skiptest))
    totest <- txwithcds
  }
  return(totest)
}