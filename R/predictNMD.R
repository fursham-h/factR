#' Predict sensitivity of mRNA transcripts to NMD
#'
#' @param x
#' Can be a GRanges object containing exon and CDS transcript features in GTF
#' format.
#' 
#' Can be a GRangesList object containing exon features for a list of transcripts.
#' 
#' Can be a GRanges object containing exon features for a transcript.
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
#' @param which
#' List containing names of transcripts from `x` to subset for analysis
#'
#' @return
#' Dataframe with prediction of NMD sensitivity and NMD features:
#'
#' is_NMD: logical value in prediciting transcript sensitivity to NMD
#'
#' dist_to_lastEJ: Integer value indicating distance of STOP codon to last EJ
#' Values are referenced from last EJ, thus a positive value indicates upstream
#' position of STOP codon while negative value indicates downstream position
#'
#' num_of_down_EJs: Number of downstream EJs
#'
#' dist_to_downEJs: Integer value indicating distance of STOP codon to each down
#' EJs. Values are referenced from last EJ, thus a positive value indicates
#' upstream position of STOP codon while negative value indicates downstream
#' position together with distances of STOP codon to last EJ
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
predictNMD <- function(x, cds = NULL, NMD_threshold = 50, which = NULL) {

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
  totest <- .prepTotest(exons, cds, which, argnames)
  
  # run NMD analysis and return result
  return(.testNMD(exons[totest], cds[totest], NMD_threshold))
}



.testNMD <- function(x, y, threshold) {
  UTRs <- IRanges::psetdiff(unlist(range(x)), (range(y)))
  fiveUTR <- lapply(UTRs, function(x){
    if (length(x) == 0) {
      return(0)
    }
    strand <- as.character(BiocGenerics::strand(x)[1])
    x <- BiocGenerics::sort(x, decreasing = strand == "-")
    return(width(x[1,]))
  })

  
  width_to_stopcodon <- unlist(fiveUTR) + sum(width(c(y))) + 3
  rel_dist_to_stopcodon <- cumsum(width(x)) - width_to_stopcodon
  
  out <- lapply(rel_dist_to_stopcodon, function(x){
    id <- ifelse(!is.null(names(x)), names(x)[1], 'transcript')
    x <- sort(x, decreasing = T)
    threeUTR <- x[1]
    dist_to_last <- x[2]
    is_NMD <- ifelse(dist_to_last > threshold, T, F)
    dist_to_eachEJ <- rev(x[-1][x[-1] > 0])
    
    return(tibble::tibble('transcript' = id, 
                          'is_NMD' = is_NMD,
                          'dist_to_lastEJ' = dist_to_last,
                          'num_of_downEJs' = length(dist_to_eachEJ),
                          'dist_to_downEJs' = paste(dist_to_eachEJ, collapse = ','),
                          'threeUTRlength' = threeUTR))
  }) %>% dplyr::bind_rows()
  return(out)
}


.prepTotest <- function(x, y, which, arg) {
  totest <- names(x) # prepare vector with names for testing
  if (!is.null(which)) {
    totest <- intersect(totest, which)
    if (length(totest) == 0) {
      rlang::abort(sprintf("transcript names in `which` is not found in `%s`",
                           arg[1]))
    } else if (length(totest) != length(which)) {
      num_unmatched <- length(which) - length(which_matched)
      rlang::warn(sprintf(
        "%s transcripts in `which` is missing from `%s%`",
        num_unmatched, arg[1]
      ))
    }
  }
  txwithcds <- intersect(totest, names(y)) # subset transcripts with cds info
  if (length(txwithcds) == 0) {
    rlang::abort("all transcripts have missing cds info")
  }
  if (length(txwithcds) < length(totest)) {
    totest <- txwithcds
    skiptest <- length(x) - length(totest)
    rlang::warn(sprintf(
      "%s transcript(s) have missing cds info and was not analyzed",
      skiptest))
  }
  return(totest)
}