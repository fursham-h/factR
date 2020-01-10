#' Predict sensitivity of mRNA transcripts to NMD
#'
#' @param exons
#' GRanges object or GRangesList object containing exons
#' for each transcript. To search for NMD-inducing features, transcripts have
#' to be coding and thus contain a cds information of the same transcript name
#' @param cds
#' GRanges object or GRangesList object containing coding regions (CDS)
#' for each transcript. Object type must match `exons` object type. GRangesList 
#' must have names that match names in exons, else transcript will not be analysed
#' @param NMD_threshold
#' Minimum distance of STOP codon to last exon junction (EJ) which triggers NMD.
#' Default = 50bp
#' @param which
#' List containing names of transcripts from `exons` to filter for analysis
#' @param return
#' If exons and cds are GRangesList, returns results for all transcripts (defailt) or
#' only NMD-sensitive (default) transcripts
#'
#' @return
#' List with prediction of NMD sensitivity and statistics:
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
#' EJs Values are referenced from last EJ, thus a positive value indicates
#' upstream position of STOP codon while negative value indicates downstream
#' position together with distances of STOP codon to last EJ
#' @export
#' @author Fursham Hamid
#'
#' @examples
#'
#' ### To visualize transcripts
#' library(wiggleplotr)
#' plotTranscripts(query_exons, query_cds)
#'
#' ### Examples with GRanges objects
#' predictNMD(query_exons$transcript1, query_cds$transcript1) # NMD-insensitive
#' predictNMD(query_exons$transcript3, query_cds$transcript3) # NMD-sensitive
#'
#' ### Examples with GRangesList object
#' predictNMD(query_exons, query_cds)
#' predictNMD(query_exons, query_cds, return = "NMD")
#' predictNMD(query_exons, query_cds, which = c("transcript1", "transcript3"))
predictNMD <- function(exons, cds, NMD_threshold = 50,
                       which = NULL, return = c("all", "NMD")) {

  # catch missing args
  mandargs <- c("exons", "cds")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    rlang::abort(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  # define global variable
  is_NMD <- NULL
  
  # retrieve input object names
  argnames <- as.character(match.call())[-1]

  # check if exons and cds are GR or GRlist
  if (all(is(exons) %in% is(cds))) {
    if (is(exons, "GRanges")) {
      intype <- "gr"
    } else if (is(exons, "GRangesList")) {
      intype <- "grl"
    } else {
      errorobj <- paste(unique(c(is(exons)[1], is(cds)[2])), collapse = ', ')
      rlang::abort(sprintf(
        "input object types %s not supported", errorobj))
    }
  } else {
    txtype <- is(exons)[1]
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

  # run testNMD_ for single GRanges object and output results
  if (intype == "gr") {
    return(testNMD(exons, cds, distance_stop_EJ = NMD_threshold))
  }

  # for GRangesList,
  if (intype == "grl") {
    totest <- names(exons) # prepare vector with names for testing
    if (!is.null(which)) {
      which_matched <- which[which %in% names(cds)]
      if (length(which_matched) == 0) {
        rlang::abort(sprintf("transcript names in `which` is not found in `%s`",
                             argnames[2]))
      } else if (length(which_matched) != length(which)) {
        num_unmatched <- length(which) - length(which_matched)
        rlang::warn(sprintf(
          "%s transcripts in `which` is missing from `%s%",
          num_unmatched, argnames[2]
        ))
      }

      totest <- totest[totest %in% which] # subset list if which list is given
      exons <- exons[names(exons) %in% which]
    }
    # check for missing cds and return warnings/errors
    totest <- totest[totest %in% names(cds)]
    if (length(totest) == 0) {
      rlang::abort("all transcripts have missing cds info. please ensure exons and cds names match")
    }
    if (length(totest) < length(exons)) {
      skiptest <- length(exons) - length(totest)
      rlang::warn(sprintf(
        "%s transcript(s) have missing cds info and have been skipped",
        skiptest
      ))
    }

    # running brlapply and testNMD
    out <- BiocParallel::bplapply(totest, function(x) {
      report <- list(
        transcript = x, is_NMD = F, dist_to_lastEJ = 0,
        num_of_down_EJs = 0, dist_to_downEJs = 0
      )
      NMDreport <- testNMD(exons[[x]], cds[[x]], distance_stop_EJ = NMD_threshold)
      report <- utils::modifyList(report, NMDreport)
      return(report)
    }, BPPARAM = BiocParallel::MulticoreParam()) %>%
      dplyr::bind_rows() %>%
      dplyr::filter(if (return[1] == "NMD") is_NMD == T else T)
    return(out)
  }
}

testNMD <- function(queryTranscript, queryCDS, distance_stop_EJ = 50) {

  # prepare output list
  output <- list(
    is_NMD = as.logical(FALSE),
    dist_to_lastEJ = as.numeric(0),
    num_of_down_EJs = as.numeric(0),
    dist_to_downEJs = as.numeric(0),
    threeUTRlength = as.numeric(0)
  )

  # define global variable
  width <- disttolastEJ <- NULL

  # sort queryCDS, by exon order (just in case)
  strand <- as.character(BiocGenerics::strand(queryCDS))[1]
  queryCDS <- BiocGenerics::sort(queryCDS, decreasing = strand == "-")
  queryTranscript <- BiocGenerics::sort(queryTranscript, decreasing = strand == "-")

  queryCDS <- resizeTranscript(queryCDS, end = -3)
  
  # test if query is NMD sensitive
  #   get distance of last EJC from start of transcript
  #   disjoin will create a new GRanges that will separate the queryTranscript
  #   into discernable 5UTR, ORF and 3UTR
  #   we can then try to use the new GRanges to infer NMD susceptibility
  lengthtolastEJ <- sum(head(BiocGenerics::width(queryTranscript), -1))
  disjoint <- BiocGenerics::append(queryCDS, queryTranscript) %>%
    GenomicRanges::disjoin(with.revmap = T) %>%
    BiocGenerics::sort(decreasing = strand == "-") %>%
    as.data.frame() %>%
    dplyr::mutate(disttolastEJ = lengthtolastEJ - cumsum(width)) %>%
    dplyr::mutate(threeUTR = dplyr::lead(rev(cumsum(rev(width))), default = 0))

  # retrieve index of last ORF segment and determine if there are
  # more than 1 exons after stop codon
  stopcodonindex <- max(which(lengths(disjoint$revmap) == 2))
  output$dist_to_lastEJ <- disjoint[stopcodonindex, ]$disttolastEJ
  output$threeUTRlength <- disjoint[stopcodonindex, ]$threeUTR

  # report number of downstream EJ and its distance to PTC
  output$num_of_down_EJs <- nrow(disjoint) - stopcodonindex - 1
  downEJCdf <- disjoint %>%
    dplyr::filter(dplyr::row_number() >= stopcodonindex + 1) %>%
    dplyr::filter(disttolastEJ >= 0) %>%
    dplyr::mutate(disttoPTC = cumsum(width))
  output$dist_to_downEJs <- paste(downEJCdf$disttoPTC, collapse = ",")

  if (output$dist_to_lastEJ > distance_stop_EJ) {
    # annotated transcript as NMD if dist_to_lastEJ is NMD triggering
    output$is_NMD <- TRUE
  }
  # return output
  return(output)
}
