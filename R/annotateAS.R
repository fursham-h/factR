#' Annotate alternative-spliced segments from a GTF annotation
#'
#' @param x
#' GRanges object containing transcript features in GTF format
#'
#' @param as.data.frame
#' If TRUE, function will return a dataframe instead of GRanges GTF object
#' (default = FALSE)
#'
#' @param append
#' If TRUE, function will append the alternative-spliced segments to x and return
#' a new GTF GRanges object
#'
#' @return
#' GRanges object or data-frame containing alternative-spliced segments found in x
#'
#' @export
#' @author Fursham Hamid

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

  seqnames <- gene_id <- transcript_id <- coord <- AStype <- NULL

  # retrieve input object names and carry out checks
  argnames <- as.character(match.call())[-1]
  .ASchecks(x, argnames)

  # run Alternative splicing annotation function
  annotatedAS <- .runAS(x[x$type == "exon"])

  # return
  if (as.data.frame) {
    return(annotatedAS %>% as.data.frame() %>%
      dplyr::mutate(coord = paste0(seqnames, ":", start, "-", end)) %>%
      dplyr::select(gene_id, transcript_id, coord, AStype))
  } else if (append) {
    return(c(x, annotatedAS))
  } else {
    return(annotatedAS)
  }
}

#' Compare and classify alternative spliced segments between two or more transcripts
#'
#' @param x
#' Can be a GRangesList object containing exons for each transcripts (Intralist mode)
#'
#' Can be a GRanges object containing exons for a transcript (Pair-wise mode). If so,
#' at least one GRanges object is to be provided in `...` for comparison
#'
#' @param ...
#' In pair-wise mode, argument is one or more GRanges object containing exons for a
#' particular transcript to compare with `x`

#'
#' @return
#' GRangesList object containing alternatively-spliced segments for each transcript
#' in comparison.
#'
#' @author Fursham Hamid
#' @export
compareAS <- function(x, ...) {

  # catch missing args
  mandargs <- c("x")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    rlang::abort(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  # define global variables
  group_name <- transcript_id <- group <- NULL

  argnames <- as.character(match.call())[-1]

  if (is(x, "GRanges")) {
    if (!"transcript_id" %in% names(S4Vectors::mcols(x))) {
      x$transcript_id <- "transcript0"
    }
    x <- S4Vectors::split(x, ~transcript_id)
  }
  if (!is(x, "GRangesList")) {
    rlang::abort(sprintf("%s is not a GRanges or GRangesList object", argnames[1]))
  }

  # if a second GRanges object is given, carry out pairwise comparison
  if (length(list(...)) > 0) {
    # multiple comparisons not supported. can be a feature in the future
    # if (length(list(...)) > 1) {
    #   # idea for multiple comparison. convert grangeslist to df for comparison
    #   rlang::warn("Multiple comparisons is not yet supported in pair-wise mode. First item in ... used")
    # }
    dots <- list(...)
    newdots <- dots[unlist(lapply(dots, function(x) {
      is(x, "GRanges")
    }))]

    if (length(newdots) == 0) {
      # return warning and proceed to intra-list mode
      rlang::warn("No GRanges object found in `...`. Comparing AS in intra-list mode")
    } else {
      if (length(newdots) < length(dots)) {
        # return warning for elements which are not GRanges
        notGR <- length(dots) - length(newdots)
        rlang::warn("Non GRanges object in `...` were removed")
      }
      newdots <- as(newdots, "GRangesList")

      newdotsmeta <- newdots %>% as.data.frame()
      if (!"transcript_id" %in% names(newdotsmeta)) {
        names(newdots) <- paste0("transcript", as.character(c(1:length(newdots))))
        newdots <- mutateeach(newdots, transcript_id = group_name)
      } else {
        newdots <- newdots %>%
          as.data.frame() %>%
          dplyr::mutate(transcript_id = ifelse(is.na(transcript_id), paste0("transcript", group), transcript_id)) %>%
          dplyr::mutate(group_name = transcript_id) %>%
          dplyr::select(-group) %>%
          GenomicRanges::makeGRangesListFromDataFrame(split.field = "group_name", keep.extra.columns = T)
      }
    }
    x <- c(x, newdots)
    if (length(x) == 1) {
      rlang::abort("Insufficient transcripts for comparison")
    }
  }
  if (!"gene_id" %in% names(S4Vectors::mcols(unlist(x)))) {
    x <- mutateeach(x, gene_id = "NA")
  }
  return(S4Vectors::split(.runAS(x), ~transcript_id))
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
  transcript_id <- pos <- seqnames <- strand <- gene_id <- NULL
  first.X.start <- second.start <- first.X.end <- second.end <- NULL
  first.pos <- AStype <- first.X.strand <- NULL

  x <- x %>%
    as.data.frame() %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(pos = dplyr::row_number()) %>%
    dplyr::mutate(pos = ifelse(pos == 1, "First", pos)) %>%
    dplyr::mutate(pos = ifelse(pos == dplyr::n(), "Last", pos)) %>%
    dplyr::select(seqnames, start, end, strand, gene_id, transcript_id, pos) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

  # get reduced intron boundaries
  exonsbytx <- S4Vectors::split(x, ~transcript_id)
  intronsbytx <- GenomicRanges::psetdiff(unlist(range(exonsbytx)), exonsbytx)
  intronsreduced <- GenomicRanges::reduce(unlist(GenomicRanges::psetdiff(unlist(range(exonsbytx)), exonsbytx)))

  # get exons that overlap with reduced introns
  altexons <- IRanges::findOverlapPairs(x, intronsreduced)

  # annotate AS exons
  altannotate <- altexons %>%
    as.data.frame() %>%
    dplyr::mutate(AStype = "CE") %>%
    dplyr::mutate(AStype = ifelse(first.X.start < second.start & first.X.end < second.end & !first.pos %in% c("First", "Last"), "SD", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start > second.start & first.X.end > second.end & !first.pos %in% c("First", "Last"), "SA", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start < second.start & first.X.end < second.end & first.pos %in% c("First", "Last"), "Te", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start > second.start & first.X.end > second.end & first.pos %in% c("First", "Last"), "Ts", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start < second.start & first.X.end > second.end, "RI", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start > second.start & first.X.end < second.end & first.pos == "First", "FE", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start > second.start & first.X.end < second.end & first.pos == "Last", "LE", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.strand == "-", chartr("DAFLes", "ADLFse", AStype), AStype))

  # get segments that fall within intron and annotate that segment
  altexons <- GenomicRanges::pintersect(altexons)
  if (length(altexons) == 0) {
    altexons$pos <- altexons$hit <- NULL
    return(altexons)
  } else {
    altexons$type <- "AS"
    altexons$AStype <- toupper(altannotate$AStype)
    altexons$pos <- altexons$hit <- NULL

    return(altexons)
  }
}
