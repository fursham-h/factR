#' Label alternative-spliced segments from a GTF annotation
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
#' a new GTF GRanges object (default = TRUE)
#'
#' @return
#' GRanges object or data-frame containing alternative-spliced segments found in x
#'
#' @export
#' @author Fursham Hamid
#'
#' @examples
#'
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING TOY DATASET
#' ## ---------------------------------------------------------------------
#' require(GenomicRanges)
#'
#' ## Create toy GRanges GTF object
#' gr <- GRanges("chr1", IRanges(start = c(1, 101, 1, 51, 101), width = c(20, 20, 20, 10, 20)), "+",
#'   type = "exon", gene_id = "geneA",
#'   transcript_id = Rle(c("transcript1", "transcript2"), lengths = c(2, 3))
#' )
#'
#' ## Annotate alternative segments
#' labelSplicedSegment(gr)
#'
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING SAMPLE DATASET
#' ## ---------------------------------------------------------------------
#' ## Using GTF GRanges as input
#' labelSplicedSegment(query_gtf)
#'
#' ## Output as dataframe
#' labelSplicedSegment(query_gtf, as.data.frame = TRUE)
#'
#' ## Append AS info to input
#' labelSplicedSegment(query_gtf, append = TRUE)
#'
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING TRANSCRIPT ANNOTATION DATABASE
#' ## ---------------------------------------------------------------------
#' \dontrun{
#' library(AnnotationHub)
#'
#' ## Retrieve GRCm38 trancript annotation
#' ah <- AnnotationHub()
#' GRCm38_gtf <- ah[["AH60127"]]
#'
#' ## Run tool on specific gene family
#' labelSplicedSegment(GRCm38_gtf)
#' }
#'
labelSplicedSegment <- function(x, as.data.frame = FALSE, append = FALSE,
                                group.by = "gene_id") {
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

  # Add gene_name column if absent
  if (!"gene_name" %in% names(S4Vectors::mcols(x))) {
    x$gene_name <- x$gene_id
  }

  # run Alternative splicing annotation function
  annotatedAS <- .runAS(x[x$type == "exon"], group.by)

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
#'
#' @examples
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING TOY DATASET
#' ## ---------------------------------------------------------------------
#' require(GenomicRanges)
#'
#' ## Create toy GRanges GTF object
#' gr1 <- GRanges("chr1", IRanges(start = c(1, 101), width = c(20, 20)), "+")
#' gr2 <- GRanges("chr1", IRanges(start = c(1, 51, 101), width = c(20, 10, 20)), "+")
#'
#' ## Pairwise comparison between GRanges object
#' compareSplicedSegment(gr1, gr2)
#'
#' ## Multiple comparisons can be done by providing more GRanges input
#' gr3 <- GRanges("chr1", IRanges(start = c(1, 91), width = c(20, 30)), "+")
#' compareSplicedSegment(gr1, gr2, gr3)
#'
#' ## GRangesList containing exons per transcript can be given as input
#' grl <- GRangesList(list(gr1, gr2, gr3))
#' compareSplicedSegment(grl)
#'
#'
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING SAMPLE DATASET
#' ## ---------------------------------------------------------------------
#' ## Using GRangesList as input
#' compareSplicedSegment(query_exons)
#'
#' ## Compare AS between individual GRanges object
#' compareSplicedSegment(query_exons[[1]], query_exons[[3]])
compareSplicedSegment <- function(x, ...) {

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
  group_name <- transcript_id <- group <- gene_id <- NULL

  argnames <- as.character(match.call())[-1]

  if (is(x, "GRanges")) {
    x <- GenomicRanges::GRangesList(list(x))
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

      # newdotsmeta <- newdots %>% as.data.frame()
      # if (!"transcript_id" %in% names(newdotsmeta)) {
      #   names(newdots) <- paste0("transcript", as.character(c(1:length(newdots))))
      #   newdots <- mutateeach(newdots, transcript_id = group_name)
      # } else {
      #   newdots <- newdots %>%
      #     as.data.frame() %>%
      #     dplyr::mutate(transcript_id = ifelse(is.na(transcript_id), paste0("transcript", group), transcript_id)) %>%
      #     dplyr::mutate(group_name = transcript_id) %>%
      #     dplyr::select(-group) %>%
      #     GenomicRanges::makeGRangesListFromDataFrame(split.field = "group_name", keep.extra.columns = T)
      # }
    }
    x <- c(x, newdots)
  }

  # check metadata before running AS
  if (length(x) == 1) {
    rlang::abort("Insufficient transcripts for comparison")
  }
  if (!"transcript_id" %in% names(S4Vectors::mcols(unlist(x)))) {
    x <- mutateeach(x, group_name = group, transcript_id = NA)
  }
  x <- mutateeach(x,
    transcript_id = ifelse(is.na(transcript_id), group, transcript_id),
    group_name = transcript_id
  )
  if (!"gene_id" %in% names(S4Vectors::mcols(unlist(x)))) {
    x <- mutateeach(x, gene_id = "NA")
  }
  if (!"gene_name" %in% names(S4Vectors::mcols(x))) {
    x <- mutateeach(x, gene_name = gene_id)
  }

  return(S4Vectors::split(.runAS(x, "gene_id"), ~transcript_id))
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


.runAS <- function(x, group.by) {
  # define global variables
  transcript_id <- pos <- seqnames <- strand <- gene_id <- NULL
  first.X.start <- second.X.start <- first.X.end <- second.X.end <- NULL
  first.pos <- AStype <- first.X.strand <- gene_name <- termini <- NULL

  # order exons by chromosome coord and label position
  ##### mutate group name herr
  x <- x %>%
    as.data.frame() %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(pos = dplyr::row_number()) %>%
    dplyr::mutate(pos = ifelse(pos == 1, "First", pos)) %>%
    dplyr::mutate(pos = ifelse(pos == dplyr::n(), "Last", pos)) %>%
    dplyr::mutate(grouping = get(group.by)) %>%
    dplyr::select(seqnames, start, end, strand, gene_id, gene_name, transcript_id, pos, grouping) %>%
    dplyr::group_by(grouping) %>%  
    dplyr::arrange(start, end) %>%
    dplyr::mutate(termini = dplyr::row_number()) %>%
    dplyr::mutate(termini = ifelse(termini == 1 | termini == dplyr::n(), T, F)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  

  # get gene sizes
  t2g <- x %>%
    as.data.frame() %>%
    dplyr::select(grouping, transcript_id) %>%
    dplyr::distinct()

  # get reduced intron boundaries
  exonsbytx <- S4Vectors::split(x, ~transcript_id)
  intronsbytx <- GenomicRanges::psetdiff(BiocGenerics::unlist(range(exonsbytx)), exonsbytx)
  
  # regroup introns by groups
  intronsbygroup <- t2g %>% 
    dplyr::left_join(as.data.frame(intronsbytx), by = c("transcript_id" = "group_name")) %>% 
    dplyr::filter(!is.na(seqnames)) %>% 
    GenomicRanges::makeGRangesListFromDataFrame(split.field = "grouping") %>% 
    GenomicRanges::reduce() %>% 
    as.data.frame() %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  
  # intronsreduced <- GenomicRanges::reduce(unlist(GenomicRanges::psetdiff(unlist(range(exonsbytx)), exonsbytx)))
  # intronsreduced <- GenomicRanges::reduce(unlist(GenomicRanges::psetdiff(genewidths[t2g$gene_id], exonsbytx[t2g$transcript_id])))

  # get exons that overlap with reduced introns
  altexons <- IRanges::findOverlapPairs(x, intronsbygroup)

  # annotate internal AS exons
  altannotate <- altexons %>%
    as.data.frame() %>%
    dplyr::mutate(AStype = "CE") %>%
    dplyr::mutate(AStype = ifelse(first.X.start < second.X.start & first.X.end < second.X.end, "SD", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start > second.X.start & first.X.end > second.X.end, "SA", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start < second.X.start & first.X.end <= second.X.end & first.pos == "Last", NA, AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start >= second.X.start & first.X.end > second.X.end & first.pos == "First", NA, AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start < second.X.start & first.X.end > second.X.end, "RI", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start >= second.X.start & first.X.end < second.X.end & first.pos == "First", "FE", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.start > second.X.start & first.X.end <= second.X.end & first.pos == "Last", "LE", AStype)) %>%
    dplyr::mutate(AStype = ifelse(first.X.strand == "-", chartr("DAFLes", "ADLFse", AStype), AStype)) %>%
    dplyr::mutate(same.group = ifelse(first.grouping == second.group_name,TRUE,FALSE))

  # get segments that fall within intron and annotate that segment
  altexons <- GenomicRanges::pintersect(altexons)
  if (length(altexons) == 0) {
    altexons$pos <- altexons$hit <- altexons$termini <- altexons$grouping<- NULL
    return(altexons)
  } else {
    altexons$type <- "AS"
    altexons$AStype <- toupper(altannotate$AStype)
    altexons
    altexons <- altexons[!is.na(altexons$AStype) & altannotate$same.group]


    # get distal FE and LE coordinates
    distalFE <- x %>%
      as.data.frame() %>%
      dplyr::filter(gene_id %in% altannotate[altannotate$AStype == "FE", ]$first.gene_id, pos == ifelse(strand == "-", "Last", "First"), termini) %>%
      dplyr::mutate(type = "AS", AStype = "FE") %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

    distalLE <- x %>%
      as.data.frame() %>%
      dplyr::filter(gene_id %in% altannotate[altannotate$AStype == "LE", ]$first.gene_id, pos == ifelse(strand == "-", "First", "Last"), termini) %>%
      dplyr::mutate(type = "AS", AStype = "LE") %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)

    # append FE and LE coordinates and remove unwanted meta columns
    # altexons <- c(altexons, distalFE, distalLE)
    altexons$pos <- altexons$hit <- altexons$termini <-altexons$grouping <- NULL

    return(BiocGenerics::sort(altexons))
  }
}
