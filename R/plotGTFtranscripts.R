#' Plot transcripts directly from GTF.
#'
#' A wrapper around the plotTranscripts function. See the documentation for (\code{\link[wiggleplotr]{plotTranscripts}})
#' for more information.
#'
#' @param x
#' GRanges object containing transcript annotation in GTF format
#'
#' @param ...
#' Logical conditions to pass to dplyr::filter to subset transcripts for analysis.
#' Variables are metadata information found in `x` and multiple conditions can be
#' provided delimited by comma. Example: transcript_id == "transcript1"
#'
#' @param rescale_introns
#' Specifies if the introns should be scaled to fixed length or not. (default: FALSE)
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' plotGTFtranscripts(query_gtf)
#' plotGTFtranscripts(query_gtf, transcript_id == "transcript1")
#' plotGTFtranscripts(ref_gtf)
plotGTFtranscripts <- function(x, ..., rescale_introns = F) {

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
  # is_NMD <- NULL

  # retrieve input object names
  argnames <- as.character(match.call())[-1]

  if (!is_gtf(x)) {
    rlang::abort(sprintf("`%s` is not a GTF GRanges object", argnames[1]))
  }

  if (!missing(...)) {
    x <- tryCatch(
      {
        x %>%
          as.data.frame() %>%
          dplyr::filter(...) %>%
          GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
      },
      error = function(e) {
        rlang::abort(sprintf(
          "Metadata in ... not found in `%s`",
          argnames[1]
        ))
      }
    )
    if (length(x) == 0) {
      rlang::abort("No transcripts to display")
    }
  }

  # gene_metadata <- x %>% as.data.frame() %>%
  #   dplyr::select(transcript_id, gene_name, strand)
  #
  # Fetch gene exons and cdss
  exons <- S4Vectors::split(x[x$type == "exon"], ~transcript_id)
  cdss <- S4Vectors::split(x[x$type == "CDS"], ~transcript_id)
  if (length(cdss) == 0) {
    cdss <- NULL
  }

  wiggleplotr::plotTranscripts(
    exons = exons,
    cdss = cdss,
    rescale_introns = rescale_introns
  )
}
