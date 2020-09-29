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
#' provided delimited by comma. Example: gene_name == "Ptbp1"
#'
#' @param rescale_introns
#' Specifies if the introns should be scaled to fixed length or not. (default: FALSE)
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' viewTranscripts(query_gtf)
#' viewTranscripts(query_gtf, transcript_id == "transcript1")
#' viewTranscripts(ref_gtf)
viewTranscripts <- function(x, ..., rescale_introns = F) {

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

  # Need to have a check for plotting multiple genes.....

  # Fetch gene exons and cdss
  exons <- S4Vectors::split(x[x$type == "exon"], ~transcript_id)
  cdss <- S4Vectors::split(x[x$type == "CDS"], ~transcript_id)
  as <- S4Vectors::split(x[x$type == "AS"], ~transcript_id)
  if (length(cdss) == 0) {
    cdss <- NULL
  }

  plot <- wiggleplotr::plotTranscripts(
    exons = exons,
    cdss = cdss,
    rescale_introns = rescale_introns
  )
  plot
}
