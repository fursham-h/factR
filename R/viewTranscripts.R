#' Plot transcripts directly from GTF.
#'
#' @description
#' A wrapper around wiggleplotr's plotTranscripts function. See the documentation for (\code{\link[wiggleplotr]{plotTranscripts}})
#' for more information.
#'
#' @param x
#' GRanges object containing transcript annotation in GTF format
#'
#' @param ...
#' Logical conditions to pass to dplyr::filter to subset transcripts for plotting.
#' Variables are metadata information found in `x` and multiple conditions can be
#' provided delimited by comma. Example: gene_name == "Ptbp1"
#'
#' @param rescale_introns
#' Specifies if the introns should be scaled to fixed length or not. (default: FALSE)
#'
#' @return ggplot2 object
#' @export
#'
#' @author Fursham Hamid
#'
#' @examples
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING SAMPLE DATASET
#' ## ---------------------------------------------------------------------
#' viewTranscripts(query_gtf)
#' viewTranscripts(query_gtf, transcript_id == "transcript1")
#' viewTranscripts(ref_gtf)
#'
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING TRANSCRIPT ANNOTATION
#' ## ---------------------------------------------------------------------
#' \donttest{
#' library(AnnotationHub)
#'
#' ## Retrieve GRCm38 trancript annotation
#' ah <- AnnotationHub()
#' GRCm38_gtf <- ah[["AH60127"]]
#'
#' ## Plot transcripts from Ptbp1 gene
#' viewTranscripts(GRCm38_gtf, gene_name == "Ptbp1")
#' }
#'
viewTranscripts <- function(x, ..., rescale_introns = FALSE) {

    # catch missing args
    mandargs <- c("x")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }

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
                    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
            },
            error = function(e) {
                rlang::abort(sprintf(
                    "Variables given in ... are not found in `%s`",
                    argnames[1]
                ))
            }
        )
        if (length(x) == 0) {
            rlang::abort("No transcripts to plot")
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

    # Control check for number of plotted transcripts
    if (length(exons) > 25) {
        exons <- exons[seq_len(25)]
        rlang::warn("Plotting only first 25 transcripts")
    }

    # main plot function
    plot <- suppressWarnings(wiggleplotr::plotTranscripts(
        exons = exons,
        cdss = cdss[names(cdss) %in% names(exons)],
        rescale_introns = rescale_introns
    ))
    plot
}
