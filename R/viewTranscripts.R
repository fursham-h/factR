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
#' Character value of features to plot. Multiple features can be plotted by
#' entering comma-delimited values. Features will be extracted from
#' metadata gene_name, gene_id and transcript_id of the GTF.
#'
#' @param rescale_introns
#' Specifies if the introns should be scaled to fixed length or not. (default: FALSE)
#' 
#' @param ncol
#' Number of columns to patch the output plots (default: 1)
#'
#' @return ggplot2 object. If multiple genes are detected, plots will be 
#' combined using patchwork
#' @export
#'
#' @author Fursham Hamid
#'
#' @examples
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING SAMPLE DATASET
#' ## ---------------------------------------------------------------------
#' viewTranscripts(query_gtf)
#' viewTranscripts(query_gtf, "transcript1")
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
#' viewTranscripts(GRCm38_gtf, "Ptbp1")
#' 
#' # Plot transcripts from Ptbp1 and Ptbp2 genes
#' viewTranscripts(GRCm38_gtf, "Ptbp1", "Ptbp2")
#' }
#'
viewTranscripts <- function(x, ..., rescale_introns = FALSE, ncol = 1) {
    
    gene_name <- gene_id <- transcript_id <- NULL
    transcript_id <- meta <- val <- n <- type <-  NULL

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
    
    # prepare features
    featmeta <- tryCatch(
        {
            GenomicRanges::mcols(x) %>% 
                as.data.frame() %>% 
                dplyr::select(gene_name, gene_id, transcript_id) %>% 
                dplyr::mutate(n = dplyr::row_number()) %>% 
                tidyr::gather(meta, val, -n)
        },
        error = function(e) {
            GenomicRanges::mcols(x) %>% 
                as.data.frame() %>% 
                dplyr::select(-type) %>% 
                dplyr::mutate(n = dplyr::row_number()) %>% 
                tidyr::gather(meta, val, -n)
        }
    )

    if (!missing(...)) {
        x <- tryCatch(
            {
                x[featmeta[featmeta$val %in% c(...),"n"]]
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
    ngenes <- unique(x$gene_name)
    plot <- BiocGenerics::do.call(patchwork::wrap_plots, 
                                  lapply(ngenes, function(y){
                                      # Fetch gene exons and cdss
                                      exons <- S4Vectors::split(x[x$type == "exon" & x$gene_name == y], ~transcript_id)
                                      cdss <- S4Vectors::split(x[x$type == "CDS" & x$gene_name == y], ~transcript_id)
                                      as <- S4Vectors::split(x[x$type == "AS" & x$gene_name == y], ~transcript_id)
                                      if (length(cdss) == 0) {
                                          cdss <- NULL
                                      }
                                      
                                      
                                      # Control check for number of plotted transcripts
                                      if (length(exons) > 25) {
                                          exons <- exons[seq_len(25)]
                                          rlang::warn(sprintf("Plotting only first 25 transcripts for %s gene", y))
                                      }
                                      
                                      # main plot function
                                      suppressWarnings(wiggleplotr::plotTranscripts(
                                          exons = exons,
                                          cdss = cdss[names(cdss) %in% names(exons)],
                                          rescale_introns = rescale_introns
                                      )) + ggplot2::ggtitle(y)
                                  }))

    
    plot + patchwork::plot_layout(ncol = ncol)
}
