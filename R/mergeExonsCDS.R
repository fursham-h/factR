#' Merge GRangesList exons and cds objects
#'
#' @description
#' Merges two GRangesList objects; one containing exon coordinates and one containing cds coordinates.
#'
#' @param exons
#' GRangesList object containing exons
#' for each transcript.
#' @param cds
#' GRangesList object containing cds
#' for each transcript.
#'
#' @return
#' GRanges object
#' @export
mergeExonsCDS <- function(exons, cds) {

  # catch missing args
  mandargs <- c("exons", "cds")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }

  # check seqlevels
  if (GenomeInfoDb::seqlevelsStyle(exons) != GenomeInfoDb::seqlevelsStyle(cds)) {
    stop("exons and cds has unmatched seqlevel styles. try matching using matchSeqLevels function")
  }

  # unlist GRangesList and fill important attribute columns
  exons <- unlist(exons)
  cds <- unlist(cds)
  if (any(!c("type", "transcript_id") %in% names(S4Vectors::mcols(exons)))) {
    metadata <- c("type", "transcript_id")
    missing <- metadata[!metadata %in% names(S4Vectors::mcols(exons))]
    if ("type" %in% missing) {
      S4Vectors::mcols(exons)$type <- "exon"
    }
    if ("transcript_id" %in% missing) {
      S4Vectors::mcols(exons)$transcript_id <- names(exons)
    }
  }
  if (any(!c("type", "transcript_id", "phase") %in% names(S4Vectors::mcols(cds)))) {
    metadata <- c("type", "transcript_id")
    missing <- metadata[!metadata %in% names(S4Vectors::mcols(cds))]
    if ("type" %in% missing) {
      S4Vectors::mcols(cds)$type <- "CDS"
    }
    if ("transcript_id" %in% missing) {
      S4Vectors::mcols(cds)$transcript_id <- names(cds)
    }
    if ("phase" %in% missing) {
      S4Vectors::mcols(cds)$phase <- data.table::shift(cumsum(BiocGenerics::width(exons) %% 3) %% 3, fill = 0)
    }
  }
  names(exons) <- NULL
  names(cds) <- NULL

  return(c(exons, cds))
}
