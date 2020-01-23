#' Match seqlevels of input GRanges to reference GRanges or BioString objects
#'
#' @description
#' A convenient wrapper to match seqlevels of a query GRanges object to a reference
#' object that contain seqlevels information. Reference can be a GRanges, GRangesList,
#' BioString or DNAString object. Seqlevels which fail to match will be dropped.
#'
#' @param x GRanges object with seqnames to change
#' @param to GRanges object from which seqnames is referenced
#'
#' @return Corrected input GRanges
#' @export
#' @author Fursham Hamid
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomeInfoDb mapSeqlevels
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb renameSeqlevels
#'
matchSeqLevels <- function(x, to) {
  suppressWarnings(
    if (any(!GenomeInfoDb::seqlevels(x) %in% GenomeInfoDb::seqlevels(to))){
      newStyle <- mapSeqlevels(seqlevels(x), (seqlevelsStyle(to)[1]))
      newStyle <- newStyle[!is.na(newStyle)]
      x <- renameSeqlevels(x, newStyle)

      if (any(!GenomeInfoDb::seqlevels(x) %in% GenomeInfoDb::seqlevels(to))) {
        GenomeInfoDb::seqlevels(x, pruning.mode = "tidy") <- as.vector(newStyle)
      }
    }
  )
  return(x)
}
