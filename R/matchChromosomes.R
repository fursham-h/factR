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
#' @examples
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING TOY DATASET
#' ## ---------------------------------------------------------------------
#' require(GenomicRanges)
#'
#' ## Create toy GRanges objects
#' gr1 <- GRanges("1", IRanges(start = c(1, 101), width = c(20, 20)), "+")
#' gr2 <- GRanges("chr1", IRanges(start = c(1, 101), width = c(20, 20)), "+")
#'
#' ## Match Ensembl-style chromosomes from gr1 to UCSC-style gr2
#' matchChromosomes(gr1, gr2)
#'
#' ## Possible to match chrosomomes from GRanges object to a Biostrings object containing seqlevels
#' x0 <- c("chr2" = "CTCACCAGTAT", "chr3" = "TGTCAGTCGA")
#' dna <- Biostrings::DNAStringSet(x0)
#'
#' ## Match gr1 to dna
#' matchChromosomes(gr1, dna)
matchChromosomes <- function(x, to) {
    nseqlevelsbefore <- length(GenomeInfoDb::seqlevels(x))
    suppressWarnings(
        if (!has_consistentSeqlevels(x, to)) {
            # attempt to match style first
            newStyle <- mapSeqlevels(seqlevels(x), (seqlevelsStyle(to)[1]))
            newStyle <- newStyle[!is.na(newStyle)]
            x <- renameSeqlevels(x, newStyle)

            # # prune if there are remaining unmatched levels
            # if (any(!GenomeInfoDb::seqlevels(x) %in% GenomeInfoDb::seqlevels(to))) {
            #   GenomeInfoDb::seqlevels(x, pruning.mode = "tidy") <- as.vector(newStyle)
            # }
        }
    )
    if (length(GenomeInfoDb::seqlevels(x)) < nseqlevelsbefore) {
        nseqlevelsafter <- nseqlevelsbefore - length(GenomeInfoDb::seqlevels(x))
        rlang::warn(sprintf("%s seqlevels were removed after matching", nseqlevelsafter))
    }
    msg <- has_consistentSeqlevels(x, to)
    return(x)
}
