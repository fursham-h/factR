#' Test consistency in seqlevels across objects
#'
#' @param ... Two or more objects with seqlevels information
#'
#' @export
#' @author Fursham Hamid
#' @return Logical value as to whether all objects have consistent seqlevels
#'
#' @importFrom dplyr %>%
#' @importFrom IRanges %over%
#' @importFrom IRanges %within%
#' @importFrom data.table :=
#' @importFrom graphics frame
#' @importFrom methods as is new slot slotNames
#' @importFrom stats dist end start
#' @importFrom utils head tail
#' @importFrom rlang .data
#' 
#' @examples 
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING TOY DATASET
#' ## ---------------------------------------------------------------------
#' require(GenomicRanges)
#' 
#' ## Create toy GTF objects
#' gr1 <- GRanges("1", IRanges(start = c(1,101), width=c(20,20)), '+')
#' gr2 <- GRanges("chr1", IRanges(start = c(1,101), width=c(20,20)), '+')
#' 
#' ## Test for seqlevels consistency
#' has_consistentSeqlevels(gr1, gr2)
#' 
#' ## Input can be Biostrings object with seqlevels
#' x0 <- c("CTCACCAGTAT", "TGTCAGTCGA")
#' dna <- Biostrings::DNAStringSet(x0)
#' seqlevels(dna) <- c("chr1", "chr2")
#' 
#' ## Test for seqlevels consistency
#' has_consistentSeqlevels(gr1, dna)
#' has_consistentSeqlevels(gr2, dna)
#' 
has_consistentSeqlevels <- function(...) {
  dots <- list(...)
  argnames <- as.character(match.call())[-1]

  if (length(dots) < 2) {
    rlang::abort("Insufficient input")
  }

  consistent <- c()
  for (i in seq(1, length(dots) - 1)) {
    for (j in seq(i + 1, length(dots))) {
      if (identical(dots[[i]], dots[[j]])) {
        next
      }
      test <- all(GenomeInfoDb::seqlevels(dots[[i]]) %in% GenomeInfoDb::seqlevels(dots[[j]]))
      consistent <- c(consistent, test)
      if (!test) {
        rlang::warn(sprintf(
          "Try running: %s <- matchSeqLevels(%s, %s)", argnames[i], argnames[i], argnames[j]
        ))
        break
      }
    }
  }
  return(all(consistent))
}
