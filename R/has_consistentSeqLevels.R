#' Test consistency of chromosome naming styles 
#' (aka seqlevels; e.g. "chr1" vs "1") across multiple objects
#'
#' @description
#' This function will determine if all input ranges objects have the 
#' same chromosome naming convention. Input objects can be GenomicRanges, 
#' BSgenome or Biostrings object with seqlevel information.
#'
#' @param ... Two or more objects with seqlevels information
#' @param verbose Whether to print out message
#'
#' @export
#' @author Fursham Hamid
#' @return Logical value as to whether all objects have consistent 
#' seqlevel styles
#'
#' @importFrom dplyr %>%
#' @importFrom IRanges %over%
#' @importFrom IRanges %within%
#' @importFrom data.table :=
#' @importFrom graphics frame
#' @importFrom methods as is new slot slotNames
#' @importFrom stats dist end start offset
#' @importFrom utils head tail
#' @importFrom rlang .data
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
#' ## Test for seqlevels consistency
#' has_consistentSeqlevels(gr1, gr2)
#'
#' ## Input can be a Biostrings object with seqlevels information
#' x0 <- c("chr2" = "CTCACCAGTAT", "chr3" = "TGTCAGTCGA")
#' dna <- Biostrings::DNAStringSet(x0)
#'
#' ## Test for seqlevels consistency
#' has_consistentSeqlevels(gr1, dna)
#' has_consistentSeqlevels(gr2, dna)
has_consistentSeqlevels <- function(..., verbose = TRUE) {
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
            test <- any(GenomeInfoDb::seqlevelsStyle(dots[[i]]) %in% 
                            GenomeInfoDb::seqlevelsStyle(dots[[j]]))
            testlevels <- sum(!GenomeInfoDb::seqlevels(dots[[i]]) %in% 
                                  GenomeInfoDb::seqlevels(dots[[j]]))
            consistent <- c(consistent, test)
            if (!test) {
                if(verbose){
                    rlang::warn(sprintf(
                        "Try running: %s <- matchChromosomes(%s, %s)", 
                        argnames[i], argnames[i], argnames[j]
                    )) 
                }
                
            } else if (testlevels > 0) {
                if(verbose){
                    rlang::warn(sprintf(
                        "%s seqlevel(s) in `%s` are not found in `%s`", 
                        testlevels, argnames[i], argnames[j]
                    ))
                }
                
            }
            break
        }
    }
    return(all(consistent))
}
