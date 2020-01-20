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
has_consistentSeqlevels <- function(...) {
  dots <- list(...)
  argnames <- as.character(match.call())[-1]
  
  if (length(dots) < 2) {
    rlang::abort("Insufficient input")
  }
  
  consistent <- c()
  for (i in seq(1, length(dots)-1)) {
    for (j in seq(i+1, length(dots))) {
      if (identical(dots[[i]], dots[[j]])) {
        next
      }
      test <- all(GenomeInfoDb::seqlevels(dots[[i]]) %in% GenomeInfoDb::seqlevels(dots[[j]]))
      consistent <- c(consistent, test)
      if (!test) {
        rlang::warn(sprintf(
          "Try running: %s <- matchSeqLevels(%s, %s)", argnames[i], argnames[i], argnames[j]))
        break
      }
      
    }
  }
  return(all(consistent))
}
  
