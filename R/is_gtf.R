#' Test for GTF GRanges object
#'
#' @param ... One or more objects to test
#'
#' @return Logical values as to whether each object is a GTF GRanges
#' @author Fursham Hamid
#'
#' @examples
#' is_gtf(query_gtf, ref_gtf)
is_gtf <- function(...) {
  type <- gene_id <- transcript_id <- NULL
  return(unlist(lapply(list(...), function(x) {
    if (is(x, "GRanges")) {
      x <- as.data.frame(x)
      if (all(c("type", "gene_id", "transcript_id") %in% names(x))) {
        x <- x %>%
          dplyr::select(type, gene_id, transcript_id) %>%
          dplyr::distinct()
        if (nrow(x) > 1) {
          return(TRUE)
        }
      }
    }
    return(FALSE)
  })))
}
