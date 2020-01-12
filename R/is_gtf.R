is_gtf <- function(...) {
  return(unlist(lapply(list(...), function(x){
    if (is(x, "GRanges")){
      x <- as.data.frame(x)
      if(all(c('type', 'gene_id', 'transcript_id') %in% names(x))){
        x <- x %>%
          dplyr::select(type, gene_id, transcript_id) %>%
          dplyr::distinct()
        if (nrow(x) >1){
          return(TRUE)
        }
      }
    }
    return(FALSE)
  })))
}