sorteach <- function(x, ...) {
  stopifnot(is(x, "GRangesList"))
  expr <- quos(...)
  if ("~exonorder" %in% as.character(expr)) {
    index <- which('~exonorder' %in% as.character(expr))
    expr[[index]] <- quo(ifelse(strand == '-', dplyr::desc(start), start))
  }
  return(x %>% as.data.frame() %>%
           dplyr::arrange(!!!expr) %>%
           dplyr::select(-group) %>%
           GenomicRanges::makeGRangesListFromDataFrame(split.field = 'group_name', 
                                                       keep.extra.columns = T))
}