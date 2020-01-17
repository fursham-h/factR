
filtereach <- function(x, ...) {
  stopifnot(is(x, "GRangesList"))
  return(x %>% as.data.frame() %>%
    dplyr::filter(...) %>%
    dplyr::select(-group) %>%
    GenomicRanges::makeGRangesListFromDataFrame(split.field = 'group_name', 
                                                keep.extra.columns = T))
}
 # mutateeach?