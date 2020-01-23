mutateeach <- function(x, ...) {
  group <- NULL
  stopifnot(is(x, "GRangesList"))
  return(x %>% as.data.frame() %>%
           dplyr::mutate(...) %>%
           dplyr::select(-group) %>%
           GenomicRanges::makeGRangesListFromDataFrame(
             split.field = "group_name",
             keep.extra.columns = T
           ))
}
