identifyNMDexons <- function(x, NMD.result = NULL) {
  
  # checks etc, will add later
  ## check for x format, check for CDS
  ## check for NMD.result, if yes check if all tx have NMD annotation
  
  # get reference CDS transcript for each gene
  
  # get AS segments between NMD transcript and reference transcript
  
  # recreate hypothetical tx by inserting/removing AS segments
  #### difficult bit
  
  # test NMD again
  
}




#### Play with strand and chromsome ####
addExonstoTx <- function(x, y) {
  return(GenomicRanges::reduce(IRanges::punion(x,y)))
}


removeExonsfromTx <- function(x, y) {
  
  # combine x and y by its group, also annotate which pair could be an IR
  temp.y <- as.data.frame(y) %>% 
    dplyr::mutate(group=dplyr::row_number()) 
  
  x.y <- as.data.frame(x) %>% 
    dplyr::left_join(temp.y, by = "group") %>% 
    dplyr::mutate(is.internal = ifelse(start.x < start.y & end.x > end.y,
                                       TRUE, FALSE))
  
  # extract all IR events and change the end coord
  ## this will serve as a duplicate
  ir2 <- x.y[x.y$is.internal,] %>% 
    dplyr::mutate(end.y = end.x)
  
  # change start coord of IR events and combine with ir2
  x.y <- x.y %>% 
    dplyr::mutate(start.y = ifelse(is.internal, start.x, start.y)) %>% 
    dplyr::bind_rows(ir2)
  
  # Re-generate x and y GRanges
  x <- x.y %>%
    dplyr::select(group, group_name, ends_with("x")) %>%
    dplyr::rename_with(~ stringr::str_sub(.x,end=-3), ends_with("x")) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  y <- x.y %>%
    dplyr::select(group, ends_with("y")) %>%
    dplyr::rename_with(~ stringr::str_sub(.x,end=-3), ends_with("y")) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  
  # Remove exons
  diff <- GenomicRanges::psetdiff(x,y)
  diff$group_name <- x$group_name
  
  # prepare output grangeslist
  diff <- diff %>% as.data.frame() %>% 
    dplyr::mutate(group_name = x$group_name) %>% 
    dplyr::filter(width > 0) %>% 
    GenomicRanges::makeGRangesListFromDataFrame(split.field = "group_name") %>% 
    sorteach(exonorder)
  
  return(GenomicRanges::reduce(diff))
}