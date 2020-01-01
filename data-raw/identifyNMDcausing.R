identifyNMDcausing <- function(ASgranges, cds) {
  output <- list(
    NMDcausing = NA,
    NMDcausing.coord = NA
  )

  tx1index <- c(1:length(ASgranges))
  tx2index <- c((length(ASgranges) + 1):(length(ASgranges) + length(cds)))
  lastcodingexonstart <- BiocGenerics::start(cds[length(cds)])
  strand <- BiocGenerics::strand(ASgranges)[1] %>% as.character()

  combinedGRanges <- BiocGenerics::append(ASgranges, cds)
  disjoint <- GenomicRanges::disjoin(combinedGRanges, with.revmap = T)
  S4Vectors::mcols(disjoint)$AS <- IRanges::extractList(S4Vectors::mcols(combinedGRanges)$AS, disjoint$revmap) %>%
    as.list() %>%
    purrr::map(1) %>%
    unlist()

  disjoint <- disjoint %>%
    as.data.frame() %>%
    dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(as = ifelse(any(revmap %in% tx1index), T, F)) %>%
    dplyr::mutate(cds = ifelse(any(revmap %in% tx2index), T, F)) %>%
    dplyr::mutate(stop = ifelse(any(revmap %in% tail(tx2index, 1)), 1, 0)) %>%
    dplyr::mutate(outframe = ifelse(width %% 3 > 0, 1, 0)) %>%
    dplyr::mutate(dist = ifelse(strand == "-", start - lastcodingexonstart, lastcodingexonstart - start)) %>%
    dplyr::mutate(append = "") %>%
    as.data.frame()
  disjoint <- disjoint[(min(which(disjoint$cds == T))):nrow(disjoint), ]
  disjoint[((max(which(disjoint$cds == T))) + 1):nrow(disjoint), ]$append <- "3UTR-"
  disjoint$dist <- rank(replace(disjoint$dist, which(disjoint$dist < 0), NA), na.last = T)

  disjoint <- disjoint %>%
    dplyr::filter(as == T) %>%
    dplyr::arrange(dplyr::desc(stop), dplyr::desc(outframe), dist) %>%
    dplyr::mutate(AS = paste0(append, AS)) %>%
    dplyr::select(seqnames:strand, AS)

  output$NMDcausing <- disjoint[1, ]$AS
  output$NMDcausing.coord <- paste0(disjoint[1, 2:3], collapse = "-")

  return(output)
}
