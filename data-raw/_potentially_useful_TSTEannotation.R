
# terminal AS exon analysis
# termini <- c(trimTranscripts(exonsbytx, start =0, end = sum(width(exonsbytx))-1),
#              trimTranscripts(exonsbytx, start =sum(width(exonsbytx))-1, end = 0)) %>%
#   as.data.frame() %>%
#   dplyr::select(seqnames:pos) %>%
#   dplyr::distinct() %>%
#   GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
#
#
# disjointed_exons <- disjoin(x)
# disjointed_termini_exons <- mergeByOverlaps(disjointed_exons, termini)
# disjointed_internal_exons <- GenomicRanges::setdiff(disjointed_exons, disjointed_termini_exons$disjointed_exons)
#
# disjointed_termini_exons$fullexon <- overlapsAny(disjointed_termini_exons$disjointed_exons, x, type = "equal")
# disjointed_termini_exons$intronic <- overlapsAny(disjointed_termini_exons$disjointed_exons, c(intronsreduced,x[x$termini]), type = "within")
#
# disjointed_termini_exons <- disjointed_termini_exons %>%
#   as.data.frame() %>%
#   dplyr::select(seqnames = disjointed_exons.seqnames, start = disjointed_exons.start, end = disjointed_exons.end,
#                 strand = disjointed_exons.strand, gene_id:intronic) %>%
#   dplyr::mutate(end = ifelse((!fullexon | !intronic)  & pos == "First", start, end)) %>%
#   dplyr::mutate(start = ifelse((!fullexon | !intronic) & pos == "Last", end, start)) %>%
#   makeGRangesFromDataFrame(keep.extra.columns = T)
#
#
# disjointed_termini_exons$distal <- as.logical(overlapsAny(disjointed_termini_exons, x[x$termini], type = "start") + overlapsAny(disjointed_termini_exons, x[x$termini], type = "end"))
# disjointed_termini_exons$tandemts <- overlapsAny(disjointed_termini_exons, x[x$termini], type = "within")
# disjointed_termini_exons$exonic <- overlapsAny(disjointed_termini_exons, disjointed_internal_exons, type = "within")
#
# firstlastintron <- subsetByOverlaps(intronsreduced, x[x$termini], maxgap = 0)
# disjointed_termini_exons$tandemle <- overlapsAny(disjointed_termini_exons, firstlastintron, type = "within")
#
# termaltexons <- disjointed_termini_exons %>%
#   as.data.frame() %>%
#   mutate(type = "AS") %>%
#   mutate(AStype = ifelse(width > 1,"FE","Ts")) %>%
#   mutate(AStype = ifelse(pos=="Last" & AStype == "FE", "LE", AStype)) %>%
#   mutate(AStype = ifelse(pos=="Last" & AStype == "Ts", "Te", AStype)) %>%
#   mutate(prefix = "intronic") %>%
#   mutate(prefix = ifelse(distal, "distal", prefix)) %>%
#   mutate(prefix = ifelse(width == 1 & tandemts & !distal, "tandem", prefix)) %>%
#   mutate(prefix = ifelse(width == 1 & !tandemts & !distal & exonic, "exonic", prefix)) %>%
#   mutate(prefix = ifelse(width >1 & tandemle, "tandem", prefix)) %>%
#   mutate(prefix = ifelse(AStype %in% c("FE","LE") & prefix == "intronic", "internal", prefix)) %>%
#   dplyr::mutate(AStype = ifelse(strand == "-", chartr("FLse", "LFes", AStype), AStype)) %>%
#   dplyr::mutate(AStype = toupper(AStype)) %>%
#   dplyr::select(seqnames:transcript_id, type, AStype, AStype.ext = prefix) %>%
#   GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
