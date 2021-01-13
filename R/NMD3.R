identifyNMDexons <- function(x, fasta, NMD.result = NULL) {
  
  # checks etc, will add later
  ## check for x format, check for CDS
  ## check for NMD.result, if yes check if all tx have NMD annotation
  if(is.null(NMD.result)) {
    NMD.result <- suppressMessages(predictNMD(x, progress_bar = FALSE))
  }
  
  # get reference CDS transcript for each gene
  ref <- .getbestref(x, NMD.result)

  # get AS segments between NMD transcript and reference transcript
  ## shortlist NMD transcripts from GTF
  x.NMD <- x[x$transcript_id %in% NMD.result[NMD.result$is_NMD,]$transcript]
  
  ## get AS segments and annotate its splicing nature
  AS.exons <- labelSplicedSegment(c(x.NMD, ref))
  AS.exons$splice <- ifelse(AS.exons$transcript_id %in% ref$transcript_id,
                            "skipped", "spliced")
  
  # simplify df by removing redundant exons
  ## in cases where ref and NMD tx contain said exon, the skipped form will be retained
  AS.exons <- AS.exons %>%
    as.data.frame() %>% 
    dplyr::mutate(coord = paste0(seqnames, ":", start, "-", 
                                 end, "_", strand, "_", AStype)) %>% 
    dplyr::arrange(splice) %>%
    dplyr::select(gene_id, transcript_id, coord, splice) %>% 
    dplyr::distinct(coord, .keep_all = TRUE)
    
  # recreate hypothetical tx by inserting/removing AS segments
  ## create grl of ref trancripts
  ref.grl <- S4Vectors::split(ref[ref$type == "exon"], ~gene_id)
  mod.tx <- tibble::tibble(
    "seqnames" = as.character(),
    "start" = as.integer(),
    "end" = as.integer(),
    "strand" = as.character(),
    "coord" = as.character()
  )
  
  ## Add exons to ref transcripts
  if("spliced" %in% AS.exons$splice) {
    toadd <- dplyr::filter(AS.exons, splice == "spliced")
    addedtx <- ref.grl[toadd$gene_id]
    names(addedtx) <- toadd$coord
    exons.toadd <- toadd["coord"] %>%
      tidyr::separate(coord, c("seqnames", "start", "end", "strand", "AStype"),
                      sep = ":|-|_") %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    addedtx <- addExonstoTx(addedtx, exons.toadd) %>% 
      as.data.frame() %>% 
      dplyr::select(seqnames, start, end, strand, coord = group_name)
    mod.tx <- dplyr::bind_rows(mod.tx, addedtx)
  }
  
  
  ## Remove exons from ref transcripts
  if("skipped" %in% AS.exons$splice) {
    toremove <- dplyr::filter(AS.exons, splice == "skipped")
    removedtx <- ref.grl[toremove$gene_id]
    names(removedtx) <- toremove$coord
    exons.toremove <- toremove["coord"] %>%
      tidyr::separate(coord, c("seqnames", "start", "end", "strand", "AStype"),
                      sep = ":|-|_") %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    removedtx <- removeExonsfromTx(removedtx, exons.toremove) %>% 
      as.data.frame() %>% 
      dplyr::select(seqnames, start, end, strand, coord = group_name)
    mod.tx <- dplyr::bind_rows(mod.tx, removedtx)
  }
  
  
  # cleanup mod.tx by adding gene_id info and create GRanges
  mod.tx <- mod.tx %>% 
    dplyr::left_join(AS.exons %>% dplyr::select(coord, gene_id), by = "coord") %>% 
    dplyr::mutate(transcript_id = coord, type = "exon") %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  # test NMD again
  mod.tx <- suppressMessages(buildCDS(mod.tx, ref, fasta))
  mod.NMD <- suppressMessages(predictNMD(mod.tx, progress_bar = FALSE))
  
  # Report NMD exons
  AS.exons <- AS.exons %>% 
    dplyr::left_join(mod.NMD %>% dplyr::select(transcript, is_NMD), 
              by = c("coord"="transcript")) %>% 
    dplyr::mutate(NMDtype = ifelse(splice == "skipped", "ORF-maintain", "Poison")) %>% 
    dplyr::select(coord, NMDtype, is_NMD) %>% 
    tidyr::separate(coord, c("seqnames", "start", "end", "strand", "AStype"),
                    sep = ":|-|_") %>% 
    dplyr::filter(is_NMD) %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  return(AS.exons)
}

.getbestref <- function(x, NMD.result) {
  # get reference CDS transcript for each gene
  ## get sizes of all CDSs
  cds.sizes <- sum(BiocGenerics::width(S4Vectors::split(x[x$type == "CDS"], 
                                                        ~transcript_id)))
  ## shortlist non NMD transcripts
  NMD.pos <- NMD.result[!NMD.result$is_NMD,]
  
  ## select best reference
  cds.reference <- x %>% 
    as.data.frame() %>% 
    dplyr::select(gene_id, transcript_id) %>% 
    dplyr::left_join(as.data.frame(cds.sizes) %>% tibble::rownames_to_column("transcript_id"),
                     by = "transcript_id") %>% 
    dplyr::filter(transcript_id %in% NMD.pos$transcript) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::slice_max(order_by = ".", n = 1)
  
  return(x[x$transcript_id %in% cds.reference$transcript_id])
}



addExonstoTx <- function(x, y, drop.unmodified = FALSE, 
                         allow.external.exons = FALSE) {
  # check for length of x and y
  if (length(x) != length(y)){
    rlang::abort("x and y are not of the same length")
  }
  
  # retain pairs with same chr and strand
  chr.x <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(x)))
  strand.x <- as.character(S4Vectors::runValue(BiocGenerics::strand(x)))
  chr.y <- as.character(GenomeInfoDb::seqnames(y))
  strand.y <- as.character(BiocGenerics::strand(y))
  
  samechr.strand <-  chr.x == chr.y & strand.x == strand.y
  if(sum(samechr.strand) == 0) {
    rlang::abort("All x-y pairs are on different chromosome and/or strand")
  } else if(sum(samechr.strand) < length(x) & drop.unmodified) {
    rlang::warn(sprintf("%s x-y pairs are on different chromosome and/or strand and 
  these transcripts were not modified",
                         length(x)-sum(samechr.strand)))
  }
  
  # check if exons in y are within transcripts in x
  internal <- as.data.frame(IRanges::pintersect(range(x), y))$hit
  if(allow.external.exons) {
    external.addition <- sum(samechr.strand & !internal)
    if(external.addition > 0) {
      rlang::warn(sprintf("%s external exons were added", external.addition))
    }
  } else {
    samechr.strand <- samechr.strand & internal
  }
  
  if(drop.unmodified){
    return(GenomicRanges::reduce(IRanges::punion(x,y))[samechr.strand])
  } else {
    modified <- GenomicRanges::reduce(IRanges::punion(x,y))[samechr.strand] %>% 
      mutateeach(modified = TRUE)
    unmodified <- x[!samechr.strand] %>% mutateeach(modified = FALSE)
    
    return(c(modified, unmodified)[names(x)])
  }
  
}



removeExonsfromTx <- function(x, y, drop.unmodified = FALSE,
                              ignore.strand = FALSE) {
  
  # checks for length and overlaps. returns nonoverlaps granges
  nonoverlaps.x <- .checklengthandoverlap(x, y , drop.unmodified, 
                                          ignore.strand)

  # merge x and y and prepare for intron removals
  x.y <- .mergexyandcorrectIR(x, y)
  # x.y$modified <- ifelse(x.y$group_name %in% names(nonoverlaps.x),
  #                        FALSE, TRUE)
  
  # Re-generate x and y GRanges
  x <- x.y %>%
    dplyr::select(group, group_name, ends_with("x")) %>%
    dplyr::rename_with(~ stringr::str_sub(.x,end=-3), ends_with("x")) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  y <- x.y %>%
    dplyr::select(group, ends_with("y")) %>%
    dplyr::rename_with(~ stringr::str_sub(.x,end=-3), ends_with("y")) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  
  # Core part: Removal of exon segments from transcripts in x
  diff <- suppressWarnings(
    GenomicRanges::psetdiff(x,y, ignore.strand = ignore.strand))
  diff$group_name <- x$group_name

  # prepare output grangeslist
  diff <- diff %>% as.data.frame() %>% 
    dplyr::mutate(group_name = x$group_name) %>% 
    dplyr::filter(width > 0) %>% 
    GenomicRanges::makeGRangesListFromDataFrame(split.field = "group_name",
                                                keep.extra.columns=TRUE) %>% 
    sorteach(exonorder) %>% 
    GenomicRanges::reduce() %>% 
    mutateeach(modified = ifelse(group_name %in% names(nonoverlaps.x),
                                 FALSE, TRUE))
  
  # drop non-overlapping pairs if requested, or annotate modified tx
  if(drop.unmodified) {
    diff <- filtereach(diff, modified) %>% 
      mutateeach(modified = NULL)
  }

  return(diff)
}

.checklengthandoverlap <- function(x, y, drop, ignore.strand) {
  # check for length of x and y
  if (length(x) != length(y)){
    rlang::abort("x and y are not of the same length")
  }
  
  # determine poverlaps between x and y and report
  temp.x <- IRanges::pintersect(x, y, ignore.strand = ignore.strand) %>% 
    filtereach(sum(hit) == 0)
  
  if(length(temp.x) == length(x)) {
    rlang::abort("All x-y pairs do not overlap")
  } else if(length(temp.x) > 0){
    if(drop) {
      rlang::warn(sprintf("%s x-y pairs do not overlap and these transcripts were not modified",
                          length(temp.x)))
    } 
  }
  return(temp.x)
}

.mergexyandcorrectIR <- function(x,y) {
  
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
  ## also, filter out entries of different chr/strand
  x.y <- x.y %>% 
    dplyr::mutate(start.y = ifelse(is.internal, start.x, start.y)) %>% 
    dplyr::bind_rows(ir2)
  return(x.y)
}