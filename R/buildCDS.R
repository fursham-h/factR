#' Construct query CDS using reference as guide
#'
#' @description
#'
#' buildCDS will select the #########
#'
#' buildCDS will next attempt to construct query CDS information by
#' firstly deriving its start_codon as guided by the reference transcript.
#' This is done in three successive stages:
#' (1) search transcript for a start_codon as annotated in reference
#' (2) if above fails, search for an internal ATG codon as annotated in reference
#' (3) if above fails, align the frame of query to reference
#'
#' After the start_codon/coding frame have been established, the program
#' will search for an in-frame stop_codon for each transcript and return its
#' CDS GRanges entries if successfull.
#'
#' @param query
#' GRanges object containing query GTF data.
#' @param ref
#' GRanges object containing reference GTF data.
#' @param fasta
#' BSgenome or Biostrings object containing genomic sequence
#'
#' @return
#' GRanges object containing query exon entries and newly-constructed cds
#' information
#' @export
#' @author Fursham Hamid
#'
buildCDS <- function(query, ref, fasta) {

  # catch missing args
  mandargs <- c("query", "ref", "fasta")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    rlang::abort(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }

  # retrieve input object names
  argnames <- as.character(match.call())[-1]

  # carry out input checks
  .checkObjects(query, ref, fasta, argnames)
  
  # run core building function
  builtCDS <- .runbuildCDS(query, ref, fasta, argnames)
  
  # return appended CDS info
  return(c(query,builtCDS))
}

.checkObjects <- function(query, ref, fasta, argnames) {
  # check inputs are gtf
  if (any(!is_gtf(query, ref))) {
    rlang::abort(sprintf(
      "%s is/are not gtf GRanges",
      paste(argnames[1:2][!is_gtf(query, ref)], collapse = ",")
    ))
  }
  
  # check if ref have CDS info
  if (!"CDS" %in% S4Vectors::mcols(ref)$type) {
    rlang::abort(sprintf(
      "`%s` have missing CDS info", argnames[2]
    ))
  }
  
  # catch unmatched seqlevels
  if (suppressWarnings(!has_consistentSeqlevels(query, fasta))) {
    rlang::abort(sprintf(
      "`%s` and `%s` has unmatched seqlevel styles. 
Try running: %s <- matchSeqLevels(%s, %s)",
      argnames[1], argnames[3], argnames[1], argnames[1], argnames[3]
    ))
  }
  if (suppressWarnings(!has_consistentSeqlevels(ref, fasta))) {
    rlang::abort(sprintf(
      "`%s` and `%s` has unmatched seqlevel styles. 
Try running: %s <- matchSeqLevels(%s, %s)",
      argnames[2], argnames[3], argnames[2], argnames[2], argnames[3]
    ))
  }
}

.runbuildCDS <- function(query, ref, fasta, argnames) {
  
  transcript_id <- gene_id <- gene_name <- type <- width <- strand <- NULL
  phase <- tailphase <- coverage <- group <- seanames <- seqnames <- NULL
  
  # Add gene_name column if absent
  if (!"gene_name" %in% names(S4Vectors::mcols(query))) {
    query$gene_name <- query$gene_id
  }
  
  # extract important information
  totaltx <- query$transcript_id %>% unique() %>% length()
  genelist <- query %>% as.data.frame() %>%
    dplyr::select(transcript_id, gene_id, gene_name) %>%
    dplyr::distinct()
  
  # Create lists of cds and exons
  query_exons <- S4Vectors::split(query[query$type == "exon"], ~transcript_id)
  ref_cds <- S4Vectors::split(ref[ref$type == "CDS"], ~transcript_id)
  ref_exons <- S4Vectors::split(ref[ref$type == "exon"], ~transcript_id)
  nc_ref_exons <- ref_exons[!names(ref_exons) %in% names(ref_cds)]
  ref_exons <- ref_exons[names(ref_exons) %in% names(ref_cds)]
  
  # create a list of CDS blocks
  codons_gr <- ref %>% as.data.frame() %>% 
    dplyr::filter(type == "CDS") %>% 
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(tailphase = (cumsum(width) %% 3) %% 3) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(start = ifelse(strand == "+", start + phase, start + tailphase)) %>% 
    dplyr::mutate(end = ifelse(strand == "-", end - phase, end- tailphase)) %>% 
    dplyr::filter(end > start) %>%
    dplyr::arrange(start) %>%
    GenomicRanges::makeGRangesFromDataFrame() %>% unique()
  codons_gr <- IRanges::subsetByOverlaps(codons_gr, query_exons, type = "within")
  
  # run pairing of query to reference
  q2r <- .pairq2r(query_exons, ref_exons, nc_ref_exons, codons_gr)
  
  # split q2r pairs into full covs and the rest
  fullcovs <- q2r %>%
    dplyr::filter(coverage == 1)
  q2r <- dplyr::setdiff(q2r, fullcovs)
  
  # prepare outputCDS for full coverages
  if (nrow(fullcovs) > 1) {
    fulloutCDS <- ref_cds %>%
      as.data.frame() %>%
      dplyr::filter(group %in% fullcovs[[3]]) %>%
      dplyr::select(group:strand) %>%
      dplyr::mutate(type = "CDS") %>%
      dplyr::left_join(fullcovs[, c(1, 3)],
                       by = c("group" = "ref_transcript_id")
      ) %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
      dplyr::mutate(phase = rev(cumsum(rev(width) %% 3) %% 3)) %>%
      dplyr::ungroup() %>%
      dplyr::select(seqnames:phase)
  }
  
  # prepare vector for remaining comparisons
  order_query <- query_exons[q2r$transcript_index]
  order_ref <- codons_gr[q2r$ref_transcript_id]
  
  
  # run buildCDS function
  if (length(order_query) > 0){
    restoutCDS <- .getCDS(order_query, order_ref, fasta)
  } else {
    restoutCDS <- NULL
  }
  
  # combine all CDSs
  outCDS <- suppressWarnings(dplyr::bind_rows(fulloutCDS, restoutCDS) %>%
                               dplyr::arrange(transcript_id, ifelse(strand == "-", dplyr::desc(start), start)) %>%
                               dplyr::left_join(genelist, by = "transcript_id"))
  
  # print out stats and return appended GRanges GTF
  if (nrow(outCDS) > 0) {
    successtx <- outCDS$transcript_id %>% unique() %>% length()
  } else {
    successtx <- 0
  }
  message(sprintf(
    "Out of %s transcripts in `%s`, %s transcript CDSs were built",
    totaltx, argnames[1], successtx
  ))
  return(GenomicRanges::makeGRangesFromDataFrame(outCDS, keep.extra.columns = T))
}

.pairq2r <- function(query_exons, ref_exons, nc_ref_exons, codons_gr){
  
  ref_transcript_id <- strand <- transcript_id <- coverage <- NULL

  # obtain a vector of strands for each transcript
  strandList <- as.character(S4Vectors::runValue(BiocGenerics::strand(query_exons)))
  
  # prepare q2r
  fulloverlap <- GenomicRanges::findOverlaps(query_exons, ref_exons,
                                             type = "equal", select = "first")
  ncoverlap <- GenomicRanges::findOverlaps(query_exons, nc_ref_exons,
                                           type = "equal", select = "first")
  overlappos <- GenomicRanges::findOverlaps(query_exons, codons_gr, select = "first")
  overlapneg <- GenomicRanges::findOverlaps(query_exons, codons_gr, select = "last")

  fullq2r <- data.frame(
    "transcript_id" = names(query_exons),
    "transcript_index" = 1:length(query_exons),
    "ref_transcript_id" = fulloverlap,
    "strand" = strandList,
    stringsAsFactors = F
  ) %>%
    dplyr::filter(!is.na(ref_transcript_id)) %>%
    dplyr::mutate(coverage = 1)
  ncq2r <- data.frame(
    "transcript_id" = names(query_exons),
    "transcript_index" = 1:length(query_exons),
    "ref_transcript_id" = ncoverlap,
    "strand" = strandList,
    stringsAsFactors = F
  ) %>%
    dplyr::filter(!is.na(ref_transcript_id)) %>%
    dplyr::mutate(coverage = 0)
 
  startq2rpos <- data.frame(
    "transcript_id" = names(query_exons),
    "transcript_index" = 1:length(query_exons),
    "ref_transcript_id" = overlappos,
    "strand" = strandList,
    stringsAsFactors = F
  ) %>%
    dplyr::filter(!is.na(ref_transcript_id) & strand != "-") %>%
    dplyr::mutate(coverage = 2)
  startq2rneg <- data.frame(
    "transcript_id" = names(query_exons),
    "transcript_index" = 1:length(query_exons),
    "ref_transcript_id" = overlapneg,
    "strand" = strandList,
    stringsAsFactors = F
  ) %>%
    dplyr::filter(!is.na(ref_transcript_id) & strand == "-") %>%
    dplyr::mutate(coverage = 2)
  
  q2r <- dplyr::bind_rows(fullq2r, ncq2r, startq2rpos, startq2rneg) %>%
    dplyr::distinct(transcript_id, .keep_all = T) %>%
    dplyr::filter(coverage > 0)
  
  return(q2r)
}

.getCDS <- function(order_query, order_ref, fasta) {
  
  strand <- start1 <- end1 <- group <- seqnames <- newstart <- newend <- NULL
  stoppos <- width <- seqnames <- strand <- type <- transcript_id <- phase <- NULL
  
  startToend <- dplyr::bind_cols(as.data.frame(range(order_query)), as.data.frame(order_ref)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(newstart = ifelse(strand == "-", start, start1)) %>%
    dplyr::mutate(newend = ifelse(strand == "-", end1, end)) %>%
    dplyr::select(group:seqnames, start = newstart, end = newend, strand) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
  
  startToendexons <- GenomicRanges::pintersect(order_query, startToend, drop.nohit.ranges=T)
  seq <- GenomicFeatures::extractTranscriptSeqs(fasta, startToendexons)
  # prepare a dict of stop codons for pattern matching
  list_stopcodons <- Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
  pdict_stopcodons <- Biostrings::PDict(list_stopcodons)
  
  CDSstop <- BiocParallel::bplapply(seq, function(x){
    stopcodons <- Biostrings::matchPDict(pdict_stopcodons, x) %>%
      unlist() %>%
      as.data.frame() %>%
      dplyr::filter(end %% 3 == 0) %>%
      dplyr::arrange(start)
    
    # return if no in-frame stop codons are found
    if (nrow(stopcodons) == 0) {
      return(0)
    } else {
      
      # retrieve 3UTR length and update output file
      return(stopcodons[1,1])
    }
  },BPPARAM = BiocParallel::MulticoreParam())
  
  stopdf <- tibble::tibble(width = BiocGenerics::width(seq), stoppos = unlist(CDSstop)) %>%
    dplyr::mutate(threeUTR = ifelse(stoppos > 0, width - stoppos, NA))
  
  outCDS <- BiocParallel::bpmapply(function(x, y) {
    if(!is.na(y)) {
      x <- resizeTranscript(x, start = 0, end = y)
      x$type = "CDS"
      S4Vectors::mcols(x)$phase <- rev(cumsum(rev(BiocGenerics::width(x)) %% 3) %% 3)
      x <- x %>% as.data.frame() %>%
        dplyr::select(seqnames, start, end, strand, type,transcript_id, phase)
      return(x)
    } else {
      return(NULL)
    }
    
  }, startToendexons, stopdf$threeUTR,
  BPPARAM = BiocParallel::MulticoreParam(), SIMPLIFY = F
  ) %>% dplyr::bind_rows()
  
  
  if (nrow(outCDS) == 0) {
    return(NULL)
  } else {
    return(outCDS)
  }
}


