#' Reference-guided construction of CDS on GTF object
#'
#' @description
#'
#' `buildCDS()` is designed to construct CDS information on transcripts from input GTF
#' object. 
#' 
#' @details 
#' The `buildCDS()`function will first assemble a database of annotated ATG codons from
#' `ref` GTF object. Then for each transcripts in `query`, the upstream-most annotated ATG will 
#' serve as the transcript's translation start site and `buildCDS()` will search for an
#' in-frame termination codon.
#'
#' @param query
#' GRanges object containing query GTF data.
#' @param ref
#' GRanges object containing reference GTF data.
#' @param fasta
#' BSgenome or Biostrings object containing genomic sequence
#'
#' @return
#' GRanges object containing query exon entries and newly-constructed CDS
#' information
#' 
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' buildCDS(matched_query_gtf, ref_gtf, Mmusculus)
#' 
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
  return(c(query, builtCDS))
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
Try running: %s <- matchChromosomes(%s, %s)",
      argnames[1], argnames[3], argnames[1], argnames[1], argnames[3]
    ))
  }
  if (suppressWarnings(!has_consistentSeqlevels(ref, fasta))) {
    rlang::abort(sprintf(
      "`%s` and `%s` has unmatched seqlevel styles. 
Try running: %s <- matchChromosomes(%s, %s)",
      argnames[2], argnames[3], argnames[2], argnames[2], argnames[3]
    ))
  }
}

.runbuildCDS <- function(query, ref, fasta, argnames) {
  transcript_id <- gene_id <- gene_name <- type <- width <- strand <- NULL
  phase <- tailphase <- coverage <- group <- seanames <- seqnames <- NULL
  group_name <- NULL

  # Add gene_name column if absent
  if (!"gene_name" %in% names(S4Vectors::mcols(query))) {
    query$gene_name <- query$gene_id
  }

  # extract important information
  totaltx <- query$transcript_id %>%
    unique() %>%
    length()
  genelist <- query %>%
    as.data.frame() %>%
    dplyr::select(transcript_id, gene_id, gene_name) %>%
    dplyr::distinct(transcript_id, .keep_all = T)
  query <- query[BiocGenerics::strand(query) != "*"]

  # Create lists of cds and exons
  query_exons <- S4Vectors::split(query[query$type == "exon"], ~transcript_id)
  ref_cds <- S4Vectors::split(ref[ref$type == "CDS"], ~transcript_id)
  ref_exons <- S4Vectors::split(ref[ref$type == "exon"], ~transcript_id)

  # run pairing of query to reference and unpack object
  q2r_output <- .pairq2r(query_exons, ref_cds, ref_exons, ref, fasta)
  q2r <- q2r_output[[1]]
  codons_gr <- q2r_output[[2]]

  # split q2r pairs into full covs and the rest
  fullcovs <- q2r %>%
    dplyr::filter(coverage == 1)
  q2r <- dplyr::setdiff(q2r, fullcovs)

  # prepare outputCDS for full coverages
  if (nrow(fullcovs) > 1) {
    fulloutCDS <- ref_cds %>%
      as.data.frame() %>%
      dplyr::filter(group_name %in% fullcovs[[3]]) %>%
      dplyr::select(group:strand) %>%
      dplyr::mutate(type = "CDS") %>%
      dplyr::left_join(fullcovs[, c(1, 3)],
        by = c("group_name" = "ref_transcript_id")
      ) %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
      dplyr::mutate(phase = rev(cumsum(rev(width) %% 3) %% 3)) %>%
      dplyr::ungroup() %>%
      dplyr::select(seqnames:phase)
  } else {
    fulloutCDS <- NULL
  }

  # prepare vector for remaining comparisons
  order_query <- query_exons[q2r$transcript_id]
  order_ref <- codons_gr[q2r$ref_transcript_index]

  # run .getCDS function for remaining transcripts
  if (length(order_query) > 0) {
    restoutCDS <- .getCDS(order_query, order_ref, fasta)
  } else {
    restoutCDS <- NULL
  }

  # # get number of newly-found CDSs
  # if (!is.null(restoutCDS)) {
  #   newCDS <- length(unique(restoutCDS$transcript_id))
  # } else {
  #   newCDS <- "none"
  # }

  # combine all CDSs and print out stats
  outCDS <- suppressWarnings(dplyr::bind_rows(fulloutCDS, restoutCDS))
  if (nrow(outCDS) > 0) {
    outCDS <- outCDS %>%
      dplyr::arrange(transcript_id, ifelse(strand == "-", dplyr::desc(start), start)) %>%
      dplyr::left_join(genelist, by = "transcript_id")

    successtx <- length(unique(outCDS$transcript_id))
    message(sprintf(
      "Out of %s transcripts in `%s`, %s transcript CDSs were built",
      totaltx, argnames[1], successtx
    ))
    return(GenomicRanges::makeGRangesFromDataFrame(outCDS, keep.extra.columns = T))
  } else {
    message(sprintf(
      "Out of %s transcripts in `%s`, none transcript CDSs were built",
      totaltx, argnames[1]
    ))
    return(NULL)
  }
}

.pairq2r <- function(query_exons, ref_cds, ref_exons, ref, fasta) {
  ref_transcript_id <- strand <- transcript_id <- coverage <- NULL
  ref_transcript_index <- type <- width <- phase <- tailphase <- NULL
  group <- resize <- width <- NULL

  # search for exact query and ref matches
  fulloverlap <- GenomicRanges::findOverlaps(query_exons, ref_exons,
    type = "equal", select = "first"
  )

  # prepare dataframe of query to reference pairs
  q2r <- data.frame(
    "transcript_id" = names(query_exons),
    "ref_transcript_index" = fulloverlap,
    stringsAsFactors = F
  ) %>%
    dplyr::filter(!is.na(ref_transcript_index)) %>%
    dplyr::mutate(ref_transcript_id = names(ref_exons[ref_transcript_index])) %>%
    dplyr::filter(ref_transcript_id %in% names(ref_cds)) %>%
    dplyr::mutate(coverage = 1)

  # search for first ATG overlap for non-exact transcripts
  if (NA %in% fulloverlap) {
    nonexact <- query_exons[is.na(fulloverlap)]
    subsetRef <- IRanges::subsetByOverlaps(ref, nonexact)

    # prepare a GRanges of trimmed CDS segments
    codons_gr <- subsetRef %>%
      as.data.frame() %>%
      dplyr::filter(type == "CDS") %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(ifelse(strand == "-", dplyr::desc(start), start)) %>%
      dplyr::mutate(tailphase = (cumsum(width) %% 3) %% 3) %>%
      dplyr::mutate(order = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(start = ifelse(strand == "+", start + phase, start + tailphase)) %>%
      dplyr::mutate(end = ifelse(strand == "-", end - phase, end - tailphase)) %>%
      dplyr::filter(end > start) %>%
      dplyr::group_by(strand) %>%
      dplyr::arrange(order, ifelse(strand == "-", dplyr::desc(end), start)) %>%
      GenomicRanges::makeGRangesFromDataFrame() %>%
      unique()
    codons_gr <- IRanges::subsetByOverlaps(codons_gr, nonexact, type = "within")

    # further trim CDS segments to the
    codons_seq <- BSgenome::getSeq(fasta, codons_gr)
    startMatch <- Biostrings::vmatchPattern("ATG", codons_seq) %>%
      as.data.frame() %>%
      dplyr::group_by(group) %>%
      dplyr::filter(end %% 3 == 0) %>%
      dplyr::mutate(start = start - 1)
    codons_gr <- codons_gr[startMatch$group]
    codons_gr <- GenomicRanges::resize(codons_gr, width = BiocGenerics::width(codons_gr) - startMatch$start, fix = "end")
    codons_gr <- GenomicRanges::resize(codons_gr, width = 3, fix = "start")

    firstATGoverlap <- GenomicRanges::findOverlaps(nonexact, codons_gr, minoverlap = 3, select = "first")
    firstATGq2r <- data.frame(
      "transcript_id" = names(nonexact),
      "ref_transcript_index" = firstATGoverlap,
      stringsAsFactors = F
    ) %>%
      dplyr::filter(!is.na(ref_transcript_index)) %>%
      dplyr::mutate(coverage = 2)
    q2r <- dplyr::bind_rows(q2r, firstATGq2r)
  } else {
    codons_gr <- NULL
  }

  return(list(q2r, codons_gr))
}

.getCDS <- function(order_query, order_ref, fasta) {
  strand <- start1 <- end1 <- group <- seqnames <- newstart <- newend <- NULL
  stoppos <- width <- seqnames <- strand <- type <- transcript_id <- phase <- NULL
  newstrand <- exonorder <- NULL
  strand...7 <- start...4 <- start...9 <- end...10 <- end...5 <- NULL
  group_name <- seqnames...3 <- NULL

  # get range from ATG to end of transcript
  startToend <- suppressMessages(dplyr::bind_cols(as.data.frame(range(order_query)), as.data.frame(order_ref)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(newstart = list(ifelse(strand...7 == "-", start...4, start...9))) %>%
    dplyr::mutate(newend = list(ifelse(strand...7 == "-", end...10, end...5))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(strand = as.character(strand...7)) %>%
    dplyr::mutate(newstrand = strand) %>%
    dplyr::mutate(newstrand = ifelse(strand == "-" & newend > end...5, "+", newstrand)) %>%
    dplyr::mutate(newstrand = ifelse(strand != "-" & newstart < start...4, "-", newstrand)) %>%
    dplyr::select(group, group_name, seqnames = seqnames...3, start = newstart, end = newend, strand = newstrand) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T))

  # get exon coordinates and sequence from ATG to end of transcript
  startToendexons <- GenomicRanges::pintersect(order_query, startToend, drop.nohit.ranges = T) %>%
    sorteach(exonorder)
  seq <- GenomicFeatures::extractTranscriptSeqs(fasta, startToendexons)

  # search for in-frame stop codon and return a df of its position in startToendexons
  CDSstop <- lapply(c("TAG", "TGA", "TAA"), function(x) {
    a <- Biostrings::vmatchPattern(x, seq)
    as.data.frame(a)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(group) %>%
    dplyr::filter(end %% 3 == 0) %>%
    dplyr::arrange(start) %>%
    dplyr::select(group, start) %>%
    dplyr::distinct(group, .keep_all = T) %>%
    dplyr::mutate(start = start - 1)

  # prepare df with 3UTR length (if any) for each transcript
  stopdf <- tibble::tibble(group = 1:length(seq), width = BiocGenerics::width(seq)) %>%
    dplyr::left_join(CDSstop, by = "group") %>%
    dplyr::mutate(threeUTR = ifelse(!is.na(start), (width - start), width))

  # resize 3' end of transcript to just before stop codon and prepare output
  outCDS <- trimTranscripts(startToendexons, start = 0, end = stopdf$threeUTR)
  outCDS <- outCDS %>%
    as.data.frame() %>%
    dplyr::mutate(type = "CDS") %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(phase = rev(cumsum(rev(width) %% 3) %% 3)) %>%
    dplyr::ungroup() %>%
    dplyr::select(seqnames, start, end, strand, type, transcript_id, phase)


  if (nrow(outCDS) == 0) {
    return(NULL)
  } else {
    return(outCDS)
  }
}
