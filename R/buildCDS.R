#' Construct query CDS using reference as guide
#'
#' @description 
#' 
#' 
#' buildCDS will firstly construct a dataframe containing query-reference 
#' transcript pairs. In situation where a reference gene expresses more than one
#' coding transcript, the program will select a query-reference transcript
#' pair with the highest coverage score.
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
#' GRanges object containing query GTF data. Gene_id metadata have to match
#' the gene_ids in `ref` in order for the program to create a query-reference 
#' transcript pairs. See ?matchGeneIDs to match query gene_id to a reference
#' annotation
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
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' buildCDS(matched_query_gtf, ref_gtf, Mmusculus)
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
  .buildCDSchecks(query, ref, argnames)

  # makeGRangesList
  query_exons <- S4Vectors::split(query[query$type == 'exon'], ~transcript_id)
  ref_cds <- S4Vectors::split(ref[ref$type == 'CDS'], ~transcript_id)
  ref_exons <- S4Vectors::split(ref[ref$type == 'exon'], ~transcript_id)
  ref_exons <- ref_exons[names(ref_exons) %in% names(ref_cds)] 
  
  # prepare q2r
  fulloverlap <- GenomicRanges::findOverlaps(query_exons, ref_exons, 
                                             type = 'equal', select = "first")
  startoverlap <- GenomicRanges::findOverlaps(query_exons, ref_exons, 
                                              type = 'start', select = "first")
  anyoverlap <- GenomicRanges::findOverlaps(query_exons, ref_exons, 
                                            type = 'any', select = "arbitrary")

  fullq2r <- data.frame('transcript_id' = names(query_exons),
                        'transcript_index' = 1:length(query_exons),
                    'ref_transcript_id' = fulloverlap, 
                    stringsAsFactors = F) %>%
    dplyr::filter(!is.na(ref_transcript_id)) %>%
    dplyr::mutate(coverage = 1)
  startq2r <- data.frame('transcript_id' = names(query_exons),
                         'transcript_index' = 1:length(query_exons),
                        'ref_transcript_id' = startoverlap, 
                        stringsAsFactors = F) %>%
    dplyr::filter(!is.na(ref_transcript_id)) %>%
    dplyr::filter(!transcript_id %in% fullq2r$transcript_id) %>%
    dplyr::mutate(coverage = 2)
  anyq2r <- data.frame('transcript_id' = names(query_exons),
                       'transcript_index' = 1:length(query_exons),
                        'ref_transcript_id' = anyoverlap, 
                       stringsAsFactors = F) %>%
    dplyr::filter(!is.na(ref_transcript_id)) %>%
    dplyr::filter(!transcript_id %in% fullq2r$transcript_id &
                  !transcript_id %in% startq2r$transcript_id) %>%
    dplyr::mutate(coverage = 3)
  
  q2r <- dplyr::bind_rows(fullq2r, startq2r, anyq2r)
  # report unmatched tx


  # run buildCDS function
   outCDS <- .getCDSgr(query_exons, ref_cds, fasta, q2r)
   if (is.null(outCDS)) {
     successtx <- 0
   } else {
     successtx <- length(unique(S4Vectors::mcols(outCDS)$transcript_id))
   }
   
   # print out stats and return appended GRanges GTF
   message(sprintf(
     "Out of %s transcripts in `%s`, %s transcript CDSs were built", 
     length(query_exons), argnames[1], successtx))
   return(c(query, outCDS))
}



.getCDSgr <- function(query, refCDS, fasta, query2ref) {

  # define global variables
  coverage <- group_name <- strand <- width <- transcript_id <- NULL

  # get total query tx and prepare output 
  totaltx <- nrow(query2ref)
  outCDS <- NULL

  # create CDS list for tx with coverage of 1
  fullcovs <- query2ref %>%
    dplyr::filter(coverage == 1)
  query2ref <- dplyr::setdiff(query2ref, fullcovs)
    
  # prepare outputCDS for full coverages
  if (nrow(fullcovs) > 1){
    outCDS <- refCDS %>%
      as.data.frame() %>%
      dplyr::filter(group_name %in% fullcovs[[2]]) %>%
      dplyr::select(group_name:strand) %>%
      dplyr::mutate(type = "CDS") %>%
      dplyr::left_join(fullcovs[1:2],
                       by = c("group_name" = "ref_transcript_id")
      ) %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start)) %>%
      dplyr::mutate(phase = rev(cumsum(rev(width) %% 3) %% 3)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-group_name) %>%
      dplyr::mutate(built_from = 'Full coverage')
  }
  


  # create CDS list for all remaining tx
  out <- BiocParallel::bpmapply(function(x, y, z) {
    CDSreport <- .getthisCDS(y, z, fasta) %>%
      as.data.frame()
    return(CDSreport)
  }, q2r$transcript_id, 
  query_exons[q2r$transcript_index], 
  ref_cds[q2r$ref_transcript_id],
  BPPARAM = BiocParallel::MulticoreParam(), SIMPLIFY = F
  ) %>%
    dplyr::bind_rows()
  outCDS <- suppressWarnings(dplyr::bind_rows(outCDS, out) %>%
    dplyr::arrange(transcript_id, ifelse(strand == '-', dplyr::desc(start), start)))
  
  if (length(unique(outCDS$transcript_id)) < totaltx) {
    unsuccesstx <- totaltx - length(unique(outCDS$transcript_id))
    rlang::warn(sprintf(
      "%s transcripts CDS were not found", unsuccesstx))
  }
  
  if (nrow(outCDS) ==  0) {
    return(NULL)
  } else {
    return(GenomicRanges::makeGRangesFromDataFrame(outCDS, keep.extra.columns = T))
  }
}

.getthisCDS <- function(query, CDS, fasta) {

  # prepare output list
  output <- list(
    ORF_considered = as.character(NA),
    ORF_start = as.character("Not found"),
    fiveUTRlength = 0,
    threeUTRlength = 0,
    ORF_found = FALSE
  )
  
  # Extract info and sort all GRanges first
  strand <- as.character(BiocGenerics::strand(query))[1]
  queryTx <- BiocGenerics::sort(query, decreasing = strand == "-")
  knownCDS <- BiocGenerics::sort(CDS, decreasing = strand == "-")
  S4Vectors::mcols(queryTx)$transcript_id <- names(query)
  # attempt to find an aligned start codon
  report <- .getCDSstart(queryTx, knownCDS, fasta)
  output <- utils::modifyList(output, report) # update output

  # return if no start codon is found
  if (output$ORF_start == "Not found") {
    return()
  }

  # attempt to search for an in-frame stop codon
  report <- .getCDSstop(queryTx, fasta, output$fiveUTRlength)
  output <- utils::modifyList(output, report) # update output

  # return if no stop codon is found
  if (output$ORF_found == FALSE) {
    return()
  }

  # build new ORF Granges
  report <- .getCDSranges(queryTx, output$fiveUTRlength, output$threeUTRlength,
                          output$ORF_start)
  #output <- utils::modifyList(output, report) # update output
  # return(output[c('ORF_considered', 'ORF_start', 'ORF_found')])
  return(report)
}

.getCDSstart <- function(query, refCDS, fasta) {

  # prepare output list
  output <- list(
    ORF_start = as.character("Not found"),
    fiveUTRlength = 0
  )
  # define global variables
  phase <- width <- NULL
  
  # return if query and ref do not overlap at all
  if (all(suppressWarnings(!refCDS %within% query))) {
    return(output)
  }
  
  # sort all GRanges, in case
  strand <- as.character(BiocGenerics::strand(query))[1]
  query <- BiocGenerics::sort(query, decreasing = strand == "-")
  refCDS <- BiocGenerics::sort(refCDS, decreasing = strand == "-")

  # get coord of start codon on reference and strand info
  startcodon <- resizeTranscript(refCDS, end = sum(BiocGenerics::width(refCDS)) - 3)
  
  # if query containg annotated start codon:
  if (startcodon %within% query & length(startcodon) == 1) {

    # this part of the code will calculate the length of 5'UTR
    #   disjoin will break query GRanges into sub GRanges, based
    #   on the position of the start codon.
    #   we can then find length of the upstream/downstream segments
    #   of the break
    disjoint <- BiocGenerics::append(query, startcodon) %>%
      GenomicRanges::disjoin(with.revmap = T) %>%
      sort(decreasing = strand == "-") %>%
      as.data.frame() %>%
      dplyr::mutate(cumsum = cumsum(width))


    # retrieve index of segment upstream of start codon and return its cumsumwidth
    startcodonindex <- min(which(lengths(disjoint$revmap) == 2))
    if (startcodonindex > 1) {
      fiveUTRlength <- disjoint[startcodonindex - 1, ]$cumsum
    } else {
      fiveUTRlength <- 0
    }

    # update output list
    output$ORF_start <- "Annotated ATG"
    output$fiveUTRlength <- fiveUTRlength

    return(output)
  }
  # if annotated start is not found, attempt to find upstream-most internal ATG
  else {

    # get sequence of ref, find all internal inframe-ATG
    refsequence <- unlist(BSgenome::getSeq(fasta, refCDS)) #
    startcodons <- Biostrings::matchPattern("ATG", refsequence) %>% IRanges::ranges()
    inframestarts <- startcodons[BiocGenerics::end(startcodons) %% 3 == 0 &
      BiocGenerics::start(startcodons) != 1]

    # return if no internal ATG is found
    if (length(inframestarts) > 0) {

      # This function attempts to map the XStringViews output back to refGRanges
      inframestartsingranges <- do.call("c", base::mapply(function(x, y) {
        start <- x - 1
        end <- length(refsequence) - y
        startcodoninGRanges <- resizeTranscript(refCDS, start, end)
        return(startcodoninGRanges)
      }, BiocGenerics::start(inframestarts), BiocGenerics::end(inframestarts)))

      # obtain 5'UTR length if query contain any of the inframe ATG
      if (any(inframestartsingranges %within% query)) {
        firststartgranges <- inframestartsingranges[inframestartsingranges %within% query][1]
        disjoint <- BiocGenerics::append(query, firststartgranges) %>%
          GenomicRanges::disjoin(with.revmap = T) %>%
          sort(decreasing = strand == "-") %>%
          as.data.frame() %>%
          dplyr::mutate(cumsum = cumsum(width))

        startcodonindex <- min(which(lengths(disjoint$revmap) == 2))
        if (startcodonindex > 1) {
          fiveUTRlength <- disjoint[startcodonindex - 1, ]$cumsum
        } else {
          fiveUTRlength <- 0
        }

        output$ORF_start <- "Internal ATG"
        output$fiveUTRlength <- fiveUTRlength

        return(output)
      } 
    }
  }
  
  # final chance of assigning CDS, just assign inframeCDS
  #get tailphases on refCDS
  S4Vectors::mcols(refCDS)$tailphase <- cumsum(width(refCDS) %% 3)%%3
  
  combinedgr <- BiocGenerics::append(query, refCDS)
  disjoint <- combinedgr %>% 
    GenomicRanges::disjoin(with.revmap = T) %>% 
    sort(decreasing = strand == "-")
  revmap <- S4Vectors::mcols(disjoint)$revmap
  disjoint$phase <- S4Vectors::mcols(combinedgr)$tailphase[unlist(revmap)] %>% 
    IRanges::relist(revmap)
  
  firstcdsindex <- min(which(lengths(disjoint$revmap) == 2))
  firstscds <- disjoint[firstcdsindex] %>% 
    as.data.frame() %>% 
    dplyr::mutate(phase = phase[[1]][2]) %>%
    dplyr::mutate(fivetrim = (width - phase)%%3)
  
  firstcdsrange <- resizeTranscript(disjoint[firstcdsindex], firstscds$fivetrim)
  disjoint <- BiocGenerics::append(query, firstcdsrange) %>%
    GenomicRanges::disjoin(with.revmap = T) %>%
    sort(decreasing = strand == "-") %>%
    as.data.frame() %>%
    dplyr::mutate(cumsum = cumsum(width))
  
  
  # retrieve index of segment upstream of start codon and return its cumsumwidth
  startcodonindex <- min(which(lengths(disjoint$revmap) == 2))
  if (startcodonindex > 1) {
    fiveUTRlength <- disjoint[startcodonindex - 1, ]$cumsum
  } else {
    fiveUTRlength <- 0
  }
  
  # update output list
  output$ORF_start <- "Inferred frame"
  output$fiveUTRlength <- fiveUTRlength
  
  return(output)
}

.getCDSstop <- function(query, fasta, fiveUTRlength) {

  # prepare output list
  output <- list(
    ORF_found = FALSE,
    threeUTRlength = 0
  )

  # append query GRanges to start from star codon, and retrieve seq
  queryCDS <- resizeTranscript(query, start = fiveUTRlength)
  queryseq <- unlist(BSgenome::getSeq(fasta, queryCDS))
  
  # prepare a dict of stop codons for pattern matching
  list_stopcodons <- Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
  pdict_stopcodons <- Biostrings::PDict(list_stopcodons)

  # search for in-frame stop codons
  stopcodons <- Biostrings::matchPDict(pdict_stopcodons, queryseq) %>%
    unlist() %>%
    as.data.frame() %>%
    dplyr::filter(end %% 3 == 0) %>%
    dplyr::arrange(start)

  # return if no in-frame stop codons are found
  if (nrow(stopcodons) == 0) {
    return(output)
  } else {

    # retrieve 3UTR length and update output file
    firststopcodon <- stopcodons[1, ]
    threeUTRlength <- length(queryseq) - firststopcodon$end
    output$threeUTRlength <- threeUTRlength
    output$ORF_found <- TRUE

    return(output)
  }
}

.getCDSranges <- function(query, fiveUTRlength, threeUTRlength, starttype) {

  # prepare output list
  output <- list("ORF_considered" = NA)

  # define global variables
  width <- seqnames <- strand <- phase <- type <- transcript_id <- NULL

  # resize query GRanges to ORF and renew metadata info
  CDSranges <- resizeTranscript(query, fiveUTRlength, (threeUTRlength + 3))
  if (length(CDSranges) == 0) {
    return(NULL)
  }
  
  CDSranges <- CDSranges %>%
    as.data.frame() %>% 
    dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start)) %>%
    dplyr::mutate(
      type = "CDS",
      transcript_id = S4Vectors::mcols(query)$transcript_id[1]
    ) %>%
    dplyr::mutate(phase = rev(cumsum(rev(width) %% 3) %% 3)) %>%
    dplyr::select(seqnames:end, strand, type, phase, transcript_id)
  CDSranges$built_from <- starttype
  output$ORF_considered <- CDSranges
  return(CDSranges)
}

.buildCDSchecks <- function(query, ref, argnames) {
  # check inputs are gtf
  if (any(!is_gtf(query, ref))){
    rlang::abort(sprintf(
      "%s is/are not gtf GRanges", 
      paste(argnames[1:2][!is_gtf(query, ref)], collapse = ',')))
  }
  
  # check if ref have CDS info
  if (!'CDS' %in% S4Vectors::mcols(ref)$type){
    rlang::abort(sprintf(
      "`%s` have missing CDS info",argnames[2])
    )
  }
  
  # catch unmatched seqlevels
  if (GenomeInfoDb::seqlevelsStyle(query)[1] != GenomeInfoDb::seqlevelsStyle(ref)[1]) {
    rlang::abort(sprintf(
      "`%s` and `%s` has unmatched seqlevel styles. try running: 
      \t\t%s <- matchSeqLevels(%s, %s)",
      argnames[1], argnames[2], argnames[1], argnames[1], argnames[2])
    )}
}

.prepq2r <- function(query, ref, query_exons, ref_exons, ref_cds, argnames) {
  
  # define global variables
  gene_id <- ref_transcript_id <- coverage <- transcript_id <- NULL
  # prepare q2r df
  query_ids = query %>% as.data.frame() %>%
    dplyr::filter(type == 'exon') %>% 
    dplyr::select(gene_id, transcript_id) %>%
    dplyr::distinct()
  ref_ids = ref %>% as.data.frame() %>%
    dplyr::filter(type == 'CDS') %>% 
    dplyr::select(gene_id, ref_transcript_id = transcript_id) %>%
    dplyr::distinct()
  q2r <- dplyr::left_join(query_ids, ref_ids, by = 'gene_id') %>%
    dplyr::select(-gene_id)
  
  # remove unmatched
  totaltx <- length(unique(q2r$transcript_id))
  unmatchedtx <- q2r %>%
    dplyr::filter(is.na(ref_transcript_id))
  if (nrow(q2r) == nrow(unmatchedtx)) {
    rlang::abort(sprintf(
      "all gene_ids in `%s` are not found in %s. try running:
      \t%s <- matchGeneIDs(%s, %s)",
      argnames[1], argnames[2], argnames[1], argnames[1], argnames[2]))
  } else if (nrow(unmatchedtx) > 0) {
    q2r <- dplyr::setdiff(q2r, unmatchedtx)
    rlang::warn(sprintf(
      "\t%s transcripts have unmatched gene_ids", 
      nrow(unmatchedtx)))
  }
  
  
  
  # simplify q2r if transcript_id and ref_transcript_id has same name
  q2rident <- q2r %>%
    dplyr::filter(transcript_id == ref_transcript_id)
  q2r <- q2r %>% 
    dplyr::filter(!transcript_id %in% q2rident$transcript_id) %>%
    dplyr::bind_rows(q2rident)
    
    
  
  # calculate coverage values
  out <- BiocParallel::bpmapply(function(x, y) {
    covrep <- .calcCoverage (query_exons[[x]], ref_exons[[y]], 'mean')
    return(covrep)
  }, q2r$transcript_id, q2r$ref_transcript_id,
  BPPARAM = BiocParallel::MulticoreParam()
  )
  q2r$coverage <- out
  
  # select best coevrage
  q2rcovs <- q2r %>%
    dplyr::arrange(transcript_id, dplyr::desc(coverage)) %>%
    dplyr::distinct(transcript_id, .keep_all = T)
  
  # check for transcript-pairs with no coverage
  nocovtx <- q2rcovs %>%
    dplyr::filter(coverage == 0)
  if (nrow(q2rcovs) == nrow(nocovtx)){
    rlang::abort(sprintf(
      "all transcripts in `%s` have no overlap with transcripts in %s",
      argnames[1], argnames[2]))
  } else if (nrow(nocovtx) > 0) {
    q2rcovs <- dplyr::setdiff(q2rcovs, nocovtx)
    rlang::warn(sprintf(
      "\t%s transcripts have no coverage with reference transcripts", 
      nrow(nocovtx)))
  }
  return(q2rcovs)
}

.calcCoverage <- function(tx1, tx2, over) {
  chrom <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(tx1)))
  cov <- suppressWarnings(GenomicRanges::coverage(c(tx1, tx2)))
  index <- which(names(cov) == chrom)
  cov <- cov[[index]]
  cov_val <- S4Vectors::runValue(cov)
  cov_len <- S4Vectors::runLength(cov)
  
  if (over == "mean") {
    denom <- sum(BiocGenerics::width(tx1), BiocGenerics::width(tx2)) / 2
  } else if (over == "query") {
    denom <- sum(BiocGenerics::width(tx1))
  } else if (over == "ref") {
    denom <- sum(BiocGenerics::width(tx2))
  }
  return(sum(cov_len[cov_val == 2]) / denom)
}

