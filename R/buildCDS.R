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
  
  # fill query with phase info
  mergemetadata <- IRanges::mergeByOverlaps(query[query$type == 'exon'], 
                                           ref[ref$type == 'CDS'], type = 'equal')
  filled_query <- mergemetadata$`query[query$type == "exon"]`
  filled_query$phase <- mergemetadata$`ref[ref$type == "CDS"]`$phase
  
  new_query <- query %>%
    as.data.frame() %>%
    dplyr::bind_rows(as.data.frame(filled_query)) %>%
    dplyr::distinct(transcript_id, start, end)

  # makeGRangesList
  query_exons <- S4Vectors::split(query[query$type == 'exon'], ~transcript_id)
  ref_cds <- S4Vectors::split(ref[ref$type == 'CDS'], ~transcript_id)
  ref_exons <- S4Vectors::split(ref[ref$type == 'exon'], ~transcript_id)
  nc_ref_exons <- ref_exons[!names(ref_exons) %in% names(ref_cds)] 
  ref_exons <- ref_exons[names(ref_exons) %in% names(ref_cds)] 
  
  # prepare q2r
  q2r <- .prepq2r(query_exons, ref_exons, ref_cds, nc_ref_exons)


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
      dplyr::filter(group %in% fullcovs[[3]]) %>%
      dplyr::select(group:strand) %>%
      dplyr::mutate(type = "CDS") %>%
      dplyr::left_join(fullcovs[,c(1,3)],
                       by = c("group" = "ref_transcript_id")
      ) %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start)) %>%
      dplyr::mutate(phase = rev(cumsum(rev(width) %% 3) %% 3)) %>%
      dplyr::ungroup() %>%
      dplyr::select(seqnames:transcript_id) %>%
      dplyr::mutate(built_from = 'Full coverage')
  }
  

  # create CDS list for all remaining tx
  out <- BiocParallel::bpmapply(function(x, y, z) {
    CDSreport <- .getthisCDS(y, z, fasta)
    if (!is.null(CDSreport)) {
      CDSreport$transcript_id <- x
    }
    return(CDSreport)
  }, query2ref$transcript_id, 
  query_exons[query2ref$transcript_index], 
  ref_cds[query2ref$ref_transcript_id],
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
  #S4Vectors::mcols(queryTx)$transcript_id <- names(query)
  # attempt to find an aligned start codon
  report <- .getCDSstart(queryTx, knownCDS, fasta)
  output <- utils::modifyList(output, report) # update output

  # return if no start codon is found
  if (output$ORF_start == "Not found") {
    return(NULL)
  }

  # attempt to search for an in-frame stop codon
  report <- .getCDSstop(queryTx, fasta, output$fiveUTRlength)
  output <- utils::modifyList(output, report) # update output

  # return if no stop codon is found
  if (output$ORF_found == FALSE) {
    return(NULL)
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
  #return(output)
  
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
    dplyr::select(seqnames:end, strand, type, phase)
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

.prepq2r <- function(query_exons, ref_exons, ref_cds, nc_ref_exons) {
  # prepare q2r
  fulloverlap <- GenomicRanges::findOverlaps(query_exons, ref_exons, 
                                             type = 'equal', select = "first")
  startoverlappos <- GenomicRanges::findOverlaps(range(query_exons), range(ref_exons), 
                                                 type = 'start', select = "first")
  startoverlapneg <- GenomicRanges::findOverlaps(range(query_exons), range(ref_exons), 
                                                 type = 'end', select = "first")
  anyoverlap <- GenomicRanges::findOverlaps(query_exons, ref_cds, 
                                            type = 'any', select = "first")
  ncoverlap <- GenomicRanges::findOverlaps(query_exons, nc_ref_exons, 
                                           type = 'equal', select = "first")
  
  strandList <- as.character(runValue(strand(query_exons)))
  
  fullq2r <- data.frame('transcript_id' = names(query_exons),
                        'transcript_index' = 1:length(query_exons),
                        'ref_transcript_id' = fulloverlap, 
                        'strand' = strandList,
                        stringsAsFactors = F) %>%
    dplyr::filter(!is.na(ref_transcript_id)) %>%
    dplyr::mutate(coverage = 1)
  ncq2r <- data.frame('transcript_id' = names(query_exons),
                      'transcript_index' = 1:length(query_exons),
                      'ref_transcript_id' = ncoverlap, 
                      'strand' = strandList,
                      stringsAsFactors = F) %>%
    dplyr::filter(!is.na(ref_transcript_id)) %>%
    dplyr::mutate(coverage = 0)
  
  startq2rpos <- data.frame('transcript_id' = names(query_exons),
                            'transcript_index' = 1:length(query_exons),
                            'ref_transcript_id' = startoverlappos, 
                            'strand' = strandList,
                            stringsAsFactors = F) %>%
    dplyr::filter(!is.na(ref_transcript_id) & strand != '-') %>%
    dplyr::mutate(coverage = 2)
  startq2rneg <- data.frame('transcript_id' = names(query_exons),
                            'transcript_index' = 1:length(query_exons),
                            'ref_transcript_id' = startoverlapneg, 
                            'strand' = strandList,
                            stringsAsFactors = F) %>%
    dplyr::filter(!is.na(ref_transcript_id) & strand == '-') %>%
    dplyr::mutate(coverage = 2)
  
  
  anyq2r <- data.frame('transcript_id' = names(query_exons),
                       'transcript_index' = 1:length(query_exons),
                       'ref_transcript_id' = anyoverlap, 
                       'strand' = strandList,
                       stringsAsFactors = F) %>%
    dplyr::filter(!is.na(ref_transcript_id)) %>%
    dplyr::mutate(coverage = 3)
  
  q2r <- dplyr::bind_rows(fullq2r, ncq2r, startq2rpos, startq2rneg, anyq2r) %>%
    dplyr::distinct(transcript_id, .keep_all  =T) %>%
    dplyr::filter(coverage > 0)
  # report unmatched tx
  
  return(q2r)
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


# 
# # prepare a dict of stop codons for pattern matching
# list_stopcodons <- Biostrings::DNAStringSet(c("TAA", "TAG", "TGA"))
# pdict_stopcodons <- Biostrings::PDict(list_stopcodons)

# # testingggg
# out <- BiocParallel::bpmapply(function(x, y, z) {
#   
#   strand = as.character(strand(y)[1])
#   y <- sort(y, decreasing = strand == '-')
#   z <- sort(z, decreasing = strand == '-')
#   z <- resize(z, width = width(z) - z$phase, fix = 'end')
#   a <- subsetByOverlaps(z, y, type = 'within')
#   
#   if (length(a) == 0) {
#     return(NULL)
#   }
#   
#   starttoend <- dplyr::bind_cols(as.data.frame(range(y)), as.data.frame(range(a[1]))) %>%
#     dplyr::mutate(newstart = ifelse(strand == '-', start ,start1)) %>%
#     mutate(newend = ifelse(strand == '-', end1 ,end)) %>% 
#     dplyr::select(seqnames, start = newstart, end = newend, strand) %>% 
#     GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
#   
#   starttoend <- GenomicRanges::intersect(y, starttoend)
#   queryseq <- unlist(BSgenome::getSeq(fasta, starttoend))
#   
#   
#   # search for in-frame stop codons
#   stopcodons <- Biostrings::matchPDict(pdict_stopcodons, queryseq) %>%
#     unlist() %>%
#     as.data.frame() %>%
#     dplyr::filter(end %% 3 == 0) %>%
#     dplyr::arrange(start)
#   
#   # return if no in-frame stop codons are found
#   if (nrow(stopcodons) == 0) {
#     return(NULL)
#   } 
#   threeUTRlength <- length(queryseq) - stopcodons[1,2] + 3
#   cds <- resizeTranscript(starttoend, end = threeUTRlength)
#   
#   # return if no in-frame stop codons are found
#   if (length(cds) == 0) {
#     return(NULL)
#   } 
#   cds$type <- 'CDS'
#   cds$transcript_id <- x
#   cds$phase <- rev(cumsum(rev(width(cds)) %% 3) %% 3) 
#   return(as.data.frame(cds))
#   
# }, query2ref$transcript_id, 
# query_exons[query2ref$transcript_index], 
# ref_cds[query2ref$ref_transcript_id],
# BPPARAM = BiocParallel::MulticoreParam(), SIMPLIFY = F
# ) %>% dplyr::bind_rows()
