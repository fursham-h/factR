#' Construct query CDS using reference as guide
#'
#' @description buildCDS will attempt to construct query CDS information by
#' firstly deriving its start codon. This is done searching transcripts for an
#' annotated start codon, provided by its corresponding reference CDS.
#' If this fails, function will find the next common internal ATG codon and sets
#' that as the start codon. buildCDS will then search for an in-frame stop codon
#' and will return the CDS GRanges if found.
#' @param query
#' GRangesList object containing exons for each query transcript. Transcripts
#' have to be listed in query2ref dataframe, else CDS will not be constructed
#' for query missing from query2ref
#' @param refCDS
#' GRangesList object containing CDS for each reference transcript.
#' @param fasta
#' BSgenome or Biostrings object containing genomic sequence
#' @param query2ref
#' Dataframe with at least 2 columns: query transcript_id and its paired
#' reference transcript_id. Query and ref transcript_ids have to match transcript
#' names in query and refCDS objects. Transcripts with missing corrresponding
#' GRanges object will return error
#' @param ids
#' Numeric vector stating which columns of query2ref dataframe contain the
#' query and reference transcript_ids respectively.
#' @param coverage
#' Integer stating which column of query2ref dataframe contain percent coverage
#' between query and reference transcripts. Providing a column with coverage
#' values will speed up CDS building process. Query transcripts that share 100%
#' coverage with reference CDS will be assigned the reference CDS and skip the
#' CDS searching process. See getCoverages function to calculate coverage values
#' (default:NULL)
#'
#' @return
#' GRangesList object containing CDS for each query transcript
#' @export
#' @author Fursham Hamid
#'
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' buildCDS(query_exons, ref_cds, Mmusculus, q2rcovs, coverage = 3)
buildCDS <- function(query, refCDS, fasta, query2ref,
                     ids = c(1, 2), coverage = NULL) {

  # catch missing args
  mandargs <- c("query", "refCDS", "fasta", "query2ref")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }

  # define global variables
  group_name <- strand <- width <- phase <- transcript_id <- NULL

  # check if query and cds are GRlist
  if (!is(query, "GRangesList") | !is(refCDS, "GRangesList")) {
    txtype <- is(query)[1]
    cdstype <- is(refCDS)[1]

    error <- c(!is(query, "GRangesList"), !is(refCDS, "GRangesList"))
    obj <- paste(c("query", "refCDS")[error], collapse = ",")
    types <- paste(c(is(query)[1], is(refCDS)[1])[error], collapse = ",")

    stop(sprintf(
      "Incompatile input object. %s is type %s respectively",
      obj, types
    ))
  }

  # catch unmatched seqlevels
  if (GenomeInfoDb::seqlevelsStyle(query) != GenomeInfoDb::seqlevelsStyle(refCDS)) {
    stop("query and refCDS has unmatched seqlevel styles. try matching using matchSeqLevels function")
  }

  # extract colnames and try catching wrong indices
  outCDS <- NULL
  tryCatch(
    {
      txname <- names(query2ref[ids[1]])
      refname <- names(query2ref[ids[2]])
    },
    error = function(e) {
      stop("column indices in `ids` are not found in query2ref")
    }
  )
  if (txname == refname) {
    stop("`ids` contain duplicate indices")
  }

  # sanity check if all tx in q2r have GRanges object
  if (!all(query2ref[[ids[1]]] %in% names(query))) {
    missing <- sum(!query2ref[[ids[1]]] %in% names(query))
    stop(sprintf(
      "%s query transcripts have missing GRanges object",
      missing
    ))
  }
  if (!all(query2ref[[ids[2]]] %in% names(refCDS))) {
    missing <- sum(!query2ref[[ids[2]]] %in% names(refCDS))
    stop(sprintf(
      "%s reference CDSs have missing GRanges object",
      missing
    ))
  }

  # # sanity check if query and ref names are in q2r df
  if (all(!names(query) %in% query2ref[[ids[1]]])) {
    unannotatedq <- sum((!names(query) %in% query2ref[ids[1]]))
    rlang::warn(sprintf(
      "%s query transcript ids were missing from query2ref df",
      unannotatedq
    ))
  }
  if (all(!names(refCDS) %in% query2ref[[ids[2]]])) {
    unannotatedr <- sum((!names(refCDS) %in% query2ref[ids[2]]))
    rlang::warn(sprintf(
      "%s reference CDS ids were missing from query2ref df",
      unannotatedr
    ))
  }

  # create CDS list for tx with coverage of 1
  if (!is.null(coverage)) {
    covname <- names(query2ref)[coverage] # extract cov colname
    # subset q2r for full coverages
    fullcovs <- query2ref %>%
      dplyr::filter(!!as.symbol(covname) == 1)
    query2ref <- query2ref %>%
      dplyr::filter(!!as.symbol(covname) != 1)
    # prepare outputCDS for full coverages
    outCDS <- refCDS %>%
      as.data.frame() %>%
      dplyr::filter(group_name %in% fullcovs[[ids[2]]]) %>%
      dplyr::select(group_name:strand) %>%
      dplyr::mutate(type = "CDS") %>%
      dplyr::left_join(fullcovs[ids],
        by = c("group_name" = refname)
      ) %>%
      dplyr::group_by(group_name) %>%
      dplyr::arrange(ifelse(strand == '-', dplyr::desc(start), start)) %>%
      dplyr::mutate(phase = rev(cumsum(rev(width) %% 3) %% 3)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-group_name) %>%
      dplyr::mutate(built_from = 'Full coverage')
    
    # remove no coverage comparisons
    if (0 %in% query2ref[[coverage]]) {
      nocov <- sum(query2ref[[coverage]] == 0)
      rlang::warn(sprintf(
        "%s comparisons in query2ref have 0 coverage. These were not analyzed",
        nocov
      ))
      
      query2ref <- query2ref %>%
        dplyr::filter(!!as.symbol(covname) > 0)
    }
  }

  # create CDS list for all remaining tx
  out <- BiocParallel::bpmapply(function(x, y) {
    CDSreport <- getCDS_(query[x], refCDS[y], fasta) %>%
      as.data.frame()
    return(CDSreport)
  }, query2ref[[txname]], query2ref[[refname]],
  BPPARAM = BiocParallel::MulticoreParam(), SIMPLIFY = F
  ) %>%
    dplyr::bind_rows()
  outCDS <- suppressWarnings(dplyr::bind_rows(outCDS, out) %>%
    dplyr::mutate(group_name = transcript_id) %>%
    dplyr::arrange(group_name, ifelse(strand == '-', dplyr::desc(start), start)) %>%
    GenomicRanges::makeGRangesListFromDataFrame(
      split.field = "group_name",
      keep.extra.columns = T
    ))

  # warn users if program fails to find CDS for some transcripts
  if (length(outCDS) < length(query)) {
    missingCDS <- length(query) - length(outCDS)
    rlang::warn(sprintf("Unable to find CDS for %s transcripts", missingCDS))
  }

  return(outCDS)
}

getCDS_ <- function(query, CDS, fasta) {

  # prepare output list
  output <- list(
    ORF_considered = as.character(NA),
    ORF_start = as.character("Not found"),
    fiveUTRlength = 0,
    threeUTRlength = 0,
    ORF_found = FALSE
  )
  
  # Extract info and sort all GRanges first
  strand <- as.character(BiocGenerics::strand(query[[1]]))[1]
  queryTx <- BiocGenerics::sort(query[[1]], decreasing = strand == "-")
  knownCDS <- BiocGenerics::sort(CDS[[1]], decreasing = strand == "-")
  S4Vectors::mcols(queryTx)$transcript_id <- names(query)
  # attempt to find an aligned start codon
  report <- getCDSstart_(queryTx, knownCDS, fasta)
  output <- utils::modifyList(output, report) # update output

  # return if no start codon is found
  if (output$ORF_start == "Not found") {
    return()
  }

  # attempt to search for an in-frame stop codon
  report <- getCDSstop_(queryTx, fasta, output$fiveUTRlength)
  output <- utils::modifyList(output, report) # update output

  # return if no stop codon is found
  if (output$ORF_found == FALSE) {
    return()
  }

  # build new ORF Granges
  report <- getCDSranges_(queryTx, output$fiveUTRlength, output$threeUTRlength,
                          output$ORF_start)
  output <- utils::modifyList(output, report) # update output
  # return(output[c('ORF_considered', 'ORF_start', 'ORF_found')])
  return(output$ORF_considered)
}

getCDSstart_ <- function(query, refCDS, fasta) {

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
    output$ORF_start <- "Shared ATG"
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
    utils::relist(revmap)
  
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

getCDSstop_ <- function(query, fasta, fiveUTRlength) {

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

getCDSranges_ <- function(query, fiveUTRlength, threeUTRlength, starttype) {

  # prepare output list
  output <- list("ORF_considered" = NA)

  # define global variables
  width <- seqnames <- strand <- phase <- type <- transcript_id <- NULL

  # resize query GRanges to ORF and renew metadata info
  CDSranges <- resizeTranscript(query, fiveUTRlength, (threeUTRlength + 3))
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
  return(output)
}
