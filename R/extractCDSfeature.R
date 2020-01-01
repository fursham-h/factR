#' Extract protein features from coding mRNAs
#'
#' @param cds 
#' GRanges object or GRangesList object containing coding regions (CDS)
#' for each transcript. 
#' @param fasta 
#' BSgenome or Biostrings object containing genomic sequence
#' @param which 
#' List containing names of transcripts from cds to filter for analysis
#'
#' @return
#' Dataframe containing protein features for each cds entry
#' Current features tested:
#' (1) Presence and type of protein domains
#' (2) Probability of signal peptide sequence
#' @export
#'
#' @examples
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' extractCDSfeature(query_cds, Mmusculus)
extractCDSfeature <- function(cds, fasta, which = NULL) {
  
  # catch missing args
  mandargs <- c("cds", "fasta")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  # define global variables
  x<- startendaa <- nincluded <- desc <- id <- NULL
  
  # catch unmatched seqlevels
  if (GenomeInfoDb::seqlevelsStyle(cds) != GenomeInfoDb::seqlevelsStyle(fasta)) {
    stop("cds and fasta has unmatched seqlevel styles. try matching using matchSeqLevels function")
  }
  
  # catch wrong cds class
  if (!is(cds, "GRangesList")) {
    stop('cds class type is not GRangesList')
  }
  
  # subset cds by which
  if (!is.null(which)){
    which_matched <- which[which %in% names(cds)]
    if (length(which_matched) == 0){
      stop('transcript names in `which` is not found in cds')
    } else {
      cds <- cds[names(cds) %in% which]
    }
    
    if (length(which_matched) != length(which)){
      num_unmatched <- length(which) - length(which_matched)
      warning(sprintf('%s transcript names in `which` is missing from cds names',
                      num_unmatched))
    }
  }
  
  # get sequence
  cdsSeq <- GenomicFeatures::extractTranscriptSeqs(fasta, cds)
  aaSeq <- suppressWarnings(Biostrings::translate(cdsSeq)) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column('id') %>%
    dplyr::rowwise() %>%
    dplyr::mutate(y = strsplit(substr(x,1,nchar(x)-1), split = '')) %>%
    dplyr::mutate(startendaa = paste0(substr(x,1,1), substr(x,nchar(x),nchar(x)))) %>%
    dplyr::ungroup()
  
  if (!'M*' %in% aaSeq$startendaa) {
    aaSeqORF <- aaSeq %>%
      dplyr::filter(startendaa == 'M*')
    if (nrow(aaSeqORF) == 0) {
      stop('All cds are not ORFs. Please check cds GRangesList or fasta input')
    } else if (nrow(aaSeqORF) < nrow(aaSeq)){
      num_nonCDS <- nrow(aaSeq) - nrow(aaSeqORF)
      warning(sprintf('%s cds entries are not ORFs and have not been analysed'))
    }
    aaSeq <- aaSeqORF
  }

  # run hmmer function
  outDomain <- BiocParallel::bplapply(aaSeq$x, function(x) {
    
    # account for return errors
    report <- tryCatch(
      bio3d::hmmer(x, type = 'hmmscan', 
                   verbose = F, db = 'superfamily'),
      error = function(e) NULL
    )
    
    if (is.null(report)){
      return('NA')
    } else {
      report <- report$hit.tbl %>% 
        dplyr::rowwise() %>%
        dplyr::mutate(id = paste0(nincluded,'x',desc)) %>%
        dplyr::ungroup()
      return(report$id %>% paste(collapse = ';'))
    }
    
  }, BPPARAM = BiocParallel::MulticoreParam()) 
  
  # run signalHsmm function
  outSP <- BiocParallel::bplapply(aaSeq$y, function(x) {
    report <- signalHsmm::run_signalHsmm(x)
    return(report[[1]]$sp_probability)
  }, BPPARAM = BiocParallel::MulticoreParam())
  
  # prepare output df
  aaSeq <- aaSeq %>%
    dplyr::select(id) %>% 
    dplyr::bind_cols(list('Domains' = unlist(outDomain))) %>%
    dplyr::bind_cols(list('SignalP_prob' = unlist(outSP)))
  return(aaSeq)
}