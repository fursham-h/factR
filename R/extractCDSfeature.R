#' Extract protein features from coding mRNAs
#'
#' @param cds
#' GRanges object or GRangesList object containing coding regions (CDS)
#' for each transcript.
#' @param fasta
#' BSgenome or Biostrings object containing genomic sequence
#' @param which
#' List containing names of transcripts from cds to filter for analysis
#' @param hmmer
#' List of fields to report from hmmer protein domain search. 'default' will
#' return 'desc', 'evalue' and 'nincluded' columns. See ?bio3d::hmmer for a list
#' of returned fields. Set argument to NULL to skip hmmer analysis
#' @param ...
#' Other features to test and its fields to return. Written as named list with
#' features as names and list of fields as arguements (eg. signalHsmm = 'default')
#' Current supported features:
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
extractCDSfeature <- function(cds, fasta, which = NULL,
                              hmmer = 'default', ...) {

  # catch missing args
  mandargs <- c("cds", "fasta")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  
  # extract features to test as list and check for supported features
  features <- c(as.list(environment())['hmmer'], list(...))
  if (any(!names(features) %in% c('hmmer', 'signalHsmm'))) {
    features_support <- names(features) %in% c('hmmer', 'signalHsmm')
    warning(sprintf('feature(s) `%s` is not supported and its analysis is not run', 
                    paste(features[!features_support], collapse = ',')))
    features <- features[features_support]
  }
  
  # check for correct feature column names for return
  hmmer.header <- c("name","acc","bias","desc","evalue","flags","hindex","ndom",
                    "nincluded","nregions","nreported","pvalue","score","taxid",
                    "pdb.id","bitscore","mlog.evalue")
  if ('default' %in% features$hmmer) {   # convert default into list of headers
    features$hmmer <- c('desc', 'evalue', 'nincluded',
                         features$hmmer[features$hmmer != 'default'] 
                         ) %>% unique()
  }
  if (any(!features$hmmer %in% hmmer.header)){
    fields_support <- features$hmmer %in% hmmer.header
    warning(sprintf('`%s` are not supported hmmer reported fields', 
                    features$hmmer[!fields_support]))
    if (sum(fields_support) == 0) {
      features$hmmer <- c('desc', 'evalue', 'nincluded')
      warning('returning default hmmer fields')
    } else {
      features$hmmer <- features$hmmer[fields_support]
    }
  }
  
  
  # define global variables
  . <- x <- startendaa <- nincluded <- desc <- id <- NULL

  # catch unmatched seqlevels
  if (GenomeInfoDb::seqlevelsStyle(cds) != GenomeInfoDb::seqlevelsStyle(fasta)) {
    stop("cds and fasta has unmatched seqlevel styles. try matching using matchSeqLevels function")
  }

  # catch wrong cds class
  if (!is(cds, "GRangesList")) {
    stop("cds class type is not GRangesList")
  }

  # subset cds by which
  if (!is.null(which)) {
    which_matched <- which[which %in% names(cds)]
    if (length(which_matched) == 0) {
      stop("transcript names in `which` is not found in cds")
    } else {
      cds <- cds[names(cds) %in% which]
    }

    if (length(which_matched) != length(which)) {
      num_unmatched <- length(which) - length(which_matched)
      warning(sprintf(
        "%s transcript names in `which` is missing from cds names",
        num_unmatched
      ))
    }
  }

  # get sequence
  cdsSeq <- GenomicFeatures::extractTranscriptSeqs(fasta, cds)
  aaSeq <- suppressWarnings(Biostrings::translate(cdsSeq, if.fuzzy.codon = 'solve')) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("id") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(y = strsplit(x, split = "")) %>%
    dplyr::mutate(noATG = ifelse(y[[1]] != 'M',T,F)) %>%
    dplyr::mutate(instop = ifelse('*' %in% y,T,F)) %>% 
    dplyr::ungroup()
  
  # check for ATG and internal stop_codon, truncate proteins with internal stop codon
  ## and remove entries without proteins after truncation
  if (T %in% aaSeq$noATG) {
    warning(sprintf("%s cds entries do not start with ATG", sum(aaSeq$noATG)))
  }
  if (T %in% aaSeq$instop) {
    aaSeq <- suppressWarnings(aaSeq %>% 
      dplyr::rowwise() %>%
      dplyr::mutate(x = ifelse(instop == T, 
                               paste(y[1:which(y == '*')-1], collapse = ''),
                               x)) %>%
      dplyr::mutate(y = strsplit(x, split = "")))
      
    warning(sprintf("%s cds entries contain internal stop_codon. These proteins have been truncated", sum(aaSeq$instop)))
    if ('' %in% aaSeq$x) {
      warning(sprintf("After truncation, %s cds have no coding sequences. These entries were not analyzed", sum(aaSeq$x == '')))
      aaSeq <- aaSeq[aaSeq$x != '',]
    }
  }
  
  # prepare output
  output <- aaSeq %>% dplyr::select(id)

  # run hmmer function if requested
  if (!is.null(features$hmmer)) {
    output <- BiocParallel::bplapply(seq_len(nrow(aaSeq)), function(y) {
      # account for return errors
      report <- tryCatch(
        bio3d::hmmer(aaSeq[y,]$x,
                     type = "hmmscan",
                     verbose = F, db = "superfamily"
        ),
        error = function(e) NULL
      )
      
      if (is.null(report)) {
        return(NULL)
      } else {
        report <- report$hit.tbl %>%
          dplyr::mutate(id = aaSeq[y,]$id)

        return(report)
      }
    }, BPPARAM = BiocParallel::MulticoreParam()) %>%
      dplyr::bind_rows()
    
    # subset columns based on requested fields
    output <- output %>% 
      dplyr::select(id, features$hmmer) %>% 
      dplyr::rename_at(dplyr::vars(features$hmmer), 
                       ~ paste0('hmmer.', features$hmmer)) %>%
      dplyr::right_join(aaSeq %>% dplyr::select(id), by = 'id')
  }
  

  
  # run signalHsmm function
  if (!is.null(features$signalHsmm)) {
    # check for correct field names for return
    signalHsmm.header <- c("sp_probability","sp_start","sp_end",
                        "struc","prot","name","str_approx")
    if ('default' %in% features$signalHsmm) {   # convert default into list of headers
      features$signalHsmm <- c("sp_probability",
                            features$signalHsmm[features$signalHsmm != 'default'] 
      ) %>% unique()
    }
    if (any(!features$signalHsmm %in% signalHsmm.header)){
      fields_support <- features$signalHsmm %in% signalHsmm.header
      warning(sprintf('`%s` are not supported signalHsmm reported fields', 
                      features$signalHsmm[!fields_support]))
      if (sum(fields_support) == 0) {
        features$signalHsmm <- "sp_probability"
        warning('returning default signalHsmm fields')
      } else {
        features$signalHsmm <- features$signalHsmm[fields_support]
      }
    }
    
    # run signalHsmm
    outSP <- BiocParallel::bplapply(seq_len(nrow(aaSeq)), function(x) {
      report <- signalHsmm::run_signalHsmm(aaSeq[x,]$y)[[1]]
      report$struc <- paste(report$struc, collapse = '')
      report$prot <- paste(report$prot, collapse = '')
      
      report <- report %>% unlist() %>% dplyr::bind_rows() %>%
        dplyr::mutate(id = aaSeq[x,]$id)
      return(report)
    }, BPPARAM = BiocParallel::MulticoreParam()) %>% dplyr::bind_rows()
    
    # combine outSP with output
    output <- outSP %>%
      dplyr::select(id, features$signalHsmm) %>% 
      dplyr::rename_at(dplyr::vars(features$signalHsmm), 
                       ~ paste0('signalHsmm.', features$signalHsmm)) %>%
      dplyr::left_join(output, ., by = 'id')
  }
  
  return(output)
}
