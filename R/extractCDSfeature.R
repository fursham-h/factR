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
  aaSeq <- suppressWarnings(Biostrings::translate(cdsSeq)) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("id") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(y = strsplit(substr(x, 1, nchar(x) - 1), split = "")) %>%
    dplyr::mutate(startendaa = paste0(substr(x, 1, 1), substr(x, nchar(x), nchar(x)))) %>%
    dplyr::ungroup()

  if (!"M*" %in% aaSeq$startendaa) {
    aaSeqORF <- aaSeq %>%
      dplyr::filter(startendaa == "M*")
    if (nrow(aaSeqORF) == 0) {
      stop("All cds are not ORFs. Please check cds GRangesList or fasta input")
    } else if (nrow(aaSeqORF) < nrow(aaSeq)) {
      num_nonCDS <- nrow(aaSeq) - nrow(aaSeqORF)
      warning(sprintf("%s cds entries are not ORFs and have not been analysed"))
    }
    aaSeq <- aaSeqORF
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
        return("NA")
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
