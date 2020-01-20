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
extractCDSfeature <- function(cds, fasta, 
                              hmmer = 'default',
                              signalHsmm = NULL, 
                              ...) {

  # catch missing args
  mandargs <- c("cds", "fasta")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
  
  # get argnames and carry out checks
  argnames <- as.character(match.call())[-1]
  cds <- .extractCDSchecks(cds, fasta, argnames)
  feature <- .featurechecks(hmmer, signalHsmm)

  # define global variables
  . <- x <- y<- instop <- nincluded <- desc <- id <- NULL

  # get sequence
  aaSeq <- .getSequence(cds, fasta)
  
  
  # prepare output and run analysis
  output <- aaSeq %>% dplyr::select(id)

  # run hmmer function if requested
  if (!is.null(feature$hmmer)) {
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
      dplyr::rename_at(dplyr::vars(feature$hmmer), 
                       ~ paste0('hmmer.', feature$hmmer)) %>%
      dplyr::right_join(aaSeq %>% dplyr::select(id), by = 'id')
  }
  
  # run signalHsmm function
  if (!is.null(feature$signalHsmm)) {
    
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
      dplyr::select(id, feature$signalHsmm) %>% 
      dplyr::rename_at(dplyr::vars(feature$signalHsmm), 
                       ~ paste0('signalHsmm.', feature$signalHsmm)) %>%
      dplyr::left_join(output, ., by = 'id')
  }
  
  return(output)
}


.extractCDSchecks <- function(cds, fasta, argnames) {
  if (suppressWarnings(!has_consistentSeqlevels(cds, fasta))) {
    rlang::abort(sprintf(
      "`%s` and `%s` has unmatched seqlevel styles. 
Try running: %s <- matchSeqLevels(%s, %s)",
      argnames[1], argnames[2], argnames[1], argnames[1], argnames[2]))
  }
  # catch wrong cds class
  if (is_gtf(cds)) {
    cds <- S4Vectors::split(cds[cds$type == 'CDS'], ~transcript_id)
  }
  
  if (!is(cds, "GRangesList")) {
    rlang::abort("cds class type is not GRanges GTF or GRangesList")
  }
  return(sorteach(cds, exonorder))
}

.featurechecks <- function(hmmer, signalHsmm) {
  # check for correct feature column names for return
  hmmer.header <- c("name","acc","bias","desc","evalue","flags","hindex","ndom",
                    "nincluded","nregions","nreported","pvalue","score","taxid",
                    "pdb.id","bitscore","mlog.evalue")
  if ('default' %in% hmmer) {   # convert default into list of headers
    hmmer <- c('desc', 'evalue', 'nincluded',
               hmmer[hmmer != 'default'] 
    ) %>% unique()
  }
  if (any(!hmmer %in% hmmer.header)){
    fields_support <- hmmer %in% hmmer.header
    rlang::warn(sprintf('`%s` are not supported hmmer reported fields', 
                        hmmer[!fields_support]))
    if (sum(fields_support) == 0) {
      hmmer <- c('desc', 'evalue', 'nincluded')
      rlang::warn('returning default hmmer fields')
    } else {
      hmmer <- hmmer[fields_support]
    }
  }
  # check for correct field names for return
  signalHsmm.header <- c("sp_probability","sp_start","sp_end",
                         "struc","prot","name","str_approx")
  if ('default' %in% signalHsmm) {   # convert default into list of headers
    signalHsmm <- c("sp_probability",
                    signalHsmm[signalHsmm != 'default'] 
    ) %>% unique()
  }
  if (any(!signalHsmm %in% signalHsmm.header)){
    fields_support <- signalHsmm %in% signalHsmm.header
    rlang::warn(sprintf('`%s` are not supported signalHsmm reported fields', 
                        signalHsmm[!fields_support]))
    if (sum(fields_support) == 0) {
      signalHsmm <- "sp_probability"
      rlang::warn('returning default signalHsmm fields')
    } else {
      signalHsmm <- signalHsmm[fields_support]
    }
  }
  return(list(hmmer = hmmer, signalHsmm = signalHsmm))
}

.getSequence <- function(cds, fasta) {
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
    rlang::warn(sprintf("%s cds entries do not start with ATG", sum(aaSeq$noATG)))
  }
  if (T %in% aaSeq$instop) {
    aaSeq <- suppressWarnings(aaSeq %>% 
                                dplyr::rowwise() %>%
                                dplyr::mutate(x = ifelse(instop == T, 
                                                         paste(y[1:which(y == '*')-1], collapse = ''),
                                                         x)) %>%
                                dplyr::mutate(y = strsplit(x, split = "")) %>%
                                dplyr::ungroup())
    
    rlang::warn(sprintf("%s cds entries contain internal stop_codon. These proteins have been truncated", sum(aaSeq$instop)))
    if ('' %in% aaSeq$x) {
      rlang::warn(sprintf("After truncation, %s cds have no coding sequences. These entries were not analyzed", sum(aaSeq$x == '')))
      aaSeq <- aaSeq[aaSeq$x != '',]
    }
  }
  return(aaSeq)
}

