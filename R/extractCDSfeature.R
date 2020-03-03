#' Extract protein features from coding mRNAs
#'
#' @param x
#' Can be a GRanges object containing 'CDS' features in GTF format
#'
#' Can be a GRangesList object containing CDS ranges for each transcript
#' @param fasta
#' BSgenome or Biostrings object containing genomic sequence
#' @param ...
#' Logical conditions to pass to dplyr::filter to subset transcripts for analysis.
#' Variables are metadata information found in `x` and multiple conditions can be
#' provided delimited by comma. Example: transcript_id == "transcript1"
#' @param hmmer
#' List of fields to report from hmmer protein domain search. 'default' will
#' return 'desc', 'evalue' and 'nincluded' columns. See ?bio3d::hmmer for a list
#' of returned fields. Set argument to NULL to skip hmmer analysis
#' @param signalHsmm
#' List of fields to report from signalHsmmm prediction tool. By default, argument
#' is set to NULL which skips the analysis. Seting argument to 'default' will
#' return 'sp_probability' columns. See ?signalHsmm::signalHsmm for a list
#' of returned fields.
#'
#' @return
#' Dataframe containing protein features for each cds entry
#' @export
getCDSfeature <- function(x, fasta, ..., plot = T) {

  # catch missing args
  mandargs <- c("x", "fasta")
  passed <- names(as.list(match.call())[-1])
  if (any(!mandargs %in% passed)) {
    stop(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }

  # get argnames and carry out checks
  argnames <- as.character(match.call())[-1]
  cds <- .extractCDSchecks(x, fasta, argnames, ...)
  #feature <- .featurechecks(hmmer, signalHsmm)

  # define global variables
  . <- id <- NULL

  # get sequence
  aaSeq <- .getSequence(cds, fasta)

  # prepare output and run analysis
  output <- aaSeq %>% dplyr::select(id)

  # run hmmer function if requested
  url <- paste("https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan")
  url.opts <- list(httpheader = "Expect:", httpheader = "Accept:text/xml", 
                   verbose = F, followlocation = TRUE)
  #if (!is.null(feature$hmmer)) {
    output <- BiocParallel::bplapply(seq_len(nrow(aaSeq)), function(y) {
      # account for return errors
      report <- tryCatch(
        .getdomains(url, curl.opts, aaSeq[y, ]$x, aaSeq[y,]$id, nchar(aaSeq[y, ]$x), y),
        error = function(e) NULL
      )

      if (is.null(report)) {
        return(NULL)
      } else {
        return(report)
      }
    }, BPPARAM = BiocParallel::MulticoreParam()) %>%
      dplyr::bind_rows()
    
    plot.out <- draw_canvas(output) %>% 
      draw_chains(output) %>% 
      draw_domains(output, label_domains = F) +
      theme_bw() + # white background
      theme(panel.grid.minor=element_blank(), 
            panel.grid.major=element_blank()) +
      theme(axis.ticks = element_blank(), 
            axis.text.y = element_blank()) +
      theme(panel.border = element_blank())
    
    if (plot) {
      plot.out
    }

    # subset columns based on requested fields
    table.out <- output %>%
      dplyr::select(id = entryName, description, eval, begin, end) %>%
      dplyr::filter(type == "DOMAIN") %>%
      dplyr::right_join(aaSeq %>% dplyr::select(id), by = "id")
  #}

  # run signalHsmm function
  # if (!is.null(feature$signalHsmm)) {
  # 
  #   # run signalHsmm
  #   outSP <- BiocParallel::bplapply(seq_len(nrow(aaSeq)), function(x) {
  #     report <- signalHsmm::run_signalHsmm(aaSeq[x, ]$y)[[1]]
  #     report$struc <- paste(report$struc, collapse = "")
  #     report$prot <- paste(report$prot, collapse = "")
  # 
  #     report <- report %>%
  #       unlist() %>%
  #       dplyr::bind_rows() %>%
  #       dplyr::mutate(id = aaSeq[x, ]$id)
  #     return(report)
  #   }, BPPARAM = BiocParallel::MulticoreParam()) %>% dplyr::bind_rows()
  # 
  #   # combine outSP with output
  #   output <- outSP %>%
  #     dplyr::select(id, feature$signalHsmm) %>%
  #     dplyr::rename_at(
  #       dplyr::vars(feature$signalHsmm),
  #       ~ paste0("signalHsmm.", feature$signalHsmm)
  #     ) %>%
  #     dplyr::left_join(output, ., by = "id")
  # }

  return(table.out)
}


.extractCDSchecks <- function(cds, fasta, argnames, ...) {
  # define global variables
  exonorder <- NULL

  if (suppressWarnings(!has_consistentSeqlevels(cds, fasta))) {
    rlang::abort(sprintf(
      "`%s` and `%s` has unmatched seqlevel styles. 
Try running: %s <- matchSeqLevels(%s, %s)",
      argnames[1], argnames[2], argnames[1], argnames[1], argnames[2]
    ))
  }
  # catch wrong cds class
  if (is_gtf(cds)) {
    cds <- S4Vectors::split(cds[cds$type == "CDS"], ~transcript_id)
    if (length(cds) == 0) {
      rlang::abort(sprintf(
        "`%s` do not contain CDS information", argnames[1]
      ))
    }
  }

  if (!is(cds, "GRangesList")) {
    rlang::abort("cds class type is not GRanges GTF or GRangesList")
  }

  cds <- filtereach(cds, ...)
  return(sorteach(cds, exonorder))
}

.featurechecks <- function(hmmer, signalHsmm) {
  # check for correct feature column names for return
  hmmer.header <- c(
    "name", "acc", "bias", "desc", "evalue", "flags", "hindex", "ndom",
    "nincluded", "nregions", "nreported", "pvalue", "score", "taxid",
    "pdb.id", "bitscore", "mlog.evalue"
  )
  if ("default" %in% hmmer) { # convert default into list of headers
    hmmer <- c(
      "desc", "evalue", "nincluded",
      hmmer[hmmer != "default"]
    ) %>% unique()
  }
  if (any(!hmmer %in% hmmer.header)) {
    fields_support <- hmmer %in% hmmer.header
    rlang::warn(sprintf(
      "`%s` are not supported hmmer reported fields",
      hmmer[!fields_support]
    ))
    if (sum(fields_support) == 0) {
      hmmer <- c("desc", "evalue", "nincluded")
      rlang::warn("returning default hmmer fields")
    } else {
      hmmer <- hmmer[fields_support]
    }
  }
  # check for correct field names for return
  signalHsmm.header <- c(
    "sp_probability", "sp_start", "sp_end",
    "struc", "prot", "name", "str_approx"
  )
  if ("default" %in% signalHsmm) { # convert default into list of headers
    signalHsmm <- c(
      "sp_probability",
      signalHsmm[signalHsmm != "default"]
    ) %>% unique()
  }
  if (any(!signalHsmm %in% signalHsmm.header)) {
    fields_support <- signalHsmm %in% signalHsmm.header
    rlang::warn(sprintf(
      "`%s` are not supported signalHsmm reported fields",
      signalHsmm[!fields_support]
    ))
    if (sum(fields_support) == 0) {
      signalHsmm <- "sp_probability"
      rlang::warn("returning default signalHsmm fields")
    } else {
      signalHsmm <- signalHsmm[fields_support]
    }
  }
  return(list(hmmer = hmmer, signalHsmm = signalHsmm))
}

.getSequence <- function(cds, fasta) {
  x <- y <- instop <- NULL
  cdsSeq <- GenomicFeatures::extractTranscriptSeqs(fasta, cds)
  aaSeq <- suppressWarnings(Biostrings::translate(cdsSeq, if.fuzzy.codon = "solve")) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("id") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(y = strsplit(x, split = "")) %>%
    dplyr::mutate(noATG = ifelse(y[[1]] != "M", T, F)) %>%
    dplyr::mutate(instop = ifelse("*" %in% y, T, F)) %>%
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
        paste(y[1:which(y == "*") - 1], collapse = ""),
        x
      )) %>%
      dplyr::mutate(y = strsplit(x, split = "")) %>%
      dplyr::ungroup())

    rlang::warn(sprintf("%s cds entries contain internal stop_codon. These proteins have been truncated", sum(aaSeq$instop)))
    if ("" %in% aaSeq$x) {
      rlang::warn(sprintf("After truncation, %s cds have no coding sequences. These entries were not analyzed", sum(aaSeq$x == "")))
      aaSeq <- aaSeq[aaSeq$x != "", ]
    }
  }
  return(aaSeq)
}



.getdomains <- function(url, curl.opts, seq, id, length, n) {
  hmm <- RCurl::postForm(url, hmmdb = "superfamily", seqdb = NULL, 
                         seq = seq, style = "POST", .opts = curl.opts, .contentEncodeFun = RCurl::curlPercentEncode, 
                         .checkParams = TRUE)
  xml <- XML::xmlParse(hmm)
  family <- XML::xpathSApply(xml, "///family", XML::xpathSApply,  "@*")
  segment <- XML::xpathSApply(xml, "///segments", XML::xpathSApply, "@*")
  data <- rbind(family, segment)
  data <- as.data.frame(t(data), stringsAsFactors = FALSE) %>%
    dplyr::mutate(type = "DOMAIN", begin = as.numeric(start), end = as.numeric(end)) %>%
    dplyr::select(type, description = famdesc, eval = fameval, begin, end) %>%
    dplyr::bind_rows(tibble::tibble(type = "CHAIN", description = id, begin = 1, end = length)) %>%
    dplyr::mutate(entryName = id) %>%
    dplyr::mutate(order = n)
  return(data)
}
