#' Find NMD-causing exons from GTF objects
#' 
#' @description 
#' This function will identify alternative exons that triggers NMD when 
#' spliced or skipped.
#' 
#' @details 
#' Alternative splicing of regulated exons may introduce premature termination 
#' codon which are cleared by nonsense-mediated decay. There are two major class 
#' of such exons: Poison and ORF-maintaining exons. The former triggers NMD when
#' spliced in, hence its name, while the latter triggers NMD when removed from a 
#' transcript.
#' 
#' To identify and classify NMD-causing exons, identifyNMDexons first find
#' alternative segments between NMD-sensitive transcripts and pre-determined 
#' reference transcript of the same gene. The longest, NMD-insensitive coding 
#' transcript of each gene is chosen as the best reference. For each alternative segment, 
#' a modified version of reference RNA is built by inserting exons spliced in NMD transcripts, or
#' by removing exons skipped from NMD transcripts. These modified transcripts are
#' then tested for NMD features and if positive, the alternative exon will be
#' annotated as NMD-causing.
#' 
#' @param x 
#' GRanges object containing query GTF data with CDS information.
#' @param fasta 
#' BSgenome or Biostrings object containing genomic sequence
#' @param NMD.result 
#' Output dataframe from predictNMD() [Optional]. If not provided, or if 
#' data-frame do not contain necessary variables, identifyNMDexons will run
#' predictNMD() first.
#' @param ConsScores
#' Character value of the annotation database to use to calculate mean exon 
#' conservation scores. Database can be a Bioconductor annotation package or
#' one of the available database from AnnotationHub. See
#' \code{\link[GenomicScores]{availableGScores}} function for a list of available 
#' databases. A vector of connservation score for each exon will be appended
#' as a metadata column in the output GRanges object.
#' By default ("none"), exon conservation scores will not be calculated.
#'
#' @return
#' GRanges object containing coordinates of identified NMD-causing exons, 
#' including metadata on its gene of origin, type of AS exon, type of NMD exon
#' and whether it is spliced within CDS
#' 
#' @export
#'
#' @examples
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING SAMPLE DATASET
#' ## ---------------------------------------------------------------------
#' library(BSgenome.Mmusculus.UCSC.mm10)
#' 
#' # build CDS information first
#' out.gtf <- buildCDS(matched_query_gtf, ref_gtf, Mmusculus)
#' 
#' # proceed with identifying NMD exons
#' NMD.exons_1 <- identifyNMDexons(out.gtf, Mmusculus)
#' 
#' # if output of predictNMD() has been generated, this can be provided to
#' # NMD.result argument
#' NMD.output <- predictNMD(out.gtf)
#' NMD.exons_2 <- identifyNMDexons(out.gtf, Mmusculus, NMD.result = NMD.output)
#' 
#' identical(NMD.exons_1, NMD.exons_2)
#' 
#' 
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING TRANSCRIPT ANNOTATION
#' ## ---------------------------------------------------------------------
#' \donttest{
#' library(AnnotationHub)
#'
#' ## Retrieve GRCm38 trancript annotation
#' ah <- AnnotationHub()
#' GRCm38_gtf <- ah[["AH60127"]]
#'
#' ## Identify NMD-causing exons in entire transcriptome
#' NMD.exons.all <- identifyNMDexons(GRCm38_gtf, Mmusculus)
#' }
#' 
identifyNMDexons <- function(x, fasta, 
                             NMD.result = NULL,
                             ConsScores = "none") {
  
  # catch missing args
  .catchargs(c("x", "fasta"), names(as.list(match.call())[-1]))
  
  # retrieve input object names
  argnames <- as.character(match.call())[-1]
  
  # check input objects
  NMD.result <- .identifynmdchecks(x, fasta, NMD.result, argnames) 
  phastGScore <- .GScorecheck(ConsScores)

  # get reference CDS transcript for each gene
  ref <- .getbestref(x, NMD.result)
  
  # run core function and return GRanges of predicted NMD exons
  return(.runidentifynmdexons(x, ref, fasta, NMD.result, phastGScore))

}





.catchargs <- function(mandargs, passed) {
  if (any(!mandargs %in% passed)) {
    rlang::abort(paste(
      "missing values for",
      paste(setdiff(mandargs, passed), collapse = ", ")
    ))
  }
}




.identifynmdchecks <- function(x, fasta, NMD, argnames) {
  
  # check inputs are gtf
  if (!is_gtf(x)) {
    rlang::abort(sprintf(
      "`%s` is not gtf GRanges",
      argnames[1]
    ))
  }
  
  # check if ref have CDS info
  if (!"CDS" %in% S4Vectors::mcols(x)$type) {
    rlang::abort(sprintf(
      "`%s` have missing CDS info. Run buildCDS() to construct CDSs", argnames[1]
    ))
  }
  
  # catch unmatched seqlevels
  if (suppressWarnings(!has_consistentSeqlevels(x, fasta))) {
    rlang::abort(sprintf(
      "`%s` and `%s` has unmatched seqlevel styles. 
Try running: %s <- matchChromosomes(%s, %s)",
      argnames[1], argnames[2], argnames[1], argnames[1], argnames[2]
    ))
  }
  
  
  
  # check for correct formatting of NMD df and create new one if not
  if(!is.null(NMD)) {
    # check format of NMD df
    if(all(c("transcript", "is_NMD") %in% names(NMD))) {
      if(all(x[x$type=="CDS"]$transcript_id %in% NMD$transcript)){
        return(NMD)
      } else {
        rlang::inform(sprintf("Some mRNAs in `%s` are not in `%s`, re-running predictNMD()",
                              argnames[1], argnames[3]))
      }
    } else {
      rlang::inform("`NMD.result` is of wrong format, running predictNMD()")
    }
  } else{
    rlang::inform("`NMD.result` not provided, running predictNMD()")
  }
  return(suppressMessages(predictNMD(x, progress_bar = FALSE)))
}


.GScorecheck <- function(ConsScore){
  if(ConsScore != "none"){
    
    
    ## try loading GScore package
    phast <- tryCatch(
      {
        library(ConsScore, character.only = T)
        rlang::inform(sprintf("Loaded %s package",
                              ConsScore))
        get(ConsScore)
      },
      error = function(cond){
        GScoresList <- rownames(GenomicScores::availableGScores())
        if(ConsScore %in% GScoresList){
          rlang::inform(sprintf("Retrieving %s scores",
                                ConsScore))
          GenomicScores::getGScores(ConsScore)
        } else {
          rlang::inform(sprintf("%s score database not found. Skipping conservation scoring.",
                                ConsScore))
          return(NULL)
        }
      })
    
  }
  return(phast)
}


.getbestref <- function(x, NMD.result) {
  rlang::inform("Selecting best reference mRNAs")
  # get reference CDS transcript for each gene
  ## get sizes of all CDSs
  cds.sizes <- sum(BiocGenerics::width(S4Vectors::split(x[x$type == "CDS"], 
                                                        ~transcript_id)))
  ## shortlist non NMD transcripts
  NMD.pos <- NMD.result[!NMD.result$is_NMD,]
  
  ## select best reference
  cds.reference <- x %>% 
    as.data.frame() %>% 
    dplyr::select(gene_id, transcript_id) %>% 
    dplyr::left_join(as.data.frame(cds.sizes) %>% tibble::rownames_to_column("transcript_id"),
                     by = "transcript_id") %>% 
    dplyr::filter(transcript_id %in% NMD.pos$transcript) %>% 
    dplyr::group_by(gene_id) %>% 
    dplyr::slice_max(order_by = cds.sizes, n = 1)
  
  return(x[x$transcript_id %in% cds.reference$transcript_id])
}

.runidentifynmdexons <- function(x, ref, fasta, NMD.result, phast) {
  
  rlang::inform("Finding NMD causing exons")
  
  # get AS segments between NMD transcript and reference transcript
  ## shortlist NMD transcripts from GTF
  x.NMD <- x[x$transcript_id %in% NMD.result[NMD.result$is_NMD,]$transcript]
  
  # shortlist NMD transcripts from genes containing reference mRNAs
  x.NMD <- x.NMD[x.NMD$gene_id %in% ref$gene_id]
  
  ## get AS segments and annotate its splicing nature
  AS.exons <- labelSplicedSegment(c(x.NMD, ref))
  AS.exons$splice <- ifelse(AS.exons$transcript_id %in% ref$transcript_id,
                            "skipped", "spliced")
  
  # return if no AS exons were found
  if(length(AS.exons) == 0) {
    rlang::warn("No alternatively spliced exons found")
    return(NULL)
  }
  
  # simplify df by removing redundant exons
  ## in cases where ref and NMD tx contain said exon, the skipped form will be retained
  AS.exons <- AS.exons %>%
    as.data.frame() %>% 
    dplyr::mutate(coord = paste0(seqnames, "_", start, "_", 
                                 end, "_", strand, "_", AStype)) %>% 
    dplyr::arrange(splice) %>%
    dplyr::select(gene_id, transcript_id, coord, splice) %>% 
    dplyr::distinct(coord, .keep_all = TRUE)
  
  # recreate hypothetical tx by inserting/removing AS segments 
  ## create grl of ref trancripts
  ref.grl <- S4Vectors::split(ref[ref$type == "exon"], ~gene_id)
  mod.tx <- tibble::tibble(
    "seqnames" = as.character(),
    "start" = as.integer(),
    "end" = as.integer(),
    "strand" = as.character(),
    "coord" = as.character()
  )
  
  ## Add exons to ref transcripts
  if("spliced" %in% AS.exons$splice) {
    toadd <- dplyr::filter(AS.exons, splice == "spliced")
    addedtx <- ref.grl[toadd$gene_id]
    names(addedtx) <- toadd$coord
    exons.toadd <- toadd["coord"] %>%
      tidyr::separate(coord, c("seqnames", "start", "end", "strand", "AStype"),
                      sep = "_") %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    addedtx <- addExonstoTx(addedtx, exons.toadd) %>% 
      as.data.frame() %>% 
      dplyr::select(seqnames, start, end, strand, coord = group_name)
    mod.tx <- dplyr::bind_rows(mod.tx, addedtx)
  }
  
  
  ## Remove exons from ref transcripts
  if("skipped" %in% AS.exons$splice) {
    toremove <- dplyr::filter(AS.exons, splice == "skipped")
    removedtx <- ref.grl[toremove$gene_id]
    names(removedtx) <- toremove$coord
    exons.toremove <- toremove["coord"] %>%
      tidyr::separate(coord, c("seqnames", "start", "end", "strand", "AStype"),
                      sep = "_") %>%
      GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    removedtx <- removeExonsfromTx(removedtx, exons.toremove) %>% 
      as.data.frame() %>% 
      dplyr::select(seqnames, start, end, strand, coord = group_name)
    mod.tx <- dplyr::bind_rows(mod.tx, removedtx)
  }
  
  
  # cleanup mod.tx by adding gene_id info and create GRanges
  mod.tx <- mod.tx %>% 
    dplyr::left_join(AS.exons %>% dplyr::select(coord, gene_id), by = "coord") %>% 
    dplyr::mutate(transcript_id = coord, type = "exon") %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  # test NMD again
  mod.tx <- suppressMessages(buildCDS(mod.tx, ref, fasta))
  mod.NMD <- suppressMessages(predictNMD(mod.tx, progress_bar = FALSE))
  
  # return if none of the reconstructed transcripts are NMD sensitive
  if(sum(mod.NMD$is_NMD) == 0) {
    rlang::warn("None of the alternative exons are NMD-causing")
    return(NULL)
  }
  
  # Report NMD exons
  AS.exons <- AS.exons %>% 
    dplyr::left_join(mod.NMD %>% dplyr::select(transcript, is_NMD), 
                     by = c("coord"="transcript")) %>% 
    dplyr::mutate(NMDtype = ifelse(splice == "skipped", "Repressive", "Stimulating")) %>% 
    dplyr::select(coord, gene_id, NMDtype, is_NMD) %>% 
    tidyr::separate(coord, c("seqnames", "start", "end", "strand", "AStype"),
                    sep = "_") %>% 
    dplyr::filter(is_NMD) %>% 
    dplyr::select(seqnames:strand, gene_id, AStype, NMDtype) %>% 
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  AS.exons$is_NMD <- NULL
  
  # Annotate location of NMD exons relative to reference CDS
  ## Useful to know if exons are 3'UTR introns etc
  ref.cds.grl <- S4Vectors::split(ref[ref$type == "CDS"], ~gene_id)
  AS.exons$within.CDS <- IRanges::overlapsAny(AS.exons, range(ref.cds.grl))
  
  ## run conservation scoring
  if(typeof(phast) == "S4"){
    rlang::inform("Quantifying exon conservation scores")
    AS.exons <- tryCatch(
      {
        AS.exons <- GenomicScores::gscores(phast, AS.exons)
        BiocGenerics::colnames(S4Vectors::mcols(AS.exons))[5] <- phast@data_pkgname
        return(AS.exons)
      },
      error = function(cond){
        rlang::warn("Unable to quantify exon conservation scores. Check if annotation package matches the genome used")
        return(AS.exons)
      }
    )
    

  }
  
  return(AS.exons)
}

