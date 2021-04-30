#' Predict NMD sensitivity on mRNA transcripts
#'
#' @param x
#' Can be a GRanges object containing exon and CDS transcript features in GTF
#' format.
#'
#' Can be a GRangesList object containing exon features for a list of 
#' transcripts.If so, `cds` argument have to be provided.
#'
#' Can be a GRanges object containing exon features for a transcript. 
#' If so, `cds` argument have to be provided.
#'
#' @param ...
#' Logical conditions to pass to dplyr::filter to subset transcripts for 
#' analysis. Variables are metadata information found in `x` and multiple 
#' conditions can be provided delimited by comma. 
#' Example: transcript_id == "transcript1"
#'
#' @param cds
#' If `x` is a GRangesList object, `cds` has to be a GRangesList containing 
#' CDS features for the list of transcripts in `x`. List names in `x` and 
#' `cds` have to match.
#'
#' If `x` is a GRanges object, `cds` has to be a GRanges containing CDS 
#' features for the transcript in `x`.
#'
#' @param NMD_threshold
#' Minimum distance of stop_codon to last exon junction 
#' (EJ) which triggers NMD. Default = 50bp
#' 
#' @param progress_bar
#' Whether to display progress 
#' Default = TRUE
#'
#' @return
#' Dataframe with prediction of NMD sensitivity and NMD features:
#'
#' is_NMD: logical value in prediciting transcript sensitivity to NMD
#'
#' stop_to_lastEJ: Integer value of the number of bases between the first
#' base of the stop_codon to the last base of EJ. A positive value indicates 
#' that the last EJ is downstream of the stop_codon.
#'
#' num_of_down_EJs: Number of EJs downstream of the stop_codon.
#'
#' `3_UTR_length`: Length of 3' UTR
#'
#' @export
#' @author Fursham Hamid
#'
#' @examples
#'
#' ## ---------------------------------------------------------------------
#' ## EXAMPLE USING SAMPLE DATASET
#' ## ---------------------------------------------------------------------
#' 
#' # Load datasets
#' data(new_query_gtf, query_exons, query_cds)
#' 
#' ## Using GTF GRanges as input
#' predictNMD(new_query_gtf)
#'
#' ### Transcripts for analysis can be subsetted using logical conditions
#' predictNMD(new_query_gtf, transcript_id == "transcript1")
#' predictNMD(new_query_gtf, 
#' transcript_id %in% c("transcript1", "transcript3"))
#'
#'
#' ## Using exon and CDS GRangesLists as input
#' predictNMD(query_exons, cds = query_cds)
#' predictNMD(query_exons, cds = query_cds, transcript_id == "transcript3")
#'
#'
#' ## Using exon and CDS GRanges as input
#' predictNMD(query_exons[[3]], cds = query_cds[[3]])
#'
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
#' ## Run tool on specific gene family
#' predictNMD(GRCm38_gtf, gene_name == "Ptbp1")
#' }
#'
predictNMD <- function(x, ..., cds = NULL, NMD_threshold = 50, 
                       progress_bar = TRUE) {

    # catch missing args
    mandargs <- c("x")
    passed <- names(as.list(match.call())[-1])
    if (any(!mandargs %in% passed)) {
        rlang::abort(paste(
            "missing values for",
            paste(setdiff(mandargs, passed), collapse = ", ")
        ))
    }
    # define global variable
    exonorder <- NULL

    # retrieve input object names
    argnames <- as.character(match.call())[-1]

    # check if exons is a gtf or both exons and cds are GR or GRlist
    if (is_gtf(x)) {
        exons <- sorteach(S4Vectors::split(x[x$type == "exon"], 
                                           ~transcript_id), exonorder)
        cds <- S4Vectors::split(x[x$type == "CDS"], ~transcript_id)
    } else if (all(is(x) %in% is(cds))) {
        if (is(x, "GRanges")) {
            exons <- sorteach(GenomicRanges::GRangesList("transcript" = x), 
                              exonorder)
            cds <- GenomicRanges::GRangesList("transcript" = cds)
        } else if (is(x, "GRangesList")) {
            exons <- sorteach(x, exonorder)
        } else {
            errorobj <- paste(unique(c(is(x)[1], is(cds)[2])), collapse = ", ")
            rlang::abort(sprintf(
                "input object types %s not supported", errorobj
            ))
        }
    } else {
        txtype <- is(x)[1]
        cdstype <- is(cds)[1]
        rlang::abort(sprintf(
            "`%s` is type %s but `%s` is type %s",
            argnames[2], cdstype, argnames[1], txtype
        ))
    }

    # catch unmatched seqlevels
    if (suppressWarnings(!has_consistentSeqlevels(exons, cds))) {
        rlang::abort("exons and cds have inconsistent seqlevels")
    }

    # subset list for `which` and check if all exons have a cds entry
    totest <- .prepTotest(exons, cds, argnames, ...)

    # run NMD analysis and return result
    return(.testNMD(exons[totest], cds[totest], NMD_threshold, progress_bar))
}



.testNMD <- function(x, y, threshold, progress_bar) {

    # define global variables
    exonorder <- strand <- start1 <- end1 <- group <- NULL
    seqnames <- newstart <- newend <- width <- NULL
    strand...7 <- start...1 <- start...4 <- end...5 <- end...12 <- NULL
    start...11 <- group...1 <- group_name...2 <- seqnames...3 <- NULL

    out <- tibble::tibble(
        "transcript" = as.character(),
        "stop_to_lastEJ" = as.double(),
        "num_of_downEJs" = as.integer(),
        # "stop_to_downEJs" = as.character(),
        "3'UTR_length" = as.double(),
        "is_NMD" = as.logical()
    )

    toStopRange <- suppressMessages(
        dplyr::bind_cols(as.data.frame(range(x)), as.data.frame(range(y))) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(newstart = ifelse(strand...7 == "-", 
                                        start...11, start...4)) %>%
        dplyr::mutate(newend = ifelse(
            strand...7 == "-", end...5, end...12)) %>%
        dplyr::select(group = group...1, group_name = group_name...2, 
                      seqnames = seqnames...3, start = newstart, 
                      end = newend, strand = strand...7) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE))

    toStopWidth <- sum(BiocGenerics::width(
        GenomicRanges::pintersect(x, toStopRange)))
    EJtoStop <- cumsum(BiocGenerics::width(x)) - toStopWidth
    
    # switch off progress_bar if requested
    if(!progress_bar) {
        pbo <- pbapply::pboptions(type = "none")
        on.exit(pbapply::pboptions(pbo), add = TRUE)
    }

    out <- dplyr::bind_rows(out, pbapply::pblapply(EJtoStop, function(x) {
        id <- ifelse(!is.null(names(x)), names(x)[1], "transcript")
        x <- sort(x, decreasing = TRUE)
        threeUTR <- x[1]
        dist_to_last <- x[2]
        is_NMD <- ifelse(dist_to_last > threshold, TRUE, FALSE)
        dist_to_eachEJ <- rev(x[-1][x[-1] > 0])



        return(tibble::tibble(
            "transcript" = id,
            "stop_to_lastEJ" = dist_to_last,
            "num_of_downEJs" = length(dist_to_eachEJ),
            # "stop_to_downEJs" = paste(dist_to_eachEJ, collapse = ","),
            "3'UTR_length" = threeUTR,
            "is_NMD" = is_NMD
        ))
    }) %>% dplyr::bind_rows() %>%
        tidyr::replace_na(list(is_NMD = FALSE)))
    return(out)
}


.prepTotest <- function(x, y, arg, ...) {
    totest <- names(x) # prepare vector with names for testing
    if (!missing(...)) {
        totest <- tryCatch(
            names(filtereach(x, ...)),
            error = function(e) {
                rlang::abort(sprintf(
                    "Metadata in ... not found in `%s`",
                    arg[1]
                ))
            }
        )
        if (length(totest) == 0) {
            return(NULL)
        }
    }
    txwithcds <- intersect(totest, names(y)) # subset transcripts with cds info
    if (length(txwithcds) == 0) {
        rlang::abort("All transcripts have no CDS information")
    }
    if (length(txwithcds) < length(totest)) {
        totest <- txwithcds
    }
    rlang::inform(sprintf(
        "Predicting NMD sensitivities for %s mRNAs",
        length(totest)
    ))
    return(totest)
}
