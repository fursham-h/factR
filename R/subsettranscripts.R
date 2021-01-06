#' Shortlist GTF GRanges object for new transcripts
#'
#' @description
#' `subsetNewTranscripts()` will retain transcripts in `query` that are distinct
#' from those in `ref`
#'
#' @details
#' `subsetNewTranscripts()` will compare query and reference GTF GRanges and return
#' query transcripts with different exon structures from reference transcripts.
#' Transcriptome assemblers may sometime extend 5' and 3' ends of known transcripts
#' based on experimental data. These annotated transcripts can be removed by inputting
#' "intron" to the refine.by argument. This will further compare and remove transcripts
#' of identical intron structures. Alternatively, transcripts with unique CDS coordinates
#' can be selected by typing "cds" to the refine.by argument.
#'
#'
#' @param query
#' GRanges object containing query GTF data.
#' @param ref
#' GRanges object containing reference GTF data.
#' @param refine.by
#' Whether to refine the selection process by removing query transcripts with similar
#' introns or CDS structure to reference. Default input is "none", and can be changed
#' to "intron" or "cds" respectively.
#'
#' @return
#' Filtered GRanges GTF object
#'
#' @export
#'
#' @author Fursham Hamid
#' @examples
#' subsetNewTranscripts(matched_query_gtf, ref_gtf)
subsetNewTranscripts <- function(query, ref, refine.by = c("none", "intron", "cds")) {
    
    # catch missing args
    mandargs <- c("query", "ref")
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
    .checkinputs(query, ref, argnames)
    
    # run subsetting and return new transcripts
    return(.subsetTranscripts(query, ref, refine.by[1]))
}

.checkinputs <- function(query, ref, argnames) {
    # check inputs are gtf
    if (any(!is_gtf(query, ref))) {
        rlang::abort(sprintf(
            "%s is/are not gtf GRanges",
            paste(argnames[seq_len(2)][!is_gtf(query, ref)], collapse = ",")
        ))
    }
    
    
    # catch unmatched seqlevels
    if (suppressWarnings(!has_consistentSeqlevels(query, ref))) {
        rlang::abort(sprintf(
            "`%s` and `%s` has unmatched seqlevel styles. 
Try running: %s <- matchChromosomes(%s, %s)",
            argnames[1], argnames[2], argnames[1], argnames[1], argnames[2]
        ))
    } else {
        outmsg <- has_consistentSeqlevels(query, ref)
    }
}

.subsetTranscripts <- function(query, ref, by) {
    
    # Filter transcripts by identical exon structure:
    # create exons by transcripts
    query_exons <- S4Vectors::split(query[query$type == "exon"], ~transcript_id)
    ref_exons <- S4Vectors::split(ref[ref$type == "exon"], ~transcript_id)
    
    # search for exact query and ref matches
    fulloverlap <- GenomicRanges::findOverlaps(query_exons, ref_exons,
                                               type = "equal", select = "first"
    )
    
    # return transcripts that are not found in reference
    query <- query[!query$transcript_id %in% names(query_exons)[!is.na(fulloverlap)]]
    
    # Refine list By identical introns:
    if (by == "intron") {
        # create exons by transcripts
        query_exons <- S4Vectors::split(query[query$type == "exon"], ~transcript_id)
        ref_exons <- S4Vectors::split(ref[ref$type == "exon"], ~transcript_id)
        
        # convert exon coord to intron coord
        query_introns <- GenomicRanges::psetdiff(BiocGenerics::unlist(range(query_exons)), query_exons)
        ref_introns <- GenomicRanges::psetdiff(BiocGenerics::unlist(range(ref_exons)), ref_exons)
        
        # search for exact query and ref matches
        fulloverlap_intron <- GenomicRanges::findOverlaps(query_introns, ref_introns,
                                                          type = "equal", select = "first"
        )
        # return transcripts that are not found in reference
        query <- query[!query$transcript_id %in% names(query_exons)[!is.na(fulloverlap_intron)]]
    }
    
    # Refine list By identical CDS:
    else if (by == "cds") {
        # create CDS by transcripts
        query_CDS <- S4Vectors::split(query[query$type == "CDS"], ~transcript_id)
        ref_CDS <- S4Vectors::split(ref[ref$type == "CDS"], ~transcript_id)
        
        # search for exact query and ref matches
        fullCDSoverlap <- GenomicRanges::findOverlaps(query_CDS, ref_CDS,
                                                      type = "equal", select = "first"
        )
        
        # return coding transcripts that are not found in reference
        query <- query[query$transcript_id %in% names(query_CDS)]
        query <- query[!query$transcript_id %in% names(query_CDS)[!is.na(fullCDSoverlap)]]
    }
    return(query)
}
