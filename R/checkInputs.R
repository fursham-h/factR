#' Check input files for workflow
#'
#' @param query GRanges object containing imported query GTF data 
#' @param ref GRanges object containing reference GTF data
#' @param fasta BSgenome or Biostrings object containing genomic sequence
#'
#' @export
#' @importFrom dplyr %>%
#' @importFrom IRanges %over%
#' @importFrom IRanges %within%
#' @importFrom data.table :=
#' @importFrom graphics frame
#' @importFrom methods as is new slot slotNames
#' @importFrom stats dist end start
#' @importFrom utils head tail
#' @importFrom rlang .data
checkInputs <- function(query, ref, fasta) {
  
  argnames <- as.character(match.call())[-1]

  # check seqlevels
  message("Checking inputs for consistent seqlevels...")
  
  seqreport <- unlist(lapply(list(query,ref), function(x){
    return(sum(!GenomeInfoDb::seqlevels(x) %in% GenomeInfoDb::seqlevels(fasta)))
  }))
  
  if (sum(seqreport) == 0){
    message('\tOK...')
  } else {
    lapply(seq_along(seqreport), function(x) {
      if (seqreport[x] > 0){
        message(sprintf("\t%s contain inconsistent seqlevels",
                        argnames[x]))
        message(sprintf("\tTry running: %s <- matchSeqLevels(%s, %s)",
                argnames[x], argnames[x], argnames[3]))
      }
    })
  }

  # check metadata
  message("\nChecking inputs for transcript_id and gene_id attributes...")
  attreport <- lapply(list(query,ref), function(x){
    return(!c('gene_id', 'transcript_id') %in% names(S4Vectors::mcols(x)))
  })
  if (sum(unlist(attreport)) == 0) {
    message('\tOK...')
    # check gene_id
    message("\nChecking query for mismatched gene_ids...")
    inIDs <- unique(query$gene_id)
    refIDs <- unique(ref$gene_id)
    unmatched <- sum(!inIDs %in% refIDs)
    if (unmatched == 0) {
      message("\tOK...")
    } else {
      message(sprintf(
        "\t%s out of %s query gene_ids are not found in reference",
        unmatched, length(inIDs)
      ))
      message(sprintf("\tTry running: %s <- matchGeneIDs(%s, %s)", 
                      argnames[1], argnames[1], argnames[2]))
    }
  } else {
    errorindex <- unlist(lapply(seq_along(attreport), function(x) {
      if (sum(attreport[[x]]) > 0) {
        missing <- c('gene_id', 'transcript_id')[attreport[[x]]]
        message(sprintf("\t%s have missing %s attribute(s)",
                        argnames[x], paste(missing, collapse = ',')))
        return(x)
      }
    }))
    
    message(sprintf('\tPlease check %s object(s) and rerun checkInput',
                    paste(argnames[errorindex], collapse = ',')))
  }
}
