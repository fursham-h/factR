#' Check input files for workflow
#'
#' @param input GRanges object containing imported query GTF data 
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
checkInputs <- function(input, ref, fasta) {

  # check seqlevels
  message("Report on SeqLevels:")
  styleI <- GenomeInfoDb::seqlevelsStyle(input)[1]
  styleR <- GenomeInfoDb::seqlevelsStyle(ref)[1]
  styleF <- GenomeInfoDb::seqlevelsStyle(fasta)[1]
  levels <- c(styleI, styleR, styleF)
  names <- c("input", "ref", "fasta")

  if (length(unique(levels)) == 1) {
    message("\tSeqlevels on all objects are identical")
  } else {
    lapply(unique(levels), function(x) {
      objs <- paste(names[which(levels == x)], collapse = ", ")
      message(sprintf("\t%s objects have Seqlevel style of %s", objs, x))
    })
    message("\tUse function matchSeqLevels")
  }

  # check metadata
  message("Report on input metadata")
  inputIDs <- c('gene_id', 'transcript_id') %in% names(mcols(input))
  refIDs <- c('gene_id', 'transcript_id') %in% names(mcols(ref))
  
  
  # check gene_id
  message("Report on gene_ids:")
  inIDs <- unique(input$gene_id)
  refIDs <- unique(ref$gene_id)
  unmatched <- sum(!inIDs %in% refIDs)
  if (unmatched == 0) {
    message("\tAll query gene_ids matched to reference")
  } else {
    message(sprintf(
      "\t%s out of %s query gene_ids are not matched to reference",
      unmatched, length(inIDs)
    ))
    message("\tUse function matchGeneIDs")
  }
  
}
