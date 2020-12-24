#' Wrapper to import GTF file
#'
#' @description 
#' This function is a wrapper to rtracklayer::import() function to conveniently
#' import GTF file into R as a GenomicRanges object.
#' 
#' @param con
#' Path to GTF file
#'
#' @return
#' Imported GenomicRanges object in GTF format
#'
#' @author Fursham Hamid
#' @export
importGTF <- function(con) {
  infile <- try(rtracklayer::import(con, format = "GTF"), silent = T)
  if (is_gtf(infile)) {
    return(infile)
  } else {
    rlang::abort(paste(con, " is not in GTF format"))
  }
}
