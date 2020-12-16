#' Wrapper to import GTF file
#'
#' @param con 
#' Path to GTF file
#'
#' @return
#' Imported GenomicRanges object in GTF format
#' 
#' @export
#'
importGTF <- function(con) {
  infile <- try(rtracklayer::import(con, format = "GTF"), silent = T)
  if(is_gtf(infile)) {
    return(infile)
  } else {
    rlang::abort(paste(con, " is not in GTF format"))
  }
}