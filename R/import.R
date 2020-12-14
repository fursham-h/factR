importGTF <- function(con) {
  infile <- rtracklayer::import(con)
  if(is_gtf(infile)) {
    return(infile)
  } else {
    rlang::abort(paste(con, " is not in GTF format"))
  }
}