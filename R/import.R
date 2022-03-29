#' Import GTF file into R
#'
#' @description
#' This function loads GTF files into R and converts it into
#'
#'
#' a wrapper to rtracklayer::import() function to conveniently
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
    if(file.exists(con)){
        infile <- try(rtracklayer::import(con, format = "GTF"), silent = TRUE)
        if (is_gtf(infile)) {
            return(infile)
        } else {
            rlang::abort(paste(con, " is not in GTF format"))
        }
    } else {
        rlang::abort(paste("File ", con, " does not exists"))
    }
    
}


#' Import FASTA file into R
#'
#' @description
#' This function is a wrapper to Biostrings::readDNAStringSet() 
#' function to import FASTA genome sequence file and simultaneously 
#' convert long chromosome names 
#' (e.g. 1 dna:chromosome chromosome:GRCm38:1:1:195471971:1 REF) 
#' to short names (e.g. 1)
#'
#' @param con
#' Path to FASTA file
#'
#' @return
#' Imported DNAStringSet object
#'
#' @author Fursham Hamid
#' @export
importFASTA <- function(con) {
    infile <- tryCatch(
        {
            Biostrings::readDNAStringSet(con)
        },
        error = function(x) {
            rlang::abort(paste(con, " not found or not in FASTA format"))
        }
    )
    names(infile) <- stringr::str_split(seqlevels(infile), " ") %>% 
        purrr::map_chr(`[`, 1)
    return(infile)
}
