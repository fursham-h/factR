#' Imported GTF file containing 4 transcript entries of the same gene
#'
#' A dataset containing coordinates of transcript and exons from 4 transcripts of
#' mouse Ptbp1. Transcript names and gene IDs have been modified to demonstrate
#' de novo origin of GTF
#'
#' @format A GRanges object with 56 ranges and 3 metadata columns:
#' \describe{
#'   \item{element}{Exon start and end coordinates from 4 transcripts}
#'   ...
#' }
#' @source \url{http://www.ensembl.org/}
"query_gtf"


#' Seqlevels and gene_id matched query data
#'
#' query_gtf data which have been corrected for its seqlevels and gene_ids
#'
#' @format A GRanges object with 56 ranges and 6 metadata columns:
#' \describe{
#'   \item{element}{Exon start and end coordinates from 4 transcripts}
#'   ...
#' }
"matched_query_gtf"

#' Query data containing CDS information
#'
#' matched_query_gtf data that has undergone buildCDS function and containing
#' CDS features
#'
#' @format A GRanges object with 105 ranges and 8 metadata columns:
#' \describe{
#'   \item{element}{Exon start and end coordinates from 4 transcripts}
#'   \item{element}{CDS start and end coordinates from 4 transcripts}
#'   ...
#' }
"new_query_gtf"


#' Exons from 4 transcripts entries of the same gene
#'
#' A dataset containing coordinates of exons from 4 transcripts of
#' mouse Ptbp1. Transcript names and gene IDs have been modified
#'
#' @format A GRangesList object with 4 elements:
#' \describe{
#'   \item{element}{Exon start and end coordinates of each transcript}
#'   ...
#' }
#' @source \url{http://www.ensembl.org/}
"query_exons"




#' CDS from 4 transcripts entries of the same gene
#'
#' A dataset containing coordinates of CDS from 4 transcripts of
#' mouse Ptbp1. Transcript names and gene IDs have been modified
#'
#' @format A GRangesList object with 4 elements:
#' \describe{
#'   \item{element}{CDS start and end coordinates of each transcript}
#'   ...
#' }
#' @source \url{http://www.ensembl.org/}
"query_cds"

#' Imported GTF file containing 2 reference transcript entries of the same gene
#'
#' A dataset containing coordinates of transcript and exons from 2 reference
#' transcripts of mouse Ptbp1.
#'
#' @format A GRanges object with 64 ranges and 5 metadata columns:
#' \describe{
#'   \item{element}{Exon start and end coordinates from 2 reference transcripts}
#'   ...
#' }
#' @source \url{https://www.gencodegenes.org/}
"ref_gtf"

#' Exons from 2 reference transcripts entries of the same gene
#'
#' A dataset containing coordinates of exons from 2 reference transcripts of
#' mouse Ptbp1.
#'
#' @format A GRangesList object with 2 elements:
#' \describe{
#'   \item{element}{Exon start and end coordinates of each transcript}
#'   ...
#' }
#' @source \url{https://www.gencodegenes.org/}
"ref_exons"

#' CDS from 2 reference transcripts entries of the same gene
#'
#' A dataset containing coordinates of CDS from 2 reference transcripts of
#' mouse Ptbp1.
#'
#' @format A GRangesList object with 2 elements:
#' \describe{
#'   \item{element}{CDS start and end coordinates of each transcript}
#'   ...
#' }
#' @source \url{https://www.gencodegenes.org/}
"ref_cds"
