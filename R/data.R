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

#' Dataframe containing a list of matched query and reference IDs
#'
#' Dataframe with 2 columns, each containing a list of transcript_id from
#' query_gtf and ref_gtf dataset respectively
#'
#' @format A dataframe with 8 rows
#' \describe{
#'   \item{transcript_id}{List of transcript_id from query_gtf}
#'   \item{ref_transcript_id}{List of transcript_id from ref_gtf}
#'   ...
#' }
"q2r"

#' Dataframe containing a list of matched query and reference IDs with coverage values
#'
#' Dataframe with 3 columns containing a list of transcript_id from
#' query_gtf, a list of matched ref_transcript_id and the percent coverage between
#' query and reference transcripts
#'
#' @format A dataframe with 8 rows
#' \describe{
#'   \item{transcript_id}{List of transcript_id from query_gtf}
#'   \item{ref_transcript_id}{List of transcript_id from ref_gtf}
#'   \item{coverage}{List of percent coverage values between query and ref transcripts}
#'   ...
#' }
"q2rcovs"
