% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{query_exons}
\alias{query_exons}
\title{GRangeList of exons from 4 transcripts entries from query_gtf}
\format{
A GRangesList object with 4 elements:
\describe{
  \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
  and its exons}
  \item{type}{Entry type; transcript or exon}
  \item{transcript_id}{ID given to transcripts}
  \item{gene_id}{Matched gene_id}
  \item{old_gene_id}{Original gene_id}
  \item{match_level}{Level of matching performed}
  \item{gene_name}{Name of gene}
  ...
}
}
\usage{
data(query_exons)
}
\description{
A dataset containing coordinates of exons from 4 transcripts of
mouse Ptbp1. Transcript names and gene IDs have been modified
}
\keyword{datasets}
