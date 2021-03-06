% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{query_gtf}
\alias{query_gtf}
\title{Imported GTF file containing 4 transcript entries of the same gene}
\format{
A GRanges object with 56 ranges and 3 metadata columns:
\describe{
  \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
  and its exons}
  \item{type}{Entry type; transcript or exon}
  \item{transcript_id}{Name or ID given to transcripts}
  \item{gene_id}{Name or ID given to gene origin of transcripts}
  ...
}
}
\source{
\url{http://www.ensembl.org/}
}
\usage{
data(query_gtf)
}
\description{
A dataset containing coordinates of transcript and exons from 4 transcripts 
of mouse Ptbp1. Transcript names and gene IDs have been modified to 
demonstrate de novo origin of GTF
}
\keyword{datasets}
