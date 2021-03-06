% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eachfunctions.R
\name{filtereach}
\alias{filtereach}
\title{Internally filter each element of a GenomicRangesList}
\usage{
filtereach(x, ...)
}
\arguments{
\item{x}{GRangesList object}

\item{...}{Logical conditions to filter each element in the GRanges by. 
Multiple conditions can be provided as comma-delimited inputs}
}
\value{
Filtered GRangesList object
}
\description{
Internally filter each element of a GenomicRangesList
}
\examples{
# Load dataset
data(query_exons)

# select first element of each GRangesList item
filtereach(query_exons, dplyr::row_number() == 1)
}
\author{
Fursham Hamid
}
