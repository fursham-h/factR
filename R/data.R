#' Imported GTF file containing 4 transcript entries of the same gene
#'
#' A dataset containing coordinates of transcript and exons from 4 transcripts 
#' of mouse Ptbp1. Transcript names and gene IDs have been modified to 
#' demonstrate de novo origin of GTF
#'
#' @format A GRanges object with 56 ranges and 3 metadata columns:
#' \describe{
#'   \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
#'   and its exons}
#'   \item{type}{Entry type; transcript or exon}
#'   \item{transcript_id}{Name or ID given to transcripts}
#'   \item{gene_id}{Name or ID given to gene origin of transcripts}
#'   ...
#' }
#' @source \url{http://www.ensembl.org/}
#' @usage data(query_gtf)
"query_gtf"

#' Chromosome matched version of "query_gtf"
#'
#' query_gtf data which have been corrected for its seqlevels
#'
#' @format A GRanges object with 56 ranges and 3 metadata columns:
#' \describe{
#'   \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
#'   and its exons}
#'   \item{type}{Entry type; transcript or exon}
#'   \item{transcript_id}{ID given to transcripts}
#'   \item{gene_id}{ID given to gene origin of transcripts}
#'   ...
#' }
#' @usage data(chrom_matched_query_gtf)
"chrom_matched_query_gtf"

#' Seqlevels and gene_id matched query data
#'
#' query_gtf data which have been corrected for its seqlevels and gene_ids
#'
#' @format A GRanges object with 56 ranges and 6 metadata columns:
#' \describe{
#'   \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
#'   and its exons}
#'   \item{type}{Entry type; transcript or exon}
#'   \item{transcript_id}{ID given to transcripts}
#'   \item{gene_id}{Matched gene_id}
#'   \item{old_gene_id}{Original gene_id}
#'   \item{match_level}{Level of matching performed}
#'   \item{gene_name}{Name of gene}
#'   ...
#' }
#' @usage data(matched_query_gtf)
"matched_query_gtf"

#' Query data containing CDS information
#'
#' matched_query_gtf data that has undergone buildCDS function and containing
#' CDS features
#'
#' @format A GRanges object with 105 ranges and 7 metadata columns:
#' \describe{
#'   \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
#'   and its exons}
#'   \item{type}{Entry type; transcript or exon}
#'   \item{transcript_id}{ID given to transcripts}
#'   \item{gene_id}{Matched gene_id}
#'   \item{old_gene_id}{Original gene_id}
#'   \item{match_level}{Level of matching performed}
#'   \item{gene_name}{Name of gene}
#'   \item{phase}{Phase of open-reading frame}
#'   ...
#' }
#' @usage data(new_query_gtf)
"new_query_gtf"


#' GRangeList of exons from 4 transcripts entries from query_gtf
#'
#' A dataset containing coordinates of exons from 4 transcripts of
#' mouse Ptbp1. Transcript names and gene IDs have been modified
#'
#' @format A GRangesList object with 4 elements:
#' \describe{
#'   \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
#'   and its exons}
#'   \item{type}{Entry type; transcript or exon}
#'   \item{transcript_id}{ID given to transcripts}
#'   \item{gene_id}{Matched gene_id}
#'   \item{old_gene_id}{Original gene_id}
#'   \item{match_level}{Level of matching performed}
#'   \item{gene_name}{Name of gene}
#'   ...
#' }
#' @usage data(query_exons)
"query_exons"




#' CDS from 4 transcripts entries of the same gene
#'
#' A dataset containing coordinates of CDS from 4 transcripts of
#' mouse Ptbp1. Transcript names and gene IDs have been modified
#'
#' @format A GRangesList object with 4 elements:
#' \describe{
#'   \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
#'   and its exons}
#'   \item{type}{Entry type; transcript or exon}
#'   \item{transcript_id}{ID given to transcripts}
#'   \item{phase}{Phase of open-reading frame}
#'   \item{built_from}{Method by which CDS was built}
#'   ...
#' }
#' @usage data(query_cds)
"query_cds"

#' Imported GTF file containing 2 reference transcript entries of the same gene
#'
#' A dataset containing coordinates of transcript and exons from 2 reference
#' transcripts of mouse Ptbp1.
#'
#' @format A GRanges object with 64 ranges and 5 metadata columns:
#' \describe{
#'   \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
#'   and its exons}
#'   \item{type}{Entry type; transcript or exon}
#'   \item{phase}{Phase of open-reading frame}
#'   \item{gene_id}{Matched gene_id}
#'   \item{gene_name}{Name of gene}
#'   \item{transcript_id}{ID given to transcripts}
#'   ...
#' }
#' @source \url{https://www.gencodegenes.org/}
#' @usage data(ref_gtf)
"ref_gtf"

#' Exons from 2 reference transcripts entries of the same gene
#'
#' A dataset containing coordinates of exons from 2 reference transcripts of
#' mouse Ptbp1.
#'
#' @format A GRangesList object with 2 elements:
#' \describe{
#'   \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
#'   and its exons}
#'   \item{type}{Entry type; transcript or exon}
#'   \item{phase}{Phase of open-reading frame}
#'   \item{gene_id}{Matched gene_id}
#'   \item{gene_name}{Name of gene}
#'   \item{transcript_id}{ID given to transcripts}
#'   ...
#' }
#' @source \url{https://www.gencodegenes.org/}
#' @usage data(ref_exons)
"ref_exons"

#' CDS from 2 reference transcripts entries of the same gene
#'
#' A dataset containing coordinates of CDS from 2 reference transcripts of
#' mouse Ptbp1.
#'
#' @format A GRangesList object with 2 elements:
#' \describe{
#'   \item{ranges}{Chromosome, start, end, and strand info of 4 transcripts
#'   and its exons}
#'   \item{type}{Entry type; transcript or exon}
#'   \item{phase}{Phase of open-reading frame}
#'   \item{gene_id}{Matched gene_id}
#'   \item{gene_name}{Name of gene}
#'   \item{transcript_id}{ID given to transcripts}
#'   ...
#' }
#' @source \url{https://www.gencodegenes.org/}
#' @usage data(ref_cds)
"ref_cds"

#' Example output of predictDomains()
#'
#' Output dataframe from predictDomains() function.
#'
#' @format A data.frame with 14880 rows and 5 columns:
#' \describe{
#'   \item{transcript}{Transcript ID of protein-coding RNAs}
#'   \item{description}{Name of domain families}
#'   \item{eval}{E-value score}
#'   \item{begin}{Start position of domain in protein}
#'   \item{end}{End position of domain in protein}
#'   ...
#' }
#' @usage data(domains.out)
"domains.out"

#' Example output of predictDomains()
#'
#' Output dataframe from predictDomains() function. mRNAs from GENCODE
#' mouse annotation was predicted for putative domain families.
#'
#' @format A data.frame with 85780 rows and 5 columns:
#' \describe{
#'   \item{transcript}{Transcript ID of protein-coding RNAs}
#'   \item{description}{Name of domain families}
#'   \item{eval}{E-value score}
#'   \item{begin}{Start position of domain in protein}
#'   \item{end}{End position of domain in protein}
#'   ...
#' }
#' @usage data(domains.known)
"domains.known"
