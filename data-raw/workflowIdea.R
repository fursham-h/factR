library(pondeR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(wiggleplotr)
library(GenomicFeatures)
library(dplyr)

checkInputs(query_gtf, ref_gtf, Mmusculus)

query_gtf <- matchSeqLevels(query_gtf, ref_gtf)
query_gtf <- matchGeneIDs(query_gtf, ref_gtf)

# make txdb and exonsBy
queryDB <- makeTxDbFromGRanges(query_gtf)
refDB <- makeTxDbFromGRanges(ref_gtf)
query_exons <- exonsBy(queryDB, by = "tx", use.names = TRUE)
ref_cds <- cdsBy(refDB, by = "tx", use.names = TRUE)
ref_exons <- exonsBy(refDB, by = "tx", use.names = TRUE)

# get coverage
query_ids <- query_gtf %>%
  as.data.frame() %>%
  dplyr::select(gene_id, transcript_id) %>%
  distinct()
ref_ids <- ref_gtf %>%
  as.data.frame() %>%
  dplyr::select(gene_id, ref_transcript_id = transcript_id) %>%
  distinct()
q2r <- left_join(query_ids, ref_ids) %>%
  dplyr::select(-gene_id)

q2rcovs <- getCoverages(query_exons, ref_exons, q2r)

# getCDS
query_cds <- buildCDS(query_exons, ref_cds,
  Mmusculus, q2rcovs,
  coverage = 3
)

# export exons and cds to gtf
makeGTF(query_exons, query_cds, ".")

# refine uORF and uATG

# testNMD
predictNMD(query_exons, query_cds)

# Features of translated proteins from CDS list
extractCDSfeature(query_cds, Mmusculus)

# get ensembl annotation
AnnHub <- AnnotationHub()
Mmus <- AnnHub[["AH60127"]]
