## ----setup, echo=FALSE, cache=FALSE-------------------------------------------
## Global options
options(max.print="75")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")
knitr::opts_knit$set(width=75)

## ----silentload, include=FALSE------------------------------------------------
#Silent load all dependencies for vignette
library(factR)
library(AnnotationHub)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicFeatures)
library(rtracklayer)

## ----install.factR, eval = FALSE----------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("fursham-h/factR")

## ----install.others, eval = FALSE---------------------------------------------
#  BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

## ----load.dep.workflow, eval = FALSE------------------------------------------
#  library(factR)
#  library(AnnotationHub)
#  library(BSgenome.Mmusculus.UCSC.mm10)
#  library(Biostrings)
#  library(GenomicFeatures)
#  library(rtracklayer)
#  library(purrr)
#  library(stringr)

## ----load,query, eval = FALSE-------------------------------------------------
#  # import gtf file using rtracklayer
#  query_gtf <- import('path/to/your.gtf')

## ----load.gencode, eval = FALSE-----------------------------------------------
#  # query database for mouse gencode basic annotation
#  ah <- AnnotationHub()
#  query(ah, c('Mus musculus', 'gencode', 'gff', 'basic'))
#  #> AnnotationHub with 2 records
#  #> # snapshotDate(): 2019-05-02
#  #> # $dataprovider: Gencode
#  #> # $species: Mus musculus
#  #> # $rdataclass: GRanges
#  #> # additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer,
#  #> #   rdatadateadded, preparerclass, tags, rdatapath, sourceurl, sourcetype
#  #> # retrieve records with, e.g., 'object[["AH49547"]]'
#  #>
#  #>             title
#  #>   AH49547 | gencode.vM6.basic.annotation.gff3.gz
#  #>   AH49549 | gencode.vM6.chr_patch_hapl_scaff.basic.annotation.gff3.gz
#  
#  #>  retrieve id 'AH49547'
#  ref_gtf <- ah[['AH49547']]

## ----load.ensembl, eval = FALSE-----------------------------------------------
#  # Load annotation hub and query database for ensembl annotation as in method 1 above.
#  
#  # retrieve latest GRCm38 annotation from ensembl
#  ref_gtf <- ah[['AH60127']]

## ----load.localref, eval = FALSE----------------------------------------------
#  # import reference gtf file
#  ref_gtf <- import('path/to/reference.gtf')

## ----load.BSgenome, eval = FALSE----------------------------------------------
#  Mmusculus <- BSgenome.Mmusculus.UCSC.mm10

## ----load.fastaAH, eval = FALSE-----------------------------------------------
#  # query database for mouse gencode basic annotation
#  ah <- AnnotationHub()
#  query(ah, c("Mus musculus", "release-91", 'fasta', 'GRCm38'))
#  #> AnnotationHub with 5 records
#  #> # snapshotDate(): 2019-05-02
#  #> # $dataprovider: Ensembl
#  #> # $species: Mus musculus
#  #> # $rdataclass: TwoBitFile
#  #> # additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer,
#  #> #   rdatadateadded, preparerclass, tags, rdatapath, sourceurl, sourcetype
#  #> # retrieve records with, e.g., 'object[["AH60491"]]'
#  #>
#  #>             title
#  #>   AH60491 | Mus_musculus.GRCm38.cdna.all.2bit
#  #>   AH60492 | Mus_musculus.GRCm38.dna.primary_assembly.2bit
#  #>   AH60493 | Mus_musculus.GRCm38.dna_rm.primary_assembly.2bit
#  #>   AH60494 | Mus_musculus.GRCm38.dna_sm.primary_assembly.2bit
#  #>   AH60495 | Mus_musculus.GRCm38.ncrna.2bit
#  
#  #>  retrieve id 'AH60492'
#  Mmusculus <- ah[['AH60492']]

## ----load.localfasta, eval = FALSE--------------------------------------------
#  # load package and import fasta file as a DNAStringSet
#  Mmusculus <- readDNAStringSet('path/to/fasta')

## ----convert.fastaheader, eval = FALSE----------------------------------------
#  # load package and import fasta file as a DNAStringSet
#  names(Mmusculus) <- str_split(seqlevels(Mmusculus), " ") %>% map_chr(`[`, 1)

## ----checkquery---------------------------------------------------------------
# inspect query object
head(query_gtf)

# list out transcript_id and gene_id names
unique(query_gtf$transcript_id)

unique(query_gtf$gene_id)

## ----viewquery----------------------------------------------------------------
viewTranscripts(query_gtf)

## ----viewref------------------------------------------------------------------
# inspect reference object
head(ref_gtf)

# list out transcript_id and gene_id names
unique(ref_gtf$transcript_id)

unique(ref_gtf$gene_id)

## ----inspectobjects-----------------------------------------------------------
seqlevelsStyle(query_gtf)

seqlevelsStyle(ref_gtf)

## ----matchseqlvls-------------------------------------------------------------
# matching seqlevels
matched_query_gtf <- matchChromosomes(query_gtf, Mmusculus)

head(matched_query_gtf)

## ----testseqlevels------------------------------------------------------------
has_consistentSeqlevels(query_gtf, ref_gtf)

has_consistentSeqlevels(matched_query_gtf, ref_gtf)

## ----matchgenemeta------------------------------------------------------------
# matching gene metadata
matched_query_gtf <- matchGeneInfo(matched_query_gtf, ref_gtf)

head(matched_query_gtf)

## ----subsettranscripts--------------------------------------------------------
subsetted_query_gtf <- subsetNewTranscripts(matched_query_gtf, ref_gtf)

# inspect length of subsetted object
length(subsetted_query_gtf)

length(matched_query_gtf)

# list de novo transcripts
unique(subsetted_query_gtf$transcript_id)

## ----buildcds-----------------------------------------------------------------
new_query_gtf <- buildCDS(matched_query_gtf, ref_gtf, Mmusculus)

## ----checkCDS-----------------------------------------------------------------
# check lengths of input and output GRanges
length(matched_query_gtf)

length(new_query_gtf)

# inspect CDS entries
head(new_query_gtf[new_query_gtf$type == "CDS"])

## ----viewCDS------------------------------------------------------------------
viewTranscripts(new_query_gtf)

## ----export, eval=FALSE-------------------------------------------------------
#  # get GrangesList of query exons and cds
#  rtracklayer::export(new_query_gtf, 'path/to/output.gtf', format = 'gtf')

## ----testfeature--------------------------------------------------------------
domains <- predictDomains(new_query_gtf, Mmusculus)

## ----testfeaturewhich, eval=F-------------------------------------------------
#  pd <- predictDomains(new_query_gtf, Mmusculus, transcript_id %in% c('transcript1', 'transcript3'))
#  pd
#  #>    transcript   description     eval begin end
#  #> 1 transcript1 Canonical RBD 4.96e-06    46 143
#  #> 2 transcript1 Canonical RBD 1.39e-05   357 446
#  #> 3 transcript1 Canonical RBD 1.28e-06   177 281
#  #> 4 transcript1 Canonical RBD 9.32e-06   469 553
#  #> 5 transcript3 Canonical RBD 4.96e-06    46 143
#  #> 6 transcript3 Canonical RBD 1.28e-06   177 281

## ----subsetandtest, eval=F----------------------------------------------------
#  subsetted_new_query_gtf <- subsetNewTranscripts(new_query_gtf, ref_gtf)
#  pd_subset <- predictDomains(subsetted_new_query_gtf, Mmusculus)
#  pd_subset
#  
#  #>    transcript   description     eval begin end
#  #> 1 transcript3 Canonical RBD 4.96e-06    46 143
#  #> 2 transcript3 Canonical RBD 1.28e-06   177 281
#  #> 3 transcript4 Canonical RBD 1.39e-05   291 380
#  #> 4 transcript4 Canonical RBD 1.28e-06   137 241
#  #> 5 transcript4 Canonical RBD 8.52e-05    39  95
#  #> 6 transcript4 Canonical RBD 8.32e-06   403 487

## ----testnmd1-----------------------------------------------------------------
pn <- predictNMD(new_query_gtf)

## ----testnmdfilter------------------------------------------------------------
pn <- predictNMD(new_query_gtf, transcript_id == 'transcript3')
pn

## ----testnmdens, eval = FALSE-------------------------------------------------
#  # retrieve latest GRCm38 annotation from ensembl
#  ah <- AnnotationHub()
#  mm10_gtf <- ah[['AH60127']]
#  
#  # predict transcripts from Ptbp1 gene for NMD features
#  pn_Ptbp1 <- predictNMD(mm10_gtf, gene_name == 'Ptbp1')
#  pn_Ptbp1
#  #> # A tibble: 9 x 6
#  #>   transcript         stop_to_lastEJ num_of_downEJs stop_to_downEJs threeUTRlength is_NMD
#  #>   <chr>                       <dbl>          <int> <chr>                    <dbl> <lgl>
#  #> 1 ENSMUST00000057343            364              3 "69,286,364"               644 TRUE
#  #> 2 ENSMUST00000095457           -130              0 ""                        1085 FALSE
#  #> 3 ENSMUST00000165704           -130              0 ""                        1497 FALSE
#  #> 4 ENSMUST00000165724            286              2 "69,286"                   362 TRUE
#  #> 5 ENSMUST00000168683           -145              0 ""                           0 FALSE
#  #> 6 ENSMUST00000169091            -36              0 ""                           0 FALSE
#  #> 7 ENSMUST00000169483            139              1 "139"                      152 TRUE
#  #> 8 ENSMUST00000171599            -63              0 ""                           0 FALSE
#  #> 9 ENSMUST00000172282           -130              0 ""                        1158 FALSE
#  #> Warning message:
#  #> 7 transcript(s) have missing cds info and was not analyzed
#  
#  # predict several transcripts from Ptbp1 gene for NMD features
#  pn_Ptbp1_transcripts <- predictNMD(mm10_gtf, gene_name == 'Ptbp1', transcript_id == c('ENSMUST00000057343', 'ENSMUST00000095457'))
#  pn_Ptbp1_transcripts
#  #> # A tibble: 2 x 6
#  #>   transcript         stop_to_lastEJ num_of_downEJs stop_to_downEJs threeUTRlength is_NMD
#  #>   <chr>                       <dbl>          <int> <chr>                    <dbl> <lgl>
#  #> 1 ENSMUST00000057343            364              3 "69,286,364"               644 TRUE
#  #> 2 ENSMUST00000095457           -130              0 ""                        1085 FALSE

## ----sessioninfo--------------------------------------------------------------
sessionInfo()

