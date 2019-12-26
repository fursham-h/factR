# factR

## Description

New Description

## Features
* Scans for NMD-inducing features on coding RNAs from various sources:
  * Most transcript annotation databases
  * Custom transcript annotation (GTF) generated from RNA-seq experiemnts
* Constructs coding region (CDS) on custom-generated GTF annotation using reference CDS as guide
* Extract features of coding proteins. Currently support:
  * Protein domains prediction
  * Signal peptide prediction
* Matches chromosome levels and gene_id levels between query and reference GRanges object
* Search coding mRNAs for uORFs and overlapping ORFs
* Resizes start and end of transcript-organized GRanges object
* *In development*
  * *Extract more features to be predicted from translated proteins*
    * *Transmembrane domain prediction*
    * *Disordered domain prediction*
  * *Compare alternative spliced segments between isoforms*
  * *Compare biological significance between mRNA isoforms*
  
See vignette for full instructions on how to get started

## Installation
```r
# install.packages("devtools")
devtools::install_github("fursham-h/ponder")
```

## Example
pondeR is packaged with sample GenomicRangesList objects containing 
exon and cds ranges of 4 transcripts from the same gene
```r
library(pondeR)

names(query_exons[1])
# [1] "transcript1" "transcript2" "transcript3" "transcript4"
names(query_cds[1])
# [1] "transcript1" "transcript2" "transcript3" "transcript4"

```

Transcript architecture of these transcripts can be visualized 
using _wiggleplotr_ package
```r
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#BiocManager::install("wiggleplotr")

library(wiggleplotr)
plotTranscripts(query_exons, query_cds)
```
<img src="wiggleplot_query.png" width="600">

Run predictNMD with exon and cds information as input
```r
predictNMD(query_exons, query_cds)
## A tibble: 1 x 6
#  tx          is_NMD dist_to_lastEJ num_of_down_EJs dist_to_downEJs threeUTRlength
#  <chr>       <lgl>           <int>           <dbl> <chr>                    <dbl>
#1 transcript3 TRUE              361               3 66,283,361                 641
```

## Acknowledgements
* [Kaur Alasoo](https://github.com/kauralasoo)
