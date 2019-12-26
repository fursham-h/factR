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

## Installation
```r
# install.packages("devtools")
devtools::install_github("fursham-h/factR")
```

## Getting started
See vignette for full instructions on how to get started

## Acknowledgements
* [Kaur Alasoo](https://github.com/kauralasoo)
