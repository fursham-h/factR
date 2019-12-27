# factR

factR is a robust and easy-to-use R package for annotating custom-assembled transcripts (GTF/GFF3) guided by a reference annotation. 

## Features
* Core features 
* Constructs coding (CDS) information for each transcript using a reference annotation as guide
* Searches coding transcripts for protein features
* Scans coding transcripts for NMD-inducing features
* Reanalyzes coding transcripts for upstream and overlapping ORFs
* Secondary features 
* Matches chromosome levels and gene_id levels between custom gtf and reference annotation 
* Calculates percent base coverage between two transcripts or two list-of-transcripts
* Resizes start and end of transcript-organized GRanges object
* *In development features*
  * *Extract more features to be predicted from translated proteins*
    * *Transmembrane domain prediction*
    * *Disordered domain prediction*
  * *Compare alternative spliced segments between isoforms*
  * *Compare biological significance between mRNA isoforms*

## Installation
```r
# install.packages("devtools")
devtools::install_github("fursham-h/factR"
```

## What you need

## Getting started
See vignette for full instructions on how to 

## Acknowledgements
* [Kaur Alasoo](https://github.com/kauralasoo)
