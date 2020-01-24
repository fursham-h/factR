# factR

factR is a robust and easy-to-use R package with tools to process custom-assembled transcripts (GTF)

## Features
* Core features 
  1. Constructs coding (CDS) information for each transcript entry using a reference annotation as guide
  2. Analyzes coding transcripts for protein-coding features
  3. Scans coding transcripts for NMD-inducing features
  4. Reanalyzes coding transcripts for upstream and overlapping ORFs
  5. Compare and classify alternative segments between transcripts
* Supporting features 
  1. Matches chromosome levels of query GTF/object to a reference annotation
  2. Matches gene_id and gene_names of query GTF to a reference annotation
  3. Filter or sort each element of a GRangesList object

## Installation
```r
# install.packages("devtools")
devtools::install_github("fursham-h/factR")
```

## What you need
1. Assembled transcripts (GTF/GFF3) imported as GenomicRanges object
2. Reference annotation as GenomicRanges object. Obtained from:
    * Resource database including AnnotationHub
    * Import of reference annotation assembly (GTF/GFF3)
3. Genomic sequence. Obtained from:
    * Resource database including BSGenome, AnnotationHub
    * Import of genomic fasta file

## Getting started
See [vignette](https://htmlpreview.github.io/?https://github.com/fursham-h/factR/blob/master/vignettes/factR.html) for full instructions on how to get started

## Acknowledgements
* [Kaur Alasoo](https://github.com/kauralasoo)
