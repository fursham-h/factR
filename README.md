# factR

factR is a robust and easy-to-use R package for post-processing of custom-assembled
transcripts (GTF/GFF3)

## Features
* Core features 
  1. Constructs coding (CDS) information for each transcript entry using a reference annotation as guide
  2. Analyzes coding transcripts for protein-coding features
  3. Scans coding transcripts for NMD-inducing features
  4. Reanalyzes coding transcripts for upstream and overlapping ORFs
* Supporting features 
  1. Matches chromosome levels of query GTF/object to a reference annotation
  2. Matches gene_id and gene_names of query GTF to a reference annotation
  3. Calculates percent base coverage between two transcripts or two list-of-transcripts
  4. Trims start and end of multi-exon GRanges transcript

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
See [vignette](https://htmlpreview.github.io/?https://github.com/fursham-h/factR/blob/master/vignettes/factR.html) for full instructions on how to make the most of your assembled transcripts

## Acknowledgements
* [Kaur Alasoo](https://github.com/kauralasoo)
