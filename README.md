# **factR**: Functional Annotation of Custom Transcriptomes in R

<!-- badges: start -->
  [![R build status](https://github.com/fursham-h/factR/workflows/R-CMD-check/badge.svg$branch=master)](https://github.com/fursham-h/factR/actions)
  [![Codecov test coverage](https://github.com/fursham-h/factR/workflows/test-coverage/badge.svg$branch=dev)](https://github.com/fursham-h/factR/actions)
<!-- badges: end -->
  
## General workflow
<p align="center">
  <img src="factR_workflow.png" width="450"/>
</p>

*factR* is a robust and easy-to-use R package with tools to process custom-assembled transcriptomes (GTF). Below are *factR*'s key functions:

* Core features 
  1. Construct transcript coding (CDS) information using a reference-guided process
  2. Predict protein domains on coding transcripts
  3. Predict sensitivity of coding transcripts to Nonsense-mediated decay
* Supporting features 
  1. Match chromosome levels of query GTF/object to reference annotation
  2. Match gene_id and gene_names of query GTF to reference annotation
  3. Plot transcripts from GTF GRanges object using *wiggleplotr*
  4. Select new transcripts from custom transcriptome

## How to install
```r
# install.packages("devtools")
devtools::install_github("fursham-h/factR")
```

## What you need
1. Custom-assembled transcriptome (GTF)
2. Reference annotation as GenomicRanges object. Obtained from:
    * Resource database including AnnotationHub
    * Import of reference annotation assembly (GTF/GFF3)
3. Genomic sequence. Obtained from:
    * Resource database including BSGenome, AnnotationHub
    * Import of genomic fasta file


## Getting started
See our [vignette](https://fursham-h.github.io/factR/articles/factR.html) for full instructions on how to get started

## Reporting bugs and issues
If you encounter a bug with *factR*, please file an issue on [GitHub](https://github.com/fursham-h/factR/issues).

## Acknowledgements
We thank [Kaur Alasoo](https://github.com/kauralasoo) for sharing code resources for *wiggleplotr* and for valuable discussions on the design of the package.
