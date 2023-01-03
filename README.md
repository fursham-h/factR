# **factR v.1**

## Functional Annotation of Custom Transcriptomes in R

<!-- badges: start -->
[![R build status](https://github.com/fursham-h/factR/workflows/R-CMD-check/badge.svg)](https://github.com/fursham-h/factR/actions)
[![Codecov test coverage](https://github.com/fursham-h/factR/workflows/test-coverage/badge.svg)](https://github.com/fursham-h/factR/actions)
[![Codecov test coverage](https://codecov.io/gh/fursham-h/factR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/fursham-h/factR?branch=master)
<!-- badges: end -->
  
## General workflow
<p align="center">
  <img src="man/figures/factR_workflow.png" width="450"/>
</p>

*factR* is a robust and easy-to-use R package with tools to process custom-assembled transcriptomes (GTF). Below are *factR*'s key functions:

* Core features 
  1. Construct transcript coding (CDS) information 
  using a reference-guided process
  2. Predict protein domains on coding transcripts
  3. Predict sensitivity of coding transcripts to Nonsense-mediated decay
* Supporting features 
  1. Match chromosome levels of query GTF/object to reference annotation
  2. Match gene_id and gene_names of query GTF to reference annotation
  3. Plot transcripts from GTF GRanges object using *wiggleplotr*
  4. Subset new transcripts from custom transcriptome

## How to install
The latest stable version can be installed directly from [Bioconductor](https://bioconductor.org/packages/release/bioc/html/factR.html):
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("factR")
```

Alternatively, you may install the development version of *factR* using devtools:
```r
# install.packages("devtools")
devtools::install_github("fursham-h/factR")
```

## Getting started
See our [quickstart guide](https://fursham-h.github.io/factR/articles/quickstart.html) or our 
[full vignette](https://fursham-h.github.io/factR/articles/factR.html) on how to get started

## Acknowledgements
We thank [Kaur Alasoo](https://github.com/kauralasoo) for sharing code 
resources for *wiggleplotr* and for valuable discussions on the design 
of the package.

## Citing factR
Please cite the following references if you use factR:

1. Fursham Hamid, Kaur Alasoo, Jaak Vilo, Eugene Makeyev (2022); Functional annotation of custom transcriptomes; Methods in Molecular Biology
2. [Fursham Hamid (2022); Functional Annotation of Custom Transcriptomes; Bioconductor](https://bioconductor.org/packages/release/bioc/html/factR.html)






