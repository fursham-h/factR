# **factR**: Functional Annotation of Custom Transcriptomes in R

*factR* is a robust and easy-to-use R package with tools to process custom-assembled transcriptomes (GTF). Below are *factR*'s key functions:

* Core features 
  1. Construct transcript coding (CDS) information using a reference-guided process
  2. Predict protein domains on coding transcripts
  3. Predict sensitivity of coding transcripts to Nonsense-mediated decay
  4. Compare and classify alternative segments between transcripts
* Supporting features 
  1. Match chromosome levels of query GTF/object to reference annotation
  2. Match gene_id and gene_names of query GTF to reference annotation
  3. Plot transcripts from GTF GRanges object using *wiggleplotr*
  4. Select new transcripts from custom transcriptome

### How to install
```r
# install.packages("devtools")
devtools::install_github("fursham-h/factR")
```

### What you need
1. Custom-assembled transcriptome (GTF)
2. Reference annotation as GenomicRanges object. Obtained from:
    * Resource database including AnnotationHub
    * Import of reference annotation assembly (GTF/GFF3)
3. Genomic sequence. Obtained from:
    * Resource database including BSGenome, AnnotationHub
    * Import of genomic fasta file


### Getting started
See [vignette](https://htmlpreview.github.io/?https://github.com/fursham-h/factR/blob/dev/doc/factR.html) for full instructions on how to get started


### Acknowledgements
We thank [Kaur Alasoo](https://github.com/kauralasoo) for sharing code resources for *wiggleplotr* and for valuable discussions on the design of the package.
