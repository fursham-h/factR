# Generating extdata
The extdata `sc_merged_sample.gtf.gz` contain a sample custom transcriptome
constructed from cortical neural populations from a dataset by Tasic et al. (2018).
To generate the full transcriptome, execute the download_assemble_sc.sh script:

```bash
bash download_assemble_sc.sh
```

A sample of this transcriptome was obtained by taking the first 500 genes
from chromosome 15. This is done in R:
```r
gtf <- factR::importGTF("sc_merged.gtf.gz")
chr15.genes <- unique(gtf[seqnames(gtf) == "chr15"]$gene_id)
chr15.500 <- chr15.genes[1:500]
gtf <- gtf[gtf$gene_id %in% chr15.500]
rtracklayer::export(gtf, "sc_merged_sample.gtf.gz")
```