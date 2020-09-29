library("devtools")
library("GenomicFeatures")
library("wiggleplotr")
library("AnnotationHub")
library("rtracklayer")
library("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Hsapiens.UCSC.hg38")
load_all("../NMDer/")

# import mouse and human gencode basic annotations as GRanges objects
GRCh38_basic <- rtracklayer::import(
  "/Users/Fursham/Documents/Programming/Project/020_NMD/From_HD/gencode.v28.basic.annotation.gtf",
  format = "gtf"
)
GRCm38_basic <- rtracklayer::import(
  "/Users/Fursham/Documents/Programming/Project/020_NMD/From_HD/gencode.vM17.basic.annotation.gtf",
  format = "gtf"
)
mm10 <- GRCm38_basic
hg38 <- GRCh38_basic
devtools::use_data(mm10, overwrite = TRUE)
devtools::use_data(hg38, overwrite = TRUE)


# download mouse and human genome sequence from Ensembl
AnnHub <- AnnotationHub()
# query(AnnHub, c("Mus musculus", "release-91"))
Ensembl_mm10 <- AnnHub[["AH188"]] # Fasta for GRCm38 primary assembly
## is there a way to save the twobit file??
devtools::use_data(Ensembl_mm10, overwrite = TRUE)


# Make TxDb from ensembl
txdb_mm <- makeTxDbFromEnsembl("Mus musculus", server = "ensembldb.ensembl.org")
saveDb(txdb_mm, "data-raw/txdb/mus_musculus_txdb.sqlite")

# Load TxDb and exptract exons and cdss
txdb_mm <- loadDb("data-raw/txdb/mus_musculus_txdb.sqlite")
exons <- exonsBy(txdb_mm, by = "tx", use.names = TRUE)
cds <- cdsBy(txdb_mm, by = "tx", use.names = TRUE)


newStyle <- mapSeqlevels(seqlevels(exons), seqlevelsStyle(Mmusculus))
newStyle <- newStyle[!is.na(newStyle)]
exons <- renameSeqlevels(exons, newStyle)
cds <- renameSeqlevels(cds, newStyle)

seqlevels(exons, pruning.mode = "tidy") <- as.vector(newStyle)
seqlevels(cds, pruning.mode = "tidy") <- as.vector(newStyle)


# Extract Ptbp2 transcripts
ptbp2_testTx <- c("ENSMUST00000029780", "ENSMUST00000197833")
# ptbp2_NMDTx = cds$ENSMUST00000029780[cds$ENSMUST00000029780 != cds$ENSMUST00000029780[10]]
ptbp2Data <- list(
  transcripts = exons[ptbp2_testTx],
  cds = cds[ptbp2_testTx],
  refCDS = cds$ENSMUST00000029780,
  skipE10CDS = ptbp2_NMDTx,
  afCDS = exons$ENSMUST00000198399
)
devtools::use_data(ptbp2Data, overwrite = TRUE)

# Extract Bak1 transcripts
bak1_testTx <- c("ENSMUST00000078691", "ENSMUST00000025034")
# bak1_NMDTx = sort(unlist(append(
#  reduce(cds$ENSMUST00000078691),
#  reduce(exons$ENSMUST00000025034[5]))),
#  decreasing=TRUE)

bak1Data <- list(
  transcripts = exons[bak1_testTx],
  cds = cds[bak1_testTx],
  refCDS = cds$ENSMUST00000078691,
  poisonCDS = bak1_NMDTx
)
devtools::use_data(bak1Data, overwrite = TRUE)

# Extract Psd95 transcripts
psd95_testTx <- c("ENSMUST00000108589", "ENSMUST00000123687")
# psd95_NMDTx = cds$ENSMUST00000108589[cds$ENSMUST00000108589 != cds$ENSMUST00000108589[20]]
psd95Data <- list(
  transcripts = exons[psd95_testTx],
  cds = cds[psd95_testTx],
  refCDS = cds$ENSMUST00000108589,
  poisonCDS = psd95_NMDTx
)
devtools::use_data(psd95Data, overwrite = TRUE)

# Add test gtf file into raw data
#   File found in /../data-raw/transcript.gtf
testData <- "../NMDer/data-raw/transcripts.gtf"
devtools::use_data(testData, overwrite = TRUE)
