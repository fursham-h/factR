#usr/bin/bash

# Section 3.2.3: Retrieve data
mkdir FASTQ_sc
while IFS=$'\t' read -r NAME CLASS LINK CLASSNUM;do
 	wget $LINK -P FASTQ_sc
 	tar -xvf "FASTQ_sc/$NAME.fastq.tar" -C FASTQ_sc
 	rm FASTQ_sc/$NAME.fastq.tar
done < sc_fastq_curated.txt

# Section 3.2.4: Merge FASTQ by apparent cell type
# Create arrays to store filenames and group names
declare -A fastqArray
declare -A fastqNames

while IFS=$'\t' read -r NAME CLASS LINK;do
    fastqArray["$CLASS.1"]="${fastqArray["$CLASS.1"]} FASTQ_sc/${NAME}_R1.fastq.gz"
    fastqArray["$CLASS.2"]="${fastqArray["$CLASS.2"]} FASTQ_sc/${NAME}_R2.fastq.gz"
    fastqNames["$CLASS.1"]=${CLASS}_R1
    fastqNames["$CLASS.2"]=${CLASS}_R2
done < sc_fastq_curated.txt

# Loop groups and concatenate FASTQ files
for n in Glutamatergic GABAergic Endothelial Astrocyte;do
	cat ${fastqArray["$n.1"]} > FASTQ_sc/${fastqNames["$n.1"]}.fastq.gz
	cat ${fastqArray["$n.2"]} > FASTQ_sc/${fastqNames["$n.2"]}.fastq.gz
done

# Section 3.2.5: Download index
if [ ! -d "mm10" ]; then
	wget https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
	tar -xf "mm10_genome.tar.gz"
	rm mm10_genome.tar.gz
fi

# Section 3.2.6: Align reads
cores=`nproc`
mkdir Hisat2_SAMs_sc
for fastq in Glutamatergic GABAergic Endothelial Astrocyte;do
	hisat2 -p $cores --dta -x mm10/genome -1 FASTQ_sc/"$fastq"_R1.fastq.gz -2 FASTQ_sc/"$fastq"_R2.fastq.gz -S Hisat2_SAMs_sc/$fastq.sam
done

# Section 3.2.7: Convert SAM to BAM
mkdir Hisat2_sorted_BAMs_sc
for file in Glutamatergic GABAergic Endothelial Astrocyte; do
	samtools view -@ $cores -Su Hisat2_SAMs_sc/$file.sam | samtools sort -@ $cores -o Hisat2_sorted_BAMs_sc/$file.bam
done

# Section 3.2.8: Download reference annotation
if [ ! -f "gencode.vM25.annotation.gtf" ]; then
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
	gzip -d gencode.vM25.annotation.gtf.gz
fi

# Section 3.2.9: Assemble transcriptome
mkdir Stringtie_gtf_sc
for file in Glutamatergic GABAergic Endothelial Astrocyte; do
	stringtie Hisat2_sorted_BAMs_sc/$file.bam -p $cores -o Stringtie_gtf_sc/$file.gtf -G gencode.vM25.annotation.gtf
done

# Section 3.2.10: Merge transcriptome and compress
stringtie --merge -G gencode.vM25.annotation.gtf Stringtie_gtf_sc/Glutamatergic.gtf Stringtie_gtf_sc/GABAergic.gtf Stringtie_gtf_sc/Endothelial.gtf Stringtie_gtf_sc/Astrocyte.gtf | gzip > sc_merged.gtf.gz 

