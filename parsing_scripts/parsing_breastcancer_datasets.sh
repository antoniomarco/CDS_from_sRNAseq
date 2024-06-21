#!/usr/bin/bash


# miRBase gff3 annotation hsa
wget --no-check-certificate https://www.mirbase.org/download/hsa.gff3
# Human genecode, only CHR, basic gene annottion
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz

# Data from 
# https://www.ebi.ac.uk/ena/browser/view/PRJNA494326
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/003/SRR7948463/SRR7948463.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/005/SRR7948465/SRR7948465.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/009/SRR7948469/SRR7948469.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/000/SRR7948470/SRR7948470.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/003/SRR7948473/SRR7948473.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/004/SRR7948474/SRR7948474.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/004/SRR7948484/SRR7948484.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/005/SRR7948485/SRR7948485.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/008/SRR7948468/SRR7948468.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/001/SRR7948471/SRR7948471.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/002/SRR7948472/SRR7948472.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/005/SRR7948475/SRR7948475.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/008/SRR7948478/SRR7948478.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/009/SRR7948479/SRR7948479.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/002/SRR7948482/SRR7948482.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/002/SRR7948462/SRR7948462.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/004/SRR7948464/SRR7948464.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/006/SRR7948466/SRR7948466.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/007/SRR7948467/SRR7948467.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/006/SRR7948476/SRR7948476.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/007/SRR7948477/SRR7948477.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/000/SRR7948480/SRR7948480.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/001/SRR7948481/SRR7948481.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR794/003/SRR7948483/SRR7948483.fastq.gz

#
# Quality control and adaptor trimming
#
for file in *.fastq.gz;
do
  fastqc $file
done

#
# Trim adaptor (universal)
#
for file in *.fastq.gz;
do
  name=${file%".fastq.gz"}
  cutadapt -a AGATCGGAAGAG -m 18 -o ${name}_t.fastq.gz $file > ${name}_cutadapt.log
done

#
# Mapping
#
for file in *_t.fastq.gz;
do
  echo "Mapping $file"
  name=${file%"_t.fastq.gz"}
  hisat2 -p 4 -x hg38 -U $file -S ${name}.sam 2> ${name}_hisat2.log
  echo "Compresing $file"
  samtools view -@ 4 -S -b ${name}.sam > ${name}.bam
  rm ${name}.sam
done

#
# Gene and microRNA mapping
#
featureCounts -a gencode.v43.basic.annotation.gtf.gz -t exon -g gene_name -o bc_patients_mRNA.tab *.bam
# microRNA to MATURE sequence! (not like the default to precursors)
featureCounts -t miRNA -g Name -a hsa.gff3 -o bc_patients_microRNA.tab *.bam


exit 0
