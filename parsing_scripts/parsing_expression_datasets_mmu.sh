#!/usr/bin/bash

# miRBase gff3 annotation hsa
wget --no-check-certificate https://www.mirbase.org/download/mmu.gff3
# Mouse genecode, only CHR, basic gene annottion. M25 (last in mm10)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz
# Extract gene names for protein coding ones
zcat gencode.vM25.basic.annotation.gtf.gz | grep -e'gene_type "protein_coding"' | awk -F'gene_name "' '{print $2}' | awk -F'"' '{print $1}' | uniq | sort | uniq > gencode.vM25.protein_coding.txt

# Brawand
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30352
# Fit/paired samples, they are all single-end reads
# 
mkdir GSE30352_Brawand_mmu
cd GSE30352_Brawand_mmu
# Heart Male 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306767/SRR306767.fastq.gz
# Testis Male 2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306776/SRR306776.fastq.gz
# Cerebellum Male 1. From Male C57BL
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306764/SRR306764.fastq.gz
# Brain Male 1. From Male C57BL
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306758/SRR306758.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306759/SRR306759.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306761/SRR306761.fastq.gz
# Kidney Female 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306769/SRR306769.fastq.gz

# #
# # Quality control and adaptor trimming
# #
# for file in *.fastq.gz;
# do
#   fastqc $file
# done
#
# Mapping
#
for file in *.fastq.gz;
do
  echo "Mapping $file"
  name=${file%".fastq.gz"}
  hisat2 -p 4 -x mm10 -U $file -S ${name}.sam 2> ${name}_hisat2.log
  echo "Compresing $file"
  samtools view -@ 4 -S -b ${name}.sam > ${name}.bam
#  rm ${name}.sam
done
#
# Gene and microRNA mapping
#
featureCounts -a ../gencode.vM25.basic.annotation.gtf.gz -t exon -g gene_name -o Brawand_mRNA_mmu.tab *.bam



cd ..
# Meunier
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40499
mkdir GSE40499_Meunier_mmu
cd GSE40499_Meunier_mmu
#
# DOWNLOAD
#
# Testis
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553586/SRR553586.fastq.gz
# Brain. From Male C57BL
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553582/SRR553582.fastq.gz
# Cerebellum
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553583/SRR553583.fastq.gz
# Heart
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553584/SRR553584.fastq.gz
# Kidney
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553585/SRR553585.fastq.gz
# #
# # Quality control and adaptor trimming
# #
# for file in *.fastq.gz;
# do
#   fastqc $file
# done
for file in *.fastq.gz;
do
  name=${file%".fastq.gz"}
  cutadapt -a ATCTCGTATGCC -m 18 -o ${name}_t.fastq.gz $file > ${name}_cutadapt.log
done
#
# Mapping
#
for file in *_t.fastq.gz;
do
  echo "Mapping $file"
  name=${file%"_t.fastq.gz"}
  hisat2 -p 4 -x mm10 -U $file -S ${name}.sam 2> ${name}_hisat2.log
  echo "Compresing $file"
  samtools view -@ 4 -S -b ${name}.sam > ${name}.bam
#  rm ${name}.sam
done
#
# Gene and microRNA mapping
#
featureCounts -a ../gencode.vM25.basic.annotation.gtf.gz -t exon -g gene_name -o Meunier_mRNA_mmu.tab *.bam
# microRNA to MATURE sequence! (not like the default to precursors)
featureCounts -t miRNA -g Name -a ../mmu.gff3 -o Meunier_microRNA_mmu.tab *.bam


exit 0
