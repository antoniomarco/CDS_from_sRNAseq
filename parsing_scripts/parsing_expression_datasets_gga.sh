#!/usr/bin/bash

# miRBase gff3 annotation hsa
wget --no-check-certificate https://www.mirbase.org/download/gga.gff3


# Chicken Genome 5.0
wget http://ftp.ensembl.org/pub/release-106/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa.gz
# Index for HISAT2
hisat2-build Gallus_gallus.GRCg6a.dna_sm.toplevel.fa.gz gga5

# ENSEMBLE Chicken gene annotation  106 (for 5.0)
wget http://ftp.ensembl.org/pub/release-106/gtf/gallus_gallus/Gallus_gallus.GRCg6a.106.chr.gtf.gz
# Extract gene names for protein coding ones
zcat Gallus_gallus.GRCg6a.106.chr.gtf.gz | grep -e'biotype "protein_coding"' | awk -F'gene_id "' '{print $2}' | awk -F'"' '{print $1}' | uniq | sort | uniq > Gga_106.protein_coding.txt



# TODO UPDATE FILES FR CHICKEN!!!

# Brawand
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30352
# Fit/paired samples, they are all single-end reads
# 
mkdir GSE30352_Brawand_gga
cd GSE30352_Brawand_gga

# br M 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306711/SRR306711.fastq.gz
# lv M 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306720/SRR306720.fastq.gz
# ts M 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306722/SRR306722.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306721/SRR306721.fastq.gz
# ts M 2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306723/SRR306723.fastq.gz
# cb M 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306713/SRR306713.fastq.gz
# ht M 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306715/SRR306715.fastq.gz
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
  hisat2 -p 4 -x ../gga5 -U $file -S ${name}.sam 2> ${name}_hisat2.log
  echo "Compresing $file"
  samtools view -@ 4 -S -b ${name}.sam > ${name}.bam
#  rm ${name}.sam
done
#
# Gene and microRNA mapping
#
featureCounts -a ../Gallus_gallus.GRCg6a.106.chr.gtf.gz -t exon -g gene_id -o Brawand_mRNA_gga.tab *.bam



cd ..
# Meunier
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40499
mkdir GSE40499_Meunier_gga
cd GSE40499_Meunier_gga
#
# DOWNLOAD
#
# Brain
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553597/SRR553597.fastq.gz
# Kidney
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553600/SRR553600.fastq.gz
# Testis
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553601/SRR553601.fastq.gz
# Cerebellum
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553598/SRR553598.fastq.gz
# Heart
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553599/SRR553599.fastq.gz
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
  hisat2 -p 4 -x ../gga5 -U $file -S ${name}.sam 2> ${name}_hisat2.log
  echo "Compresing $file"
  samtools view -@ 4 -S -b ${name}.sam > ${name}.bam
#  rm ${name}.sam
done
#
# Gene and microRNA mapping
#
featureCounts -a ../Gallus_gallus.GRCg6a.106.chr.gtf.gz -t exon -g gene_id -o Meunier_mRNA_gga.tab *.bam
# microRNA to MATURE sequence! (not like the default to precursors)
featureCounts -t miRNA -g Name -a ../gga.gff3 -o Meunier_microRNA_gga.tab *.bam


exit 0
