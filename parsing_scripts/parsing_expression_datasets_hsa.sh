#!/usr/bin/bash

# miRBase gff3 annotation hsa
wget --no-check-certificate https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3
# Human genecode, only CHR, basic gene annottion
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz
# Extract gene names for protein coding ones
zcat gencode.v43.basic.annotation.gtf.gz | grep -e'gene_type "protein_coding"' | awk -F'gene_name "' '{print $2}' | awk -F'"' '{print $1}' | uniq | sort | uniq > gencode.v43.protein_coding.txt

# Brawand
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30352
# Fit/paired samples, they are all single-end reads
# 
mkdir GSE30352_Brawand
cd GSE30352_Brawand
# Br Male 3
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306839/SRR306839.fastq.gz
# Cb Male 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306845/SRR306845.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306846/SRR306846.fastq.gz
# Heart Male 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306848/SRR306848.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306849/SRR306849.fastq.gz
# Kidney Male 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR306/SRR306852/SRR306852.fastq.gz
#
# Quality control and adaptor trimming
#
for file in *.fastq.gz;
do
  fastqc $file
done
#
# Mapping
#
for file in *.fastq.gz;
do
  echo "Mapping $file"
  name=${file%".fastq.gz"}
  hisat2 -p 4 -x hg38 -U $file -S ${name}.sam 2> ${name}_hisat2.log
  echo "Compresing $file"
  samtools view -@ 4 -S -b ${name}.sam > ${name}.bam
#  rm ${name}.sam
done
#
# Gene and microRNA mapping
#
featureCounts -a ../gencode.v43.basic.annotation.gtf.gz -t exon -g gene_name -o Brawand_mRNA.tab *.bam



cd ..
# Meunier
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40499
mkdir GSE40499_Meunier
cd GSE40499_Meunier
#
# DOWNLOAD
#
# Frontal Cortex
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553572/SRR553572.fastq.gz
# Cerebellum
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553573/SRR553573.fastq.gz
# Heart
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553574/SRR553574.fastq.gz
# Kidney
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553575/SRR553575.fastq.gz
# Testis
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR553/SRR553576/SRR553576.fastq.gz
#
# Quality control and adaptor trimming
#
for file in *.fastq.gz;
do
  fastqc $file
done
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
  hisat2 -p 4 -x hg38 -U $file -S ${name}.sam 2> ${name}_hisat2.log
  echo "Compresing $file"
  samtools view -@ 4 -S -b ${name}.sam > ${name}.bam
#  rm ${name}.sam
done
#
# Gene and microRNA mapping
#
featureCounts -a ../gencode.v43.basic.annotation.gtf.gz -t exon -g gene_name -o Meunier_mRNA.tab *.bam
# microRNA to MATURE sequence! (not like the default to precursors)
featureCounts -t miRNA -g Name -a ../hsa.gff3 -o Meunier_microRNA.tab *.bam


## Mapping expression reads to CDS and UTRs
# Get UTR/CDS coordinates
Rscript generate_utr_cds_gff.R

cd GSE30352_Brawand
featureCounts -a ../UTR_CDS_coordinates.gtf.gz -t CDS -g gene_id -o Brawand_mRNA_CDS.tab *.bam
featureCounts -a ../UTR_CDS_coordinates.gtf.gz -t UTR5 -g gene_id -o Brawand_mRNA_UTR5.tab *.bam
featureCounts -a ../UTR_CDS_coordinates.gtf.gz -t UTR3 -g gene_id -o Brawand_mRNA_UTR3.tab *.bam

cd ../GSE40499_Meunier
featureCounts -a ../UTR_CDS_coordinates.gtf.gz -t CDS -g gene_id -o Meunier_mRNA_CDS.tab *.bam
featureCounts -a ../UTR_CDS_coordinates.gtf.gz -t UTR5 -g gene_id -o Meunier_mRNA_UTR5.tab *.bam
featureCounts -a ../UTR_CDS_coordinates.gtf.gz -t UTR3 -g gene_id -o Meunier_mRNA_UTR3.tab *.bam

cd ../
tar -cf hsa_mapping_UTR_CDS.tar GSE30352_Brawand/Brawand_mRNA_CDS.tab GSE30352_Brawand/Brawand_mRNA_UTR5.tab GSE30352_Brawand/Brawand_mRNA_UTR3.tab GSE40499_Meunier/Meunier_mRNA_CDS.tab GSE40499_Meunier/Meunier_mRNA_UTR5.tab GSE40499_Meunier/Meunier_mRNA_UTR3.tab
gzip hsa_mapping_UTR_CDS.tar

exit 0
