library(biomaRt) 
library(tidyverse)

# Connect to ENSEMBL
ensembl <-   biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Gen ENSEMBL info
attributes <- c(
  "external_gene_name",
  "ensembl_gene_id",
  "ensembl_transcript_id",
  "chromosome_name",
  "strand",
  "5_utr_start",
  "5_utr_end",
  "3_utr_start",
  "3_utr_end",
  "genomic_coding_start",
  "genomic_coding_end"
)
ensembl_coordinates_all <- biomaRt::getBM(attributes = attributes,
                                      filters = c("transcript_biotype"),
                                      values = c("protein_coding"),
                                      mart = ensembl)

# Limit to canonical nuclear chromosomes
ensembl_coordinates <- ensembl_coordinates_all |>
  filter(chromosome_name %in% c(1:22, "X", "Y"))

# Add required GTF fields
gtf_data_cds <- ensembl_coordinates |>
  mutate(
    seqname = chromosome_name,
    source = "Ensembl",
    feature = "CDS",
    start = genomic_coding_start,
    end = genomic_coding_end,
    score = ".",
    strand = ifelse(strand == 1, "+", "-"),
    frame = ".",
    attribute = paste0(
      'gene_id "', ensembl_gene_id, '"; ',
      'gene_name "', external_gene_name, '"; '
      )
  ) |>
  filter(!is.na(genomic_coding_start)) |>
  filter(!is.na(genomic_coding_end)) |>
  select(seqname, source, feature, start, end, score, strand, frame, attribute)
gtf_data_utr5 <- ensembl_coordinates |>
  mutate(
    seqname = chromosome_name,
    source = "Ensembl",
    feature = "UTR5",
    start = `5_utr_start`,
    end = `5_utr_end`,
    score = ".",
    strand = ifelse(strand == 1, "+", "-"),
    frame = ".",
    attribute = paste0(
      'gene_id "', ensembl_gene_id, '"; ',
      'gene_name "', external_gene_name, '"; '
    )
  ) |>
  filter(!is.na(`5_utr_start`)) |>
  filter(!is.na(`5_utr_end`)) |>
  select(seqname, source, feature, start, end, score, strand, frame, attribute)
gtf_data_utr3 <- ensembl_coordinates |>
  mutate(
    seqname = chromosome_name,
    source = "Ensembl",
    feature = "UTR3",
    start = `3_utr_start`,
    end = `3_utr_end`,
    score = ".",
    strand = ifelse(strand == 1, "+", "-"),
    frame = ".",
    attribute = paste0(
      'gene_id "', ensembl_gene_id, '"; ',
      'gene_name "', external_gene_name, '"; '
    )
  ) |>
  filter(!is.na(`3_utr_start`)) |>
  filter(!is.na(`3_utr_end`)) |>
  select(seqname, source, feature, start, end, score, strand, frame, attribute)

# Save the results to a file
write.table(rbind(gtf_data_cds, gtf_data_utr5, gtf_data_utr3), gzfile("UTR_CDS_coordinates.gtf.gz"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)

