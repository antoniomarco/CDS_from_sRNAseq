library(pheatmap)
library(patchwork)
library(gprofiler2)
library(tidyverse)
# Bioconductor
library(edgeR)
library(DESeq2)
library(EnhancedVolcano)
# Set seed
set.seed(1221)

# Authors: Aygun Azadova and Antonio Marco

## Load expression read-counts files

# microRNAs
bc_microRNA <-read.table("datasets/bc_patients_microRNA.tab.gz", 
                         header = TRUE, row.names = 1)
bc_microRNA <-bc_microRNA[,-c(1:5)]
colnames(bc_microRNA) <-gsub(".bam","",colnames(bc_microRNA))
# Reorder as per patient/sample
bc_microRNA <- bc_microRNA[,c("SRR7948471", "SRR7948470", "SRR7948473",
                              "SRR7948472", "SRR7948467", "SRR7948466", 
                              "SRR7948469", "SRR7948468", "SRR7948475", 
                              "SRR7948474", "SRR7948482", "SRR7948481", 
                              "SRR7948480", "SRR7948479", "SRR7948478", 
                              "SRR7948477", "SRR7948485", "SRR7948476", 
                              "SRR7948484", "SRR7948483", "SRR7948464",
                              "SRR7948465", "SRR7948462", "SRR7948463")]
# mRNAs
bc_mRNA <-read.table("datasets/bc_patients_mRNA.tab.gz", 
                     header = TRUE, row.names = 1)
bc_mRNA <-bc_mRNA[,-c(1:5)]
colnames(bc_mRNA) <-gsub(".bam","",colnames(bc_mRNA))
# Reorder as per patient/sample
bc_mRNA <- bc_mRNA[,c("SRR7948471", "SRR7948470", "SRR7948473",
                              "SRR7948472", "SRR7948467", "SRR7948466", 
                              "SRR7948469", "SRR7948468", "SRR7948475", 
                              "SRR7948474", "SRR7948482", "SRR7948481", 
                              "SRR7948480", "SRR7948479", "SRR7948478", 
                              "SRR7948477", "SRR7948485", "SRR7948476", 
                              "SRR7948484", "SRR7948483", "SRR7948464",
                              "SRR7948465", "SRR7948462", "SRR7948463")]
# Select only protein-coding genes
coding_genes <- read.table("datasets/gencode.v43.protein_coding.txt.gz")$V1
bc_mRNA <- bc_mRNA[rownames(bc_mRNA) %in% coding_genes,]

## Differential Gene Expression

#Experiment design matrix
cData <- data.frame(sample = c("SRR7948471", "SRR7948470", "SRR7948473",
                               "SRR7948472", "SRR7948467", "SRR7948466", 
                               "SRR7948469", "SRR7948468", "SRR7948475", 
                               "SRR7948474", "SRR7948482", "SRR7948481", 
                               "SRR7948480", "SRR7948479", "SRR7948478", 
                               "SRR7948477", "SRR7948485", "SRR7948476", 
                               "SRR7948484", "SRR7948483", "SRR7948464",
                               "SRR7948465", "SRR7948462", "SRR7948463"),
                    tissue = c("normal", "tumor", "normal", "tumor", "normal",
                               "tumor", "normal", "tumor", "normal", "tumor", 
                               "normal", "tumor", "normal", "tumor", "normal",
                               "tumor", "normal", "tumor", "normal", "tumor", 
                               "normal", "tumor", "normal", "tumor"), 
                    patient = c("p3", "p3", "p5", "p5", "p6", "p6", "p7", "p7",
                                "p14", "p14", "p20", "p20", "p61", "p61", "p67",
                                "p67", "p68", "p68", "p69", "p69", "p74", "p74",
                                "p81", "p81"))

# DGE microRNAs
ddsMatmiR <-DESeqDataSetFromMatrix(countData = bc_microRNA, colData = cData, 
                                   design = ~ patient + tissue)
ddsmiR <-DESeq(ddsMatmiR)
resmiR <- results(ddsmiR, contrast=c("tissue","tumor","normal"))
hits_over_microRNA<-resmiR[which((resmiR$padj<=0.01) & (abs(resmiR$log2FoldChange) >= 1)),]
# plots
volcano_microRNA <- EnhancedVolcano(resmiR,
                lab = row.names(resmiR),
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.01,
                xlim = c(-5,5),
                labSize = 3)
vsdmiR <-varianceStabilizingTransformation(ddsmiR, blind=FALSE)
hit_genes_microRNA <- rownames(hits_over_microRNA)
vsdmiR_hits <- vsdmiR[hit_genes_microRNA,]
heatmap_microRNA <- pheatmap(as.matrix(assay(vsdmiR_hits)),
                             fontsize_row = 6,
                             labels_col = rep(c("N", "T"), 12))
ggsave(volcano_microRNA, file = "plots/volcano_miRNA.png")
ggsave(heatmap_microRNA, file = "plots/heatmap_miRNA.png")


# DGE mRNAs
ddsMatmR <-DESeqDataSetFromMatrix(countData = bc_mRNA, colData = cData, 
                                  design = ~ patient + tissue)
ddsmR <-DESeq(ddsMatmR)
resmR <- results(ddsmR, contrast=c("tissue","tumor","normal"))
hits_over_mRNA<-resmR[which((resmR$padj<=0.01) & (abs(resmR$log2FoldChange) >= 1)),]
# plots
volcano_mRNA <- EnhancedVolcano(resmR,
                                lab = row.names(resmR),
                                x = "log2FoldChange",
                                y = "padj",
                                pCutoff = 0.01,
                                xlim = c(-5,5),
                                labSize = 3)
ggsave(volcano_mRNA, file = "plots/volcano_mRNA.png")

vsdmR <-varianceStabilizingTransformation(ddsmR, blind=FALSE)
hit_genes_mRNA <- rownames(hits_over_mRNA)
vsdmR_hits <- vsdmR[hit_genes_mRNA,]
heatmap_mRNA <- pheatmap(as.matrix(assay(vsdmR_hits)),
                         fontsize_row = 6,
                         labels_col = rep(c("N", "T"), 12),
                         show_rownames = FALSE) # no names
ggsave(heatmap_mRNA, file = "plots/heatmap_mRNA_nonames.png")

heatmap_mRNA <- pheatmap(as.matrix(assay(vsdmR_hits)),
                         fontsize_row = 6,
                         labels_col = rep(c("N", "T"), 12))
ggsave(heatmap_mRNA, file = "plots/heatmap_mRNA.png", units = "cm", width = 12, height = 96)





## Functional annotation of DEG mRNA. For information only
gostres_up <- gost(query = rownames(resmR[which((resmR$padj<=0.01) & (resmR$log2FoldChange >= 1)),]),
                   organism = "hsapiens",
                   user_threshold = 0.05,
                   correction_method = "fdr",
                   domain_scope = "custom",
                   custom_bg = rownames(resmR))
gostplot(gostres_up, capped = TRUE, interactive = TRUE)
gostres_down <- gost(query = rownames(resmR[which((resmR$padj<=0.01) & (resmR$log2FoldChange <= -1)),]),
                   organism = "hsapiens",
                   user_threshold = 0.05,
                   correction_method = "fdr",
                   domain_scope = "custom",
                   custom_bg = rownames(resmR))
gostplot(gostres_down, capped = TRUE, interactive = TRUE)
write.table(resmR, file = "plots/DGE_mRNA_BC.tab")

## Compare log2FC DGE and canonical microRNA target sites

# # Get 3'UTRs of DGE and canonical microRNA targets for DGE genes
# # Slow and unstable connection to ENSEMBL via biomaRt
# # Comment out to run, or leave it to load pre-computed file
# library(biomaRt)
# library(tidyverse)
# # requires local installation of seedVicious (https://seedvicious.essex.ac.uk)
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# utr3 <- data.frame(gene = c(), utr3 = c())
# count = 0
# for(g in hit_genes_mRNA){
#   count = count + 1
#   cat(paste(count, " out of ", length(hit_genes_mRNA), "\r", sep = ""))
#   seq <- getSequence(id = g, type = "hgnc_symbol", seqType = "3utr", mart = ensembl)
#   utr3 <- rbind(utr3, c(gene = g,utr3 = seq[which.max(nchar(seq$`3utr`)), ]))
# }
# # temp transcripts
# utr3$gene <- paste(">", utr3$gene, sep = "")
# utr3 <- utr3[,-3]
# utr3 <- utr3[!grepl("Sequence unavailable", utr3$utr3.3utr),]
# write.table(utr3, file="temp_mRNA.fas", sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
# # temp microRNAs (mature)
# system("wget https://www.mirbase.org/download/mature.fa", intern = FALSE, ignore.stderr = TRUE)
# system("cat mature.fa | grep -A 1 -e '^>hsa' | grep -v -e '--' | awk '{print $1}' > hsa_mature.fa", intern = FALSE, ignore.stderr = TRUE)
# # Targets with seedVicious
# sv_run <- "seedVicious -i temp_mRNA.fas -m hsa_mature.fa -o temp_out"
# system(sv_run, intern = FALSE, ignore.stderr = TRUE)
# targets_SV <- read.table("temp_out", sep = "\t", h = T)
# # Remove temp files
# system("rm temp_mRNA.fas temp_out mature.fa hsa_mature.fa")
# save(targets_SV, file = "datasets/targets_SV.Rdata")
load("datasets/targets_SV.Rdata")

# Load MTB data
targets_MTB <- read.table("datasets/mirtarbase_hsa_9.tab.gz")
colnames(targets_MTB) <- c("miR", "tr")
# log2FC MTB of interacting pairs
MTB_pairs_l2FC <- as_tibble(hits_over_microRNA, rownames = "miR") |> 
  right_join(targets_MTB, by = "miR") |> 
  select (miR, log2FoldChange, tr) |>
  left_join(as_tibble(hits_over_mRNA, rownames = "tr"), by = "tr") |>
  select(miR, log2FoldChange.x, tr, log2FoldChange.y) |> 
  filter((!is.na(log2FoldChange.x)) & (!is.na(log2FoldChange.y))) |>
  unique()


## Analysis of targets (MTB and SV)

# log2FC SV of interacting pairs
SV_pairs_l2FC <- as_tibble(hits_over_microRNA, rownames = "miR") |> 
  right_join(targets_SV, by = "miR") |> 
  select (miR, log2FoldChange, tr) |>
  left_join(as_tibble(hits_over_mRNA, rownames = "tr"), by = "tr") |>
  select(miR, log2FoldChange.x, tr, log2FoldChange.y) |> 
  filter((!is.na(log2FoldChange.x)) & (!is.na(log2FoldChange.y)))

# Prefilter in dataset: Not in target
targets_SV_f <- targets_SV[which((targets_SV$miR %in% row.names(hits_over_microRNA)) & (targets_SV$tr %in% row.names(hits_over_mRNA))),]
targets_SV_list <- paste0(targets_SV_f$miR, targets_SV_f$tr)
SV_pairs_l2FC_NT <- data.frame(miR = character(0),
                               log2FoldChange.x = numeric(0),
                               tr = character(0),
                               log2FoldChange.y = numeric(0))
for(miR in 1:nrow(hits_over_microRNA)){
  for(tr in 1:nrow(hits_over_mRNA)){
    if(paste0(rownames(hits_over_microRNA[miR,]), rownames(hits_over_mRNA[tr,])) %in% targets_SV_list){
    }else{
      SV_pairs_l2FC_NT <- rbind(SV_pairs_l2FC_NT, data.frame(miR = rownames(hits_over_microRNA[miR,]),
                                                             log2FoldChange.x = hits_over_microRNA[miR,'log2FoldChange'],
                                                             tr = rownames(hits_over_mRNA[tr,]),
                                                             log2FoldChange.y = hits_over_mRNA[tr,'log2FoldChange']))
    }
  }
}
SV_pairs_l2FC_NT <- as_tibble(SV_pairs_l2FC_NT)

# For contingency table
in_MTB_pn <- MTB_pairs_l2FC |> 
  unique() |> 
  filter((log2FoldChange.x > 0) & (log2FoldChange.y < 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
in_SV_pn <- SV_pairs_l2FC |> 
  unique() |> 
  filter((log2FoldChange.x > 0) & (log2FoldChange.y < 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
not_in_SV_pn <- SV_pairs_l2FC_NT |> 
  unique() |> 
  filter((log2FoldChange.x > 0) & (log2FoldChange.y < 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
in_MTB_nn <- MTB_pairs_l2FC |> 
  unique() |> 
  filter((log2FoldChange.x < 0) & (log2FoldChange.y < 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
in_SV_nn <- SV_pairs_l2FC |> 
  unique() |> 
  filter((log2FoldChange.x < 0) & (log2FoldChange.y < 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
not_in_SV_nn <- SV_pairs_l2FC_NT |> 
  unique() |> 
  filter((log2FoldChange.x < 0) & (log2FoldChange.y < 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
in_MTB_pp <- MTB_pairs_l2FC |> 
  unique() |> 
  filter((log2FoldChange.x > 0) & (log2FoldChange.y > 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
in_SV_pp <- SV_pairs_l2FC |> 
  unique() |> 
  filter((log2FoldChange.x > 0) & (log2FoldChange.y > 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
not_in_SV_pp <- SV_pairs_l2FC_NT |> 
  unique() |> 
  filter((log2FoldChange.x > 0) & (log2FoldChange.y > 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
in_MTB_np <- MTB_pairs_l2FC |> 
  unique() |> 
  filter((log2FoldChange.x < 0) & (log2FoldChange.y > 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
in_SV_np <- SV_pairs_l2FC |> 
  unique() |> 
  filter((log2FoldChange.x < 0) & (log2FoldChange.y > 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)
not_in_SV_np <- SV_pairs_l2FC_NT |> 
  unique() |> 
  filter((log2FoldChange.x < 0) & (log2FoldChange.y > 0)) |>
  mutate(tar = paste0(tr, miR)) |>
  select(tar)

# Contingecy table
in_SV_pn_in_TB <- in_SV_pn[in_SV_pn$tar %in% in_MTB_pn$tar,]
in_SV_pn_not_in_TB <- in_SV_pn[!(in_SV_pn$tar %in% in_MTB_pn$tar),]
not_in_SV_pn_in_TB <- not_in_SV_pn[not_in_SV_pn$tar %in% in_MTB_pn$tar,]
not_in_SV_pn_not_in_TB <- not_in_SV_pn[!(not_in_SV_pn$tar %in% in_MTB_pn$tar),]
in_SV_nn_in_TB <- in_SV_nn[in_SV_nn$tar %in% in_MTB_nn$tar,]
in_SV_nn_not_in_TB <- in_SV_nn[!(in_SV_nn$tar %in% in_MTB_nn$tar),]
not_in_SV_nn_in_TB <- not_in_SV_nn[not_in_SV_nn$tar %in% in_MTB_nn$tar,]
not_in_SV_nn_not_in_TB <- not_in_SV_nn[!(not_in_SV_nn$tar %in% in_MTB_nn$tar),]
in_SV_pp_in_TB <- in_SV_pp[in_SV_pp$tar %in% in_MTB_pp$tar,]
in_SV_pp_not_in_TB <- in_SV_pp[!(in_SV_pp$tar %in% in_MTB_pp$tar),]
not_in_SV_pp_in_TB <- not_in_SV_pp[not_in_SV_pp$tar %in% in_MTB_pp$tar,]
not_in_SV_pp_not_in_TB <- not_in_SV_pp[!(not_in_SV_pp$tar %in% in_MTB_pp$tar),]
in_SV_np_in_TB <- in_SV_np[in_SV_np$tar %in% in_MTB_np$tar,]
in_SV_np_not_in_TB <- in_SV_np[!(in_SV_np$tar %in% in_MTB_np$tar),]
not_in_SV_np_in_TB <- not_in_SV_np[not_in_SV_np$tar %in% in_MTB_np$tar,]
not_in_SV_np_not_in_TB <- not_in_SV_np[!(not_in_SV_np$tar %in% in_MTB_np$tar),]


# Enrichment in targets (predicted) for log fold configurations against MTB
LOR_pn_nn <- log((nrow(in_SV_pn_in_TB)*nrow(in_SV_nn_not_in_TB))/(nrow(in_SV_pn_not_in_TB)*nrow(in_SV_nn_in_TB)))
SEM_pn_nn <- sqrt(nrow(in_SV_pn_in_TB)^(-1) + nrow(in_SV_nn_not_in_TB)^(-1) + nrow(in_SV_pn_not_in_TB)^(-1) + nrow(in_SV_nn_in_TB)^(-1))
ZSC_pn_nn <- LOR_pn_nn / SEM_pn_nn
PVL_pn_nn <- pnorm(ZSC_pn_nn, mean = 0, sd = 1, lower.tail = FALSE) # One-tailed
#
# np vs pp
LOR_np_pp <- log((nrow(in_SV_np_in_TB)*nrow(in_SV_pp_not_in_TB))/(nrow(in_SV_np_not_in_TB)*nrow(in_SV_pp_in_TB)))
SEM_np_pp <- sqrt(nrow(in_SV_np_in_TB)^(-1) + nrow(in_SV_pp_not_in_TB)^(-1) + nrow(in_SV_np_not_in_TB)^(-1) + nrow(in_SV_pp_in_TB)^(-1))
ZSC_np_pp <- LOR_np_pp / SEM_np_pp
PVL_np_pp <- pnorm(ZSC_np_pp, mean = 0, sd = 1, lower.tail = FALSE) # One-tailed
#
# np vs pp
LOR_nn_pp <- log((nrow(in_SV_nn_in_TB)*nrow(in_SV_pp_not_in_TB))/(nrow(in_SV_nn_not_in_TB)*nrow(in_SV_pp_in_TB)))
SEM_nn_pp <- sqrt(nrow(in_SV_nn_in_TB)^(-1) + nrow(in_SV_pp_not_in_TB)^(-1) + nrow(in_SV_nn_not_in_TB)^(-1) + nrow(in_SV_pp_in_TB)^(-1))
ZSC_nn_pp <- LOR_nn_pp / SEM_nn_pp
PVL_nn_pp <- pnorm(ZSC_nn_pp, mean = 0, sd = 1, lower.tail = FALSE) # One-tailed

out_table <- tibble(numerator = c("miR up; tr down", "miR down; tr up", "miR down tr down"),
       denominator = c("miR down; tr down", "miR up; tr up", "miR up tr up"),
       odds_ratio = c(exp(1)^LOR_pn_nn, exp(1)^LOR_np_pp, exp(1)^LOR_nn_pp),
       z_score = c(ZSC_pn_nn, ZSC_np_pp, ZSC_nn_pp),
       p_value = c(PVL_pn_nn, PVL_np_pp, PVL_nn_pp))

write_tsv(out_table, file = "plots/target_enrichment_BCa.tab")



# From paper associated to datasets
# https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-5300-6
# To assess the physiological and pathological significance of the observed
# down-regulation of miR-10b, we measured by qRTPCR the mRNA levels of 15 
# previously validated targets of miR-10b (primer sequences in 
# Additional file 1: Table S2) based on the validated miRNA target database,
# miRTarBase [25], in all the RNA samples from our 83-patient cohort. We 
# identified that miR-10b levels showed a significant inverse correlation with
# target mRNA levels in at least one subset of samples (the tumor or the benign)
# for 4 out of 15 genes tested: MAPRE1, PIEZO1, SRSF1 and TP53 (Fig. 3a-c).
# Note, not successful: BCL2L11, BDNF, CDKN1A, CDKN2A, HOXD10, KLF4, NCOR2,
# PAX6, PPARa, PTEN, TRA2B

gold <- c(rep(1,4), rep(0, 10))

test <- resmR |> 
  as.data.frame() |>
  rownames_to_column('tr') |>
  as_tibble() |>
  filter(tr %in% c("MAPRE1", "PIEZO1", "SRSF1", "TP53",
                   "BCL2L11", "BDNF",  "CDKN1A", "CDKN2A", 
                   "HOXD10", "KLF4", "NCOR2",
                   "PPARA", "PTEN", "TRA2B")) |> # PAX6 removed as expression was not detected
  select(log2FoldChange) |>
  mutate(log2FoldChange = ifelse(log2FoldChange > 0, 1, 0)) |>
  unlist() |>
  as.vector()


# Check expression of miR-10a-5p
resmiR |> 
  as.data.frame() |>
  rownames_to_column('tr') |>
  as_tibble() |>
  filter(tr %in% c("hsa-miR-10b-5p"))

TP = length(which((test + gold) == 2))
FP = length(which((test[5:14] + gold[5:14]) == 1))
TN = length(which((test + gold) == 0))
FN = length(which((test[1:4] + gold[1:4]) == 1))

recall <- TP/(TP+FN)
FPR <- FP/(FP+TN)
precision <- TP/(TP+FP)
accuracy <- (TP+TN)/(TP+TN+FP+FN)

# Preprare output
output_lines <- list(
  TP = TP, FP = FP, TN = TN, FN = FN,
  recall = recall, FPR = FPR,
  precision = precision, accuracy = accuracy
)
formatted_lines <- sapply(names(output_lines), function(name) {
  paste0(name, " = ", output_lines[[name]])
})
# Write output
writeLines(formatted_lines, "plots/mir10-accuracy_output.txt")



# Plot
custom_order <- c("BCL2L11", "BDNF", 
                  "CDKN1A", "CDKN2A", "HOXD10", "KLF4", "NCOR2", "PAX6",
                  "PPARA", "PTEN", "TRA2B",
                  "TP53","MAPRE1", "PIEZO1", "SRSF1")

recall_plot <- resmR |> 
  as.data.frame() |>
  rownames_to_column('tr') |>
  as_tibble() |>
  filter(tr %in% custom_order) |>
  mutate(tr = factor(tr, levels = custom_order)) |> 
  arrange(tr) |>
  ggplot(aes(x = tr, y = log2FoldChange)) +
  geom_col() +
  coord_flip() +  # ← make it horizontal
  theme_minimal() +
  labs(x = "Gene", y = "log2 Fold Change", title = "Horizontal Bar Plot")
ggsave(recall_plot, file = "plots/reacall_plot.png")

# exit
q()
