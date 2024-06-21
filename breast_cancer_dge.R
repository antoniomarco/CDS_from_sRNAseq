library(pheatmap)
library(patchwork)
library(gprofiler2)
library(tidyverse)
# Bioconductor
library(edgeR)
library(DESeq2)
library(EnhancedVolcano)

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
vsdmR <-varianceStabilizingTransformation(ddsmR, blind=FALSE)
hit_genes_mRNA <- rownames(hits_over_mRNA)
vsdmR_hits <- vsdmR[hit_genes_mRNA,]
heatmap_mRNA <- pheatmap(as.matrix(assay(vsdmR_hits)),
                         fontsize_row = 6,
                         labels_col = rep(c("N", "T"), 12))

ggsave(volcano_mRNA, file = "plots/volcano_mRNA.png")
ggsave(heatmap_mRNA, file = "plots/heatmap_mRNA.png")


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


## Analysis of targets (MTB and SV)

# log2FC SV of interacting pairs
SV_pairs_l2FC <- as_tibble(hits_over_microRNA, rownames = "miR") |> 
  right_join(targets_SV, by = "miR") |> 
  select (miR, log2FoldChange, tr) |>
  left_join(as_tibble(hits_over_mRNA, rownames = "tr"), by = "tr") |>
  select(miR, log2FoldChange.x, tr, log2FoldChange.y) |> 
  filter((!is.na(log2FoldChange.x)) & (!is.na(log2FoldChange.y)))
# Plot
SV_nn <- SV_pairs_l2FC |> unique() |> filter((log2FoldChange.x < 0) & (log2FoldChange.y < 0)) |> nrow()
SV_pn <- SV_pairs_l2FC |> unique() |> filter((log2FoldChange.x > 0) & (log2FoldChange.y < 0)) |> nrow()
SV_np <- SV_pairs_l2FC |> unique() |> filter((log2FoldChange.x < 0) & (log2FoldChange.y > 0)) |> nrow()
SV_pp <- SV_pairs_l2FC |> unique() |> filter((log2FoldChange.x > 0) & (log2FoldChange.y > 0)) |> nrow()

# # 1 2 or 3 or more targets code only for information
# # Uncomment if wanted
# # log2FC SV=1 (exactly one) of interacting pairs
# SV_pairs_l2FC_unique <- SV_pairs_l2FC |>
#   distinct()
# SV_pairs_l2FC_unique <- SV_pairs_l2FC_unique |>
#   mutate(count = SV_pairs_l2FC |>
#            count(miR, tr) |>
#            pull(n))
# targets_SV_unique <- targets_SV |>
#   select(miR,tr) |>
#   distinct()
# targets_SV1 <- as_tibble(targets_SV_unique) |>
#   mutate(count = targets_SV |>
#            count(miR, tr) |>
#            pull(n)) |>
#   filter(count == 1) |>
#   select(miR, tr)
# SV1_pairs_l2FC <- as_tibble(hits_over_microRNA, rownames = "miR") |> 
#   right_join(targets_SV1, by = "miR") |> 
#   select (miR, log2FoldChange, tr) |>
#   left_join(as_tibble(hits_over_mRNA, rownames = "tr"), by = "tr") |>
#   select(miR, log2FoldChange.x, tr, log2FoldChange.y) |> 
#   filter((!is.na(log2FoldChange.x)) & (!is.na(log2FoldChange.y)))
# # Plot
# SV1_pairs_l2FC |> unique() |> filter((log2FoldChange.x) < 0 & (log2FoldChange.y < 0)) |> nrow()
# SV1_pairs_l2FC |> unique() |> filter((log2FoldChange.x) > 0 & (log2FoldChange.y < 0)) |> nrow()
# SV1_pairs_l2FC |> unique() |> filter((log2FoldChange.x) < 0 & (log2FoldChange.y > 0)) |> nrow()
# SV1_pairs_l2FC |> unique() |> filter((log2FoldChange.x) > 0 & (log2FoldChange.y > 0)) |> nrow()
# 
# # log2FC SV=2 (two) of interacting pairs
# targets_SV2 <- as_tibble(targets_SV_unique) |>
#   mutate(count = targets_SV |>
#            count(miR, tr) |>
#            pull(n)) |>
#   filter(count == 2) |>
#   select(miR, tr)
# SV2_pairs_l2FC <- as_tibble(hits_over_microRNA, rownames = "miR") |> 
#   right_join(targets_SV2, by = "miR") |> 
#   select (miR, log2FoldChange, tr) |>
#   left_join(as_tibble(hits_over_mRNA, rownames = "tr"), by = "tr") |>
#   select(miR, log2FoldChange.x, tr, log2FoldChange.y) |> 
#   filter((!is.na(log2FoldChange.x)) & (!is.na(log2FoldChange.y)))
# # Plot
# SV2_pairs_l2FC |> unique() |> filter((log2FoldChange.x) < 0 & (log2FoldChange.y < 0)) |> nrow()
# SV2_pairs_l2FC |> unique() |> filter((log2FoldChange.x) > 0 & (log2FoldChange.y < 0)) |> nrow()
# SV2_pairs_l2FC |> unique() |> filter((log2FoldChange.x) < 0 & (log2FoldChange.y > 0)) |> nrow()
# SV2_pairs_l2FC |> unique() |> filter((log2FoldChange.x) > 0 & (log2FoldChange.y > 0)) |> nrow()
# 
# # log2FC SV>2 (more than two) of interacting pairs
# targets_SV3 <- as_tibble(targets_SV_unique) |>
#   mutate(count = targets_SV |>
#            count(miR, tr) |>
#            pull(n)) |>
#   filter(count > 2) |>
#   select(miR, tr)
# SV3_pairs_l2FC <- as_tibble(hits_over_microRNA, rownames = "miR") |> 
#   right_join(targets_SV3, by = "miR") |> 
#   select (miR, log2FoldChange, tr) |>
#   left_join(as_tibble(hits_over_mRNA, rownames = "tr"), by = "tr") |>
#   select(miR, log2FoldChange.x, tr, log2FoldChange.y) |> 
#   filter((!is.na(log2FoldChange.x)) & (!is.na(log2FoldChange.y)))
# # Plot
# SV3_pairs_l2FC |> unique() |> filter((log2FoldChange.x) < 0 & (log2FoldChange.y < 0)) |> nrow()
# SV3_pairs_l2FC |> unique() |> filter((log2FoldChange.x) > 0 & (log2FoldChange.y < 0)) |> nrow()
# SV3_pairs_l2FC |> unique() |> filter((log2FoldChange.x) < 0 & (log2FoldChange.y > 0)) |> nrow()
# SV3_pairs_l2FC |> unique() |> filter((log2FoldChange.x) > 0 & (log2FoldChange.y > 0)) |> nrow()

# MTB
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
# Plot
MTB_nn <- MTB_pairs_l2FC |> filter((log2FoldChange.x < 0) & (log2FoldChange.y < 0)) |> nrow()
MTB_pn <- MTB_pairs_l2FC |> filter((log2FoldChange.x > 0) & (log2FoldChange.y < 0)) |> nrow()
MTB_np <- MTB_pairs_l2FC |> filter((log2FoldChange.x < 0) & (log2FoldChange.y > 0)) |> nrow()
MTB_pp <- MTB_pairs_l2FC |> filter((log2FoldChange.x > 0) & (log2FoldChange.y > 0)) |> nrow()


# Plots
MTB <- c(MTB_nn,MTB_pn,MTB_np,MTB_pp)
MTB_plot <- ggplot() + 
  geom_bar(aes(x=seq_along(MTB),y=MTB), stat='identity') +
  xlab('') +
  ylab('Number of miR/mR interactions') +
  annotate("text", x=1, y=(MTB_nn/2)+10, label= paste("\U2193", "miR", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=1, y=(MTB_nn/2)-10, label= paste("\U2193", "gene", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=2, y=(MTB_pn/2)+10, label= paste("\U2191", "miR", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=2, y=(MTB_pn/2)-10, label= paste("\U2193", "gene", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=3, y=(MTB_np/2)+10, label= paste("\U2191", "miR", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=3, y=(MTB_np/2)-10, label= paste("\U2191", "gene", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=4, y=(MTB_pp/2)+10, label= paste("\U2193", "miR", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=4, y=(MTB_pp/2)-10, label= paste("\U2191", "gene", sep = " "), size = 3.4, color = "white") +
  theme_bw()
SV <- c(SV_nn,SV_pn,SV_np,SV_pp)
SV_plot <- ggplot() + 
  geom_bar(aes(x=seq_along(SV),y=SV), stat='identity') +
  xlab('') +
  ylab('Number of miR/mR interactions') +
  annotate("text", x=1, y=(SV_nn/2)+100, label= paste("\U2193", "miR", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=1, y=(SV_nn/2)-100, label= paste("\U2193", "gene", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=2, y=(SV_pn/2)+100, label= paste("\U2191", "miR", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=2, y=(SV_pn/2)-100, label= paste("\U2193", "gene", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=3, y=(SV_np/2)+100, label= paste("\U2191", "miR", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=3, y=(SV_np/2)-100, label= paste("\U2191", "gene", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=4, y=(SV_pp/2)+100, label= paste("\U2193", "miR", sep = " "), size = 3.4, color = "white") +
  annotate("text", x=4, y=(SV_pp/2)-100, label= paste("\U2191", "gene", sep = " "), size = 3.4, color = "white") +
  theme_bw()
ggsave(MTB_plot + SV_plot + plot_annotation(tag_levels = "A"), file = "plots/figure_bars_targets.png", height = 3, width = 6)
# Shape asexpected!!

# From paper associated to datasets
# https://bmccancer.biomedcentral.com/articles/10.1186/s12885-019-5300-6
# The paper extracts the targets from MTB as well.
# They validate (up reg and mir10b goes down in cancer)  these: 
# "We identified that miR-10b levels showed a significant inverse correlation with
# target mRNA levels in at least one subset of samples (the tumor or the benign) 
# for 4 out of 15 genes tested: MAPRE1, PIEZO1, SRSF1 and TP53"
MTB_pairs_l2FC |> filter(miR == "hsa-miR-10b-5p") |> filter(log2FoldChange.y > 0) |> pull(tr)
MTB_pairs_l2FC |> filter(miR == "hsa-miR-21-5p") |> filter(log2FoldChange.y < 0) |> pull(tr)

MTB_pairs_l2FC |> filter(log2FoldChange.x >= max(MTB_pairs_l2FC$log2FoldChange.x))
# hsa-miR-429 COL4A1 DLC1 QKI
# DLC1 seems to be a tumor-supp gene
# https://pubmed.ncbi.nlm.nih.gov/?term=DLC1+cancer
# and check QK1 in Breast cancer
# https://pubmed.ncbi.nlm.nih.gov/33569406/

MTB_pairs_l2FC |> filter(log2FoldChange.x <= min(MTB_pairs_l2FC$log2FoldChange.x))
# hsa-miR-144-5p, RUNX1
# https://pubmed.ncbi.nlm.nih.gov/?term=RUNX1+cancer

# Panel of co-expression for all four pairs above
# transform
bc_microRNA_voom <- voom(bc_microRNA)$E[c('hsa-miR-429','hsa-miR-144-5p'),]
bc_mRNA_voom <- voom(bc_mRNA)$E[c('COL4A1','DLC1','QKI','RUNX1'),]
bc_voom <- rbind(bc_microRNA_voom,bc_mRNA_voom)
# Re-order columns so first all normal and after all tumor samples
bc_voom <- bc_voom[,c(seq(1,24,2),seq(2,24,2))]

# Plots
plot_ex_1 <- as_tibble(cbind(sample = colnames(bc_voom), t(bc_voom))) |>
  mutate_at(2:7, as.numeric) |>
  rename_with( ~ gsub("-", "_", .x, fixed = TRUE)) |>
  mutate(tissue = c(rep("normal",12),rep("tumor", 12))) |>
  ggplot(aes(x = hsa_miR_429, y = COL4A1, color = tissue)) +
  geom_point() +
  theme_bw()
plot_ex_2 <-as_tibble(cbind(sample = colnames(bc_voom), t(bc_voom))) |>
  mutate_at(2:7, as.numeric) |>
  rename_with( ~ gsub("-", "_", .x, fixed = TRUE)) |>
  mutate(tissue = c(rep("normal",12),rep("tumor", 12))) |>
  ggplot(aes(x = hsa_miR_429, y = DLC1, color = tissue)) +
  geom_point()  +
  theme_bw()
plot_ex_3 <-as_tibble(cbind(sample = colnames(bc_voom), t(bc_voom))) |>
  mutate_at(2:7, as.numeric) |>
  rename_with( ~ gsub("-", "_", .x, fixed = TRUE)) |>
  mutate(tissue = c(rep("normal",12),rep("tumor", 12))) |>
  ggplot(aes(x = hsa_miR_429, y = QKI, color = tissue)) +
  geom_point()  +
  theme_bw()
plot_ex_4 <-as_tibble(cbind(sample = colnames(bc_voom), t(bc_voom))) |>
  mutate_at(2:7, as.numeric) |>
  rename_with( ~ gsub("-", "_", .x, fixed = TRUE)) |>
  mutate(tissue = c(rep("normal",12),rep("tumor", 12))) |>
  ggplot(aes(x = hsa_miR_144_5p, y = RUNX1, color = tissue)) +
  geom_point()  +
  theme_bw()

ggsave((plot_ex_1 + plot_ex_2) / (plot_ex_3 + plot_ex_4) + plot_annotation(tag_levels = "A"), file = "plots/figure_examples.png")

# exit
q()
