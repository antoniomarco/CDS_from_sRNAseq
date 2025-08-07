library(tidyverse)
library(patchwork)
library(WebGestaltR)
# Bioconductor
# BiocManager::install("biomaRt")
library(biomaRt)
# BiocManager::install("edgeR")
library(edgeR)

# Small RNA mapped to all genes
mRNA_expression_small <- read.table("datasets/Meunier_mRNA.tab.gz", header = TRUE, row.names = 1)
mRNA_expression_small <- mRNA_expression_small[,-c(1:5)]
colnames(mRNA_expression_small) <- c("brain_small", "cerebellum_small", "heart_small", "kidney_small", "testis_small")

# Total RNA mapped to all genes
mRNA_expression_total <- read.table("datasets/Brawand_mRNA.tab.gz", header = TRUE, row.names = 1)
mRNA_expression_total <- mRNA_expression_total[,-c(1:5)]
colnames(mRNA_expression_total) <- c("brain_total", "cerebellum_1_total", "cerebellum_2_total", "heart_1_total", "heart_2_total", "kidney_total")

# Voom transformation
mRNA_expression_small_voom <- as.data.frame(voom(mRNA_expression_small)$E)
mRNA_expression_total_voom <- as.data.frame(voom(mRNA_expression_total)$E)
mRNA_expression <- merge(mRNA_expression_total_voom, mRNA_expression_small_voom, by = 0)
row.names(mRNA_expression) <- mRNA_expression$Row.names
mRNA_expression <- mRNA_expression[,-1]

# Consider only Coding Sequences
CDS_names <- read.table("datasets/gencode.v43.protein_coding.txt.gz")$V1
mRNA_expression_CDS <- mRNA_expression[CDS_names,]
# Fits (linear models)
brain.cds.lm <- lm(mRNA_expression_CDS$brain_total ~ mRNA_expression_CDS$brain_small)
cerebellum.cds.lm <- lm(mRNA_expression_CDS$cerebellum_1_total ~ mRNA_expression_CDS$cerebellum_small)
heart.cds.lm <- lm(mRNA_expression_CDS$heart_1_total ~ mRNA_expression_CDS$heart_small)
kidney.cds.lm <- lm(mRNA_expression_CDS$kidney_total ~ mRNA_expression_CDS$kidney_small)

brain_CDS_plot <- ggplot(data = mRNA_expression_CDS,
       aes(x = brain_small,
           y = brain_total)) +
       ggtitle("Brain") +
       xlab("smallRNA-seq Expression level") +
       ylab("RNA-seq Expression level") +
       geom_point(alpha=0.4, col = "darkgrey") +
       geom_abline(slope = coef(brain.cds.lm)[["mRNA_expression_CDS$brain_small"]], 
                   intercept = coef(brain.cds.lm)[["(Intercept)"]],
                   col = 'darkred') +
       annotate(geom = 'text', x = 0, y = 13, size = 3, label = paste("R^2 ==", round(summary(brain.cds.lm)$adj.r.squared, 2)), parse = TRUE) +
       theme_light() +
       theme(plot.title = element_text(hjust = 0.5, size = 11)) +
       theme(axis.title = element_text(size = 9)) 

 cerebellum_CDS_plot <- ggplot(data = mRNA_expression_CDS,
       aes(x = cerebellum_small,
           y = cerebellum_1_total)) +
       ggtitle("Cerebellum") +
       xlab("smallRNA-seq Expression level") +
       ylab("RNA-seq Expression level") +
       geom_point(alpha=0.4, col = "darkgrey") +
       geom_abline(slope = coef(cerebellum.cds.lm)[["mRNA_expression_CDS$cerebellum_small"]], 
                   intercept = coef(cerebellum.cds.lm)[["(Intercept)"]],
                   col = 'darkred') +
       annotate(geom = 'text', x = 0, y = 13, size = 3, label = paste("R^2 ==", round(summary(cerebellum.cds.lm)$adj.r.squared, 2)), parse = TRUE) +
       theme_light() +    
       theme(plot.title = element_text(hjust = 0.5, size = 11)) +
       theme(axis.title = element_text(size = 9)) 
 
heart_CDS_plot <- ggplot(data = mRNA_expression_CDS,
       aes(x = heart_small,
           y = heart_1_total)) +
       ggtitle("Heart") +
       xlab("smallRNA-seq Expression level") +
       ylab("RNA-seq Expression level") +
       geom_point(alpha=0.4, col = "darkgrey") +
       geom_abline(slope = coef(heart.cds.lm)[["mRNA_expression_CDS$heart_small"]], 
                   intercept = coef(heart.cds.lm)[["(Intercept)"]],
                   col = 'darkred') +
       annotate(geom = 'text', x = -.5, y = 13, size = 3, label = paste("R^2 ==", round(summary(heart.cds.lm)$adj.r.squared, 2)), parse = TRUE) +
       theme_light() +    
       theme(plot.title = element_text(hjust = 0.5, size = 11)) +
       theme(axis.title = element_text(size = 9)) 

kidney_CDS_plot <- ggplot(data = mRNA_expression_CDS,
       aes(x = kidney_small,
           y = kidney_total)) +
       ggtitle("Kidney") +
       xlab("smallRNA-seq Expression level") +
       ylab("RNA-seq Expression level") +
       geom_point(alpha=0.4, col = "darkgrey") +
       geom_abline(slope = coef(kidney.cds.lm)[["mRNA_expression_CDS$kidney_small"]], 
                   intercept = coef(kidney.cds.lm)[["(Intercept)"]],
                   col = 'darkred') +
       annotate(geom = 'text', x = -2, y = 13, size = 3, label = paste("R^2 ==", round(summary(kidney.cds.lm)$adj.r.squared, 2)), parse = TRUE) +
       theme_light() +       
       theme(plot.title = element_text(hjust = 0.5, size = 11)) +
       theme(axis.title = element_text(size = 9)) 

# (brain_CDS_plot | cerebellum_CDS_plot)/(heart_CDS_plot | kidney_CDS_plot) + plot_annotation(tag_levels = "A")
ggsave((brain_CDS_plot | cerebellum_CDS_plot)/(heart_CDS_plot | kidney_CDS_plot) + plot_annotation(tag_levels = "A"), file = "plots/Figure_1.png")


## Linear models for Mouse datasets
#
# Small RNA mapped to all genes
mRNA_expression_small_mmu <- read.table("datasets/Meunier_mRNA_mmu.tab.gz", header = TRUE, row.names = 1)
mRNA_expression_small_mmu <- mRNA_expression_small_mmu[,-c(1:5)]
colnames(mRNA_expression_small_mmu) <- c("brain_small", "cerebellum_small", "heart_small", "kidney_small", "testis_small")

# Total RNA mapped to all genes
mRNA_expression_total_mmu <- read.table("datasets/Brawand_mRNA_mmu.tab.gz", header = TRUE, row.names = 1)
mRNA_expression_total_mmu <- mRNA_expression_total_mmu[,-c(1:5)]
colnames(mRNA_expression_total_mmu) <- c("brain_1_total", "brain_2_total", "brain_3_total", "cerebellum_total", "heart_total", "kidney_total", "testis_total")

# Voom transformation
mRNA_expression_small_voom_mmu <- as.data.frame(voom(mRNA_expression_small_mmu)$E)
mRNA_expression_total_voom_mmu <- as.data.frame(voom(mRNA_expression_total_mmu)$E)
mRNA_expression_mmu <- merge(mRNA_expression_total_voom_mmu, mRNA_expression_small_voom_mmu, by = 0)
row.names(mRNA_expression_mmu) <- mRNA_expression_mmu$Row.names
mRNA_expression_mmu <- mRNA_expression_mmu[,-1]

# Consider only Coding Sequences
CDS_names_mmu <- read.table("datasets/gencode.vM25.protein_coding.txt.gz")$V1
mRNA_expression_CDS_mmu <- mRNA_expression_mmu[CDS_names_mmu,]

# Fits (linear models)
brain.cds.lm.mmu <- lm(mRNA_expression_CDS_mmu$brain_1_total ~ mRNA_expression_CDS_mmu$brain_small)
cerebellum.cds.lm.mmu <- lm(mRNA_expression_CDS_mmu$cerebellum_total ~ mRNA_expression_CDS_mmu$cerebellum_small)
heart.cds.lm.mmu <- lm(mRNA_expression_CDS_mmu$heart_total ~ mRNA_expression_CDS_mmu$heart_small)
kidney.cds.lm.mmu <- lm(mRNA_expression_CDS_mmu$kidney_total ~ mRNA_expression_CDS_mmu$kidney_small)
testis.cds.lm.mmu <- lm(mRNA_expression_CDS_mmu$testis_total ~ mRNA_expression_CDS_mmu$testis_small)


# Output table 1.
output_table_lm <- data.frame(species = character(), 
                              tissue = character(),
                              AdjustedR2 = numeric(),
                              Pvalue = numeric())

# Compute lms
brain.cds.lm.mmu.summary <- summary(brain.cds.lm.mmu)
cerebellum.cds.lm.mmu.summary <- summary(cerebellum.cds.lm.mmu)
heart.cds.lm.mmu.summary <- summary(heart.cds.lm.mmu)
kidney.cds.lm.mmu.summary <- summary(kidney.cds.lm.mmu)
testis.cds.lm.mmu.summary <- summary(testis.cds.lm.mmu)

# add new rows to output table
output_table_lm <- rbind(output_table_lm, 
                         data.frame(species = "mouse",
                                    tissue = "brain",
                                    AdjustedR2 = brain.cds.lm.mmu.summary$adj.r.squared, 
                                    Pvalue = pf(brain.cds.lm.mmu.summary$fstatistic, 1, brain.cds.lm.mmu.summary$df[2], lower.tail = FALSE)[1])
                         )
output_table_lm <- rbind(output_table_lm, 
                         data.frame(species = "mouse",
                                    tissue = "cerebellum",
                                    AdjustedR2 = cerebellum.cds.lm.mmu.summary$adj.r.squared, 
                                    Pvalue = pf(cerebellum.cds.lm.mmu.summary$fstatistic, 1, cerebellum.cds.lm.mmu.summary$df[2], lower.tail = FALSE)[1])
)
output_table_lm <- rbind(output_table_lm, 
                         data.frame(species = "mouse",
                                    tissue = "heart",
                                    AdjustedR2 = heart.cds.lm.mmu.summary$adj.r.squared, 
                                    Pvalue = pf(heart.cds.lm.mmu.summary$fstatistic, 1, heart.cds.lm.mmu.summary$df[2], lower.tail = FALSE)[1])
)
output_table_lm <- rbind(output_table_lm, 
                         data.frame(species = "mouse",
                                    tissue = "kidney",
                                    AdjustedR2 = kidney.cds.lm.mmu.summary$adj.r.squared, 
                                    Pvalue = pf(kidney.cds.lm.mmu.summary$fstatistic, 1, kidney.cds.lm.mmu.summary$df[2], lower.tail = FALSE)[1])
)
output_table_lm <- rbind(output_table_lm, 
                         data.frame(species = "mouse",
                                    tissue = "testis",
                                    AdjustedR2 = testis.cds.lm.mmu.summary$adj.r.squared, 
                                    Pvalue = pf(testis.cds.lm.mmu.summary$fstatistic, 1, testis.cds.lm.mmu.summary$df[2], lower.tail = FALSE)[1])
)


## Linear models for Chicken datasets
#
# Small RNA mapped to all genes
mRNA_expression_small_gga <- read.table("datasets/Meunier_mRNA_gga.tab.gz", header = TRUE, row.names = 1)
mRNA_expression_small_gga <- mRNA_expression_small_gga[,-c(1:5)]
colnames(mRNA_expression_small_gga) <- c("brain_small", "cerebellum_small", "heart_small", "kidney_small", "testis_small")

# Total RNA mapped to all genes
mRNA_expression_total_gga <- read.table("datasets/Brawand_mRNA_gga.tab.gz", header = TRUE, row.names = 1)
mRNA_expression_total_gga <- mRNA_expression_total_gga[,-c(1:5)]
colnames(mRNA_expression_total_gga) <- c("brain_total", "cerebellum_total", "heart_total", "liver_total", "testis_1_total", "testis_2_total", "testis_3_total")

# Voom transformation
mRNA_expression_small_voom_gga <- as.data.frame(voom(mRNA_expression_small_gga)$E)
mRNA_expression_total_voom_gga <- as.data.frame(voom(mRNA_expression_total_gga)$E)
mRNA_expression_gga <- merge(mRNA_expression_total_voom_gga, mRNA_expression_small_voom_gga, by = 0)
row.names(mRNA_expression_gga) <- mRNA_expression_gga$Row.names
mRNA_expression_gga <- mRNA_expression_gga[,-1]

# Consider only Coding Sequences
CDS_names_gga <- read.table("datasets/Gga_106.protein_coding.txt.gz")$V1
mRNA_expression_CDS_gga <- mRNA_expression_gga[CDS_names_gga,]

# Fits (linear models)
brain.cds.lm.gga <- lm(mRNA_expression_CDS_gga$brain_total ~ mRNA_expression_CDS_gga$brain_small)
cerebellum.cds.lm.gga <- lm(mRNA_expression_CDS_gga$cerebellum_total ~ mRNA_expression_CDS_gga$cerebellum_small)
heart.cds.lm.gga <- lm(mRNA_expression_CDS_gga$heart_total ~ mRNA_expression_CDS_gga$heart_small)
testis.cds.lm.gga <- lm(mRNA_expression_CDS_gga$testis_1_total ~ mRNA_expression_CDS_gga$testis_small)

brain.cds.lm.gga.summary <- summary(brain.cds.lm.gga)
cerebellum.cds.lm.gga.summary <- summary(cerebellum.cds.lm.gga)
heart.cds.lm.gga.summary <- summary(heart.cds.lm.gga)
testis.cds.lm.gga.summary <- summary(testis.cds.lm.gga)

# add new rows to output table
output_table_lm <- rbind(output_table_lm, 
                         data.frame(species = "chicken",
                                    tissue = "brain",
                                    AdjustedR2 = brain.cds.lm.gga.summary$adj.r.squared, 
                                    Pvalue = pf(brain.cds.lm.gga.summary$fstatistic, 1, brain.cds.lm.gga.summary$df[2], lower.tail = FALSE)[1])
)
output_table_lm <- rbind(output_table_lm, 
                         data.frame(species = "chicken",
                                    tissue = "cerebellum",
                                    AdjustedR2 = cerebellum.cds.lm.gga.summary$adj.r.squared, 
                                    Pvalue = pf(cerebellum.cds.lm.gga.summary$fstatistic, 1, cerebellum.cds.lm.gga.summary$df[2], lower.tail = FALSE)[1])
)
output_table_lm <- rbind(output_table_lm, 
                         data.frame(species = "chicken",
                                    tissue = "heart",
                                    AdjustedR2 = heart.cds.lm.gga.summary$adj.r.squared, 
                                    Pvalue = pf(heart.cds.lm.gga.summary$fstatistic, 1, heart.cds.lm.gga.summary$df[2], lower.tail = FALSE)[1])
)
output_table_lm <- rbind(output_table_lm, 
                         data.frame(species = "chicken",
                                    tissue = "testis",
                                    AdjustedR2 = testis.cds.lm.gga.summary$adj.r.squared, 
                                    Pvalue = pf(testis.cds.lm.gga.summary$fstatistic, 1, testis.cds.lm.gga.summary$df[2], lower.tail = FALSE)[1])
)


write.table(output_table_lm, file = "plots/table_R2_mmugga.tab", sep = "\t", quote = FALSE, row.names = FALSE)


# Cellular Component annotation top 10%

# Brain
hsa_brain_top10p <- mRNA_expression_CDS |>
  mutate(rank = percent_rank(brain_small)) |>
  filter(rank >= 0.9) |>
  row.names() 
# Run     WebGestaltR for Human Cell Landscape enrichment
CC_brain <- WebGestaltR(
  organism = "hsapiens",                    # Human species
  enrichMethod = "ORA",                     # Over-Representation Analysis
  interestGene = hsa_brain_top10p,                # Input gene list
  interestGeneType = "genesymbol",          # Gene type (e.g., genesymbol, entrezgene, etc.)
  referenceGene = row.names(mRNA_expression_CDS),         # Background gene list
  referenceGeneType = "genesymbol",         # Background gene type
  fdrMethod = "BH",                         # Multiple testing correction (Benjamini-Hochberg)
  isOutput = FALSE,                          # Save results to file
  enrichDatabase = "geneontology_Cellular_Component"           # Database: Human Cell Landscape
)
hsa_brain_top10p_CC <- CC_brain |>
  arrange(FDR) |>
  slice_head(n = 3) |>
  dplyr::select(description, enrichmentRatio, FDR)

# Cerebellum
hsa_cerebellum_top10p <- mRNA_expression_CDS |>
  mutate(rank = percent_rank(cerebellum_small)) |>
  filter(rank >= 0.9) |>
  row.names() 
# Run     WebGestaltR for Human Cell Landscape enrichment
CC_cerebellum <- WebGestaltR(
  organism = "hsapiens",                    # Human species
  enrichMethod = "ORA",                     # Over-Representation Analysis
  interestGene = hsa_cerebellum_top10p,                # Input gene list
  interestGeneType = "genesymbol",          # Gene type (e.g., genesymbol, entrezgene, etc.)
  referenceGene = row.names(mRNA_expression_CDS),         # Background gene list
  referenceGeneType = "genesymbol",         # Background gene type
  fdrMethod = "BH",                         # Multiple testing correction (Benjamini-Hochberg)
  isOutput = FALSE,                          # Save results to file
  enrichDatabase = "geneontology_Cellular_Component"           # Database: Human Cell Landscape
)
hsa_cerebellum_top10p_CC <- CC_cerebellum |>
  arrange(FDR) |>
  slice_head(n = 3) |>
  dplyr::select(description, enrichmentRatio, FDR)

# Heart
hsa_heart_top10p <- mRNA_expression_CDS |>
  mutate(rank = percent_rank(heart_small)) |>
  filter(rank >= 0.9) |>
  row.names() 
# Run     WebGestaltR for Human Cell Landscape enrichment
CC_heart <- WebGestaltR(
  organism = "hsapiens",                    # Human species
  enrichMethod = "ORA",                     # Over-Representation Analysis
  interestGene = hsa_heart_top10p,                # Input gene list
  interestGeneType = "genesymbol",          # Gene type (e.g., genesymbol, entrezgene, etc.)
  referenceGene = row.names(mRNA_expression_CDS),         # Background gene list
  referenceGeneType = "genesymbol",         # Background gene type
  fdrMethod = "BH",                         # Multiple testing correction (Benjamini-Hochberg)
  isOutput = FALSE,                          # Save results to file
  enrichDatabase = "geneontology_Cellular_Component"           # Database: Human Cell Landscape
)
hsa_heart_top10p_CC <- CC_heart |>
  arrange(FDR) |>
  slice_head(n = 3) |>
  dplyr::select(description, enrichmentRatio, FDR)


# Kidney
hsa_kidney_top10p <- mRNA_expression_CDS |>
  mutate(rank = percent_rank(kidney_small)) |>
  filter(rank >= 0.9) |>
  row.names() 
# Run     WebGestaltR for Human Cell Landscape enrichment
CC_kidney <- WebGestaltR(
  organism = "hsapiens",                    # Human species
  enrichMethod = "ORA",                     # Over-Representation Analysis
  interestGene = hsa_kidney_top10p,                # Input gene list
  interestGeneType = "genesymbol",          # Gene type (e.g., genesymbol, entrezgene, etc.)
  referenceGene = row.names(mRNA_expression_CDS),         # Background gene list
  referenceGeneType = "genesymbol",         # Background gene type
  fdrMethod = "BH",                         # Multiple testing correction (Benjamini-Hochberg)
  isOutput = FALSE,                          # Save results to file
  enrichDatabase = "geneontology_Cellular_Component"           # Database: Human Cell Landscape
)
hsa_kidney_top10p_CC <- CC_kidney |>
  arrange(FDR) |>
  slice_head(n = 3) |>
  dplyr::select(description, enrichmentRatio, FDR)


# Association between sequencing depth and R-squared
count_r_squared <- read_table("datasets/read_counts_and_r-squared.tab")
cor.test(count_r_squared$count/count_r_squared$Genome.Mb, count_r_squared$Rsquared)
cor.test(count_r_squared$Count.over.25nt/count_r_squared$Genome.Mb, count_r_squared$Rsquared)
count_r_squared_notestis <- count_r_squared |>
  filter(tissue != "testis")
cor.test(count_r_squared_notestis$Count.over.25nt/count_r_squared_notestis$Genome.Mb, count_r_squared_notestis$Rsquared)
# 0.6629908, p 0.02618
cor.test(count_r_squared_notestis$count/count_r_squared_notestis$Genome.Mb, count_r_squared_notestis$Rsquared)


# KEGG annotation top 10%

# Brain
# Run     WebGestaltR for Human Cell Landscape enrichment
KEGG_brain <- WebGestaltR(
  organism = "hsapiens",                    # Human species
  enrichMethod = "ORA",                     # Over-Representation Analysis
  interestGene = hsa_brain_top10p,                # Input gene list
  interestGeneType = "genesymbol",          # Gene type (e.g., genesymbol, entrezgene, etc.)
  referenceGene = row.names(mRNA_expression_CDS),         # Background gene list
  referenceGeneType = "genesymbol",         # Background gene type
  fdrMethod = "BH",                         # Multiple testing correction (Benjamini-Hochberg)
  isOutput = FALSE,                          # Save results to file
  enrichDatabase = "pathway_KEGG"           # Database: Human Cell Landscape
)
hsa_brain_top10p_KEGG <- KEGG_brain |>
  arrange(FDR) |>
  slice_head(n = 3) |>
  dplyr::select(description, enrichmentRatio, FDR)

# Cerebellum
# Run     WebGestaltR for Human Cell Landscape enrichment
KEGG_cerebellum <- WebGestaltR(
  organism = "hsapiens",                    # Human species
  enrichMethod = "ORA",                     # Over-Representation Analysis
  interestGene = hsa_cerebellum_top10p,                # Input gene list
  interestGeneType = "genesymbol",          # Gene type (e.g., genesymbol, entrezgene, etc.)
  referenceGene = row.names(mRNA_expression_CDS),         # Background gene list
  referenceGeneType = "genesymbol",         # Background gene type
  fdrMethod = "BH",                         # Multiple testing correction (Benjamini-Hochberg)
  isOutput = FALSE,                          # Save results to file
  enrichDatabase = "pathway_KEGG"           # Database: Human Cell Landscape
)
hsa_cerebellum_top10p_KEGG <- KEGG_cerebellum |>
  arrange(FDR) |>
  slice_head(n = 3) |>
  dplyr::select(description, enrichmentRatio, FDR)

# Heart
# Run     WebGestaltR for Human Cell Landscape enrichment
KEGG_heart <- WebGestaltR(
  organism = "hsapiens",                    # Human species
  enrichMethod = "ORA",                     # Over-Representation Analysis
  interestGene = hsa_heart_top10p,                # Input gene list
  interestGeneType = "genesymbol",          # Gene type (e.g., genesymbol, entrezgene, etc.)
  referenceGene = row.names(mRNA_expression_CDS),         # Background gene list
  referenceGeneType = "genesymbol",         # Background gene type
  fdrMethod = "BH",                         # Multiple testing correction (Benjamini-Hochberg)
  isOutput = FALSE,                          # Save results to file
  enrichDatabase = "pathway_KEGG"           # Database: Human Cell Landscape
)
hsa_heart_top10p_KEGG <- KEGG_heart |>
  arrange(FDR) |>
  slice_head(n = 3) |>
  dplyr::select(description, enrichmentRatio, FDR)


# Kidney
# Run     WebGestaltR for Human Cell Landscape enrichment
KEGG_kidney <- WebGestaltR(
  organism = "hsapiens",                    # Human species
  enrichMethod = "ORA",                     # Over-Representation Analysis
  interestGene = hsa_kidney_top10p,                # Input gene list
  interestGeneType = "genesymbol",          # Gene type (e.g., genesymbol, entrezgene, etc.)
  referenceGene = row.names(mRNA_expression_CDS),         # Background gene list
  referenceGeneType = "genesymbol",         # Background gene type
  fdrMethod = "BH",                         # Multiple testing correction (Benjamini-Hochberg)
  isOutput = FALSE,                          # Save results to file
  enrichDatabase = "pathway_KEGG"           # Database: Human Cell Landscape
)
hsa_kidney_top10p_KEGG <- KEGG_kidney |>
  arrange(FDR) |>
  slice_head(n = 3) |>
  dplyr::select(description, enrichmentRatio, FDR)


# Output table annotation
output_annotation_table <- rbind(hsa_brain_top10p_KEGG,
                                 hsa_cerebellum_top10p_KEGG,
                                 hsa_heart_top10p_KEGG,
                                 hsa_kidney_top10p_KEGG,
                                 hsa_brain_top10p_CC,
                                 hsa_cerebellum_top10p_CC,
                                 hsa_heart_top10p_CC,
                                 hsa_kidney_top10p_CC) |>
  mutate(tissue = rep(c(rep("brain", 3), 
                    rep("cerebellum",3), 
                    rep("heart",3), 
                    rep("kidney",3)),2)) |>
  mutate(database = c(rep("KEGG", 12), rep("Cellular Component (GO)", 12)))

write.table(output_annotation_table, file = "plots/table_annotation_10p.tab", sep = "\t", quote = FALSE, row.names = FALSE)



## Relative coverage for transcript section 5UTR, CDS, 3UTR.

# Load files with CDS/UTR read mappings
untar("datasets/hsa_mapping_UTR_CDS.tar.gz",
      exdir = tempdir())
Brawand_CDS <- read.table(file.path(tempdir(), 
                                    "GSE30352_Brawand", 
                                    "Brawand_mRNA_CDS.tab"), 
                          header = TRUE, 
                          sep = "\t")[,-c(8,11)] # Subset 4 datasets from first analysis
Brawand_UTR5 <- read.table(file.path(tempdir(), 
                                     "GSE30352_Brawand", 
                                     "Brawand_mRNA_UTR5.tab"), 
                           header = TRUE,
                           sep = "\t")[,-c(8,11)] # Subset 4 datasets from first analysis
Brawand_UTR3 <- read.table(file.path(tempdir(), 
                                     "GSE30352_Brawand",
                                     "Brawand_mRNA_UTR3.tab"), 
                           header = TRUE,
                           sep = "\t")[,-c(8,11)] # Subset 4 datasets from first analysis
Meunier_CDS <- read.table(file.path(tempdir(), 
                                    "GSE40499_Meunier", 
                                    "Meunier_mRNA_CDS.tab"),
                          header = TRUE, 
                          sep = "\t")[,-11] # Subset 4 datasets from first analysis
Meunier_UTR5 <- read.table(file.path(tempdir(), 
                                     "GSE40499_Meunier",
                                     "Meunier_mRNA_UTR5.tab"), 
                           header = TRUE, 
                           sep = "\t")[,-11] # Subset 4 datasets from first analysis
Meunier_UTR3 <- read.table(file.path(tempdir(),
                                     "GSE40499_Meunier",
                                     "Meunier_mRNA_UTR3.tab"), 
                           header = TRUE,
                           sep = "\t")[,-11] # Subset 4 datasets from first analysis

# Total reads and lengths
Brawand_CDS_total <- Brawand_CDS[,c(6:10)] |> colSums() 
Brawand_UTR5_total <- Brawand_UTR5[,c(6:10)] |> colSums() 
Brawand_UTR3_total <- Brawand_UTR3[,c(6:10)] |> colSums() 
Meunier_CDS_total <- Meunier_CDS[,c(6:10)] |> colSums() 
Meunier_UTR5_total <- Meunier_UTR5[,c(6:10)] |> colSums() 
Meunier_UTR3_total <- Meunier_UTR3[,c(6:10)] |> colSums() 

# BRAWAND; summary tables
Brawand_total <- rbind(Brawand_CDS_total,
                       Brawand_UTR5_total,
                       Brawand_UTR3_total)
# Counts per kilobase
Brawand_count_kb <- t(Brawand_total[,c(2:5)]*1000/Brawand_total[,1])
# log2FC relative to CDS
Brawand_log2_CDS <- log2(Brawand_count_kb/Brawand_count_kb[,1])

# MEUNIER; summary tables
Meunier_total <- rbind(Meunier_CDS_total,
                       Meunier_UTR5_total,
                       Meunier_UTR3_total)
Meunier_count_kb <- t(Meunier_total[,c(2:5)]*1000/Meunier_total[,1])
Meunier_log2_CDS <- log2(Meunier_count_kb/Meunier_count_kb[,1])

# Plots
# Brawand
plot_B_UTR <- tibble(dataset = rep(c(1:4), 2),
                     UTR = c(rep("3' UTR", 4), rep("5' UTR", 4)),
                     log2FC = c(Brawand_log2_CDS[,2], Brawand_log2_CDS[,3])) |>
  ggplot(aes(y = log2FC,
             x = UTR)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("5' UTR", "3' UTR")) +
  geom_point(aes(y = log2FC,
                  x = UTR)) +
  geom_line(aes(group = dataset), alpha = 0.5, col = "grey") +
  geom_hline(yintercept=0) +
  labs(y = "log2 fold-change UTR/CDS",
       x = "") +
  ylim(c(-2,3)) +
  theme_bw()
# Meunier
plot_M_UTR <- tibble(dataset = rep(c(1:4), 2),
                     UTR = c(rep("3' UTR", 4), rep("5' UTR", 4)),
                     log2FC = c(Meunier_log2_CDS[,2], Meunier_log2_CDS[,3])) |>
  ggplot(aes(y = log2FC,
             x = UTR)) +
  scale_x_discrete(limits = c("5' UTR", "3' UTR")) +
  geom_boxplot() +
  geom_point(aes(y = log2FC,
                 x = UTR)) +
  geom_line(aes(group = dataset), alpha = 0.5, col = "grey") +
  geom_hline(yintercept=0) +
  labs(y = "",
       x = "") +
  ylim(c(-2,3)) +
  theme_bw()
# Brain first and cerebellum after are the ones more similar to total RNA!
ggsave((plot_B_UTR + plot_M_UTR) + plot_annotation(tag_levels = 'A'), file = "plots/Figure_UTRs.png")



## Half life analysis
# Load and pre-process half-life file
half_lives <- read_tsv("datasets/Tani_transcript_half_life.tab.gz") |>
  dplyr::rename(
    refseq_mrna = 'RepName',
    length = length,
    rpkm0 = 'RPKM (0h)',
    half_file = 't1/2 (h)'
  ) |>
  dplyr::mutate(
    half_file = as.numeric(if_else(half_file == ">24", "25", half_file))
  ) |>
  tidyr::separate_rows(refseq_mrna, sep = ",") |> # Split multiple accessions into single rows
  dplyr::filter(refseq_mrna != "")
  
# Get gene names for half-life file
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
refseq2name <- getBM(
    attributes = c("refseq_mrna", "hgnc_symbol"),
    filters = "refseq_mrna",
    values = half_lives$refseq_mrna,
    mart = mart
  ) |> as_tibble()
half_lives <- half_lives |>
  left_join(refseq2name) |>
  dplyr::select(-refseq_mrna) |>
  unique() |>
  filter(!is.na(hgnc_symbol)) |>
  distinct(hgnc_symbol, .keep_all = TRUE) # Remove multiple values for the same gene, there are only a few.

# Join expression and half live
mRNA_expression_CDS_tibble <- mRNA_expression_CDS |> rownames_to_column(var = "hgnc_symbol") |> as_tibble()
expression_halflives <- mRNA_expression_CDS_tibble |>
  dplyr::inner_join(half_lives)


# Comparisons
anova(lm(brain_small ~ brain_total * half_file, data = expression_halflives))
summary(lm(brain_small ~ brain_total * half_file, data = expression_halflives))
anova(lm(cerebellum_small ~ cerebellum_1_total * half_file, data = expression_halflives))
summary(lm(cerebellum_small ~ cerebellum_1_total * half_file, data = expression_halflives))
anova(lm(heart_small ~ heart_1_total * half_file, data = expression_halflives))
summary(lm(heart_small ~ heart_1_total * half_file, data = expression_halflives))
anova(lm(kidney_small ~ kidney_total * half_file, data = expression_halflives))
summary(lm(kidney_small ~ kidney_total * half_file, data = expression_halflives))

q()
