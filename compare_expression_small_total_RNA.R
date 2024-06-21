library(ggplot2)
library(patchwork)
# Bioconductor
library(biomaRt)
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

# Fits (linear models)
brain.lm <- lm(mRNA_expression$brain_total ~ mRNA_expression$brain_small)
cerebellum.lm <- lm(mRNA_expression$cerebellum_1_total ~ mRNA_expression$cerebellum_small)
heart.lm <- lm(mRNA_expression$heart_1_total ~ mRNA_expression$heart_small)
kidney.lm <- lm(mRNA_expression$kidney_total ~ mRNA_expression$kidney_small)

# Plots 
brain_plot <- ggplot(data = mRNA_expression,
       aes(x = brain_small,
           y = brain_total)) +
       ggtitle("Brain") +
       xlab("smallRNA-seq Expression level") +
       ylab("RNA-seq Expression level") +
       geom_point(alpha=0.4, col = "darkgrey") +
       geom_abline(slope = coef(brain.lm)[["mRNA_expression$brain_small"]], 
                   intercept = coef(brain.lm)[["(Intercept)"]],
                   col = 'darkred') +
       annotate(geom = 'text', x = 13, y = 10, size = 3, label = paste("R^2 ==", round(summary(brain.lm)$adj.r.squared, 2)), parse = TRUE) +
       theme_light()  +
       theme(plot.title = element_text(hjust = 0.5, size = 11)) +
       theme(axis.title = element_text(size = 9)) 

cerebellum_plot <- ggplot(data = mRNA_expression,
       aes(x = cerebellum_small,
           y = cerebellum_1_total)) +
       ggtitle("Cerebellum") +
       xlab("smallRNA-seq Expression level") +
       ylab("RNA-seq Expression level") +
       geom_point(alpha=0.4, col = "darkgrey") +
       geom_abline(slope = coef(cerebellum.lm)[["mRNA_expression$cerebellum_small"]], 
                   intercept = coef(cerebellum.lm)[["(Intercept)"]],
                   col = 'darkred') +
      annotate(geom = 'text', x = 15, y = 10, size = 3, label = paste("R^2 ==", round(summary(cerebellum.lm)$adj.r.squared, 2)), parse = TRUE) +
      theme_light()  +
      theme(plot.title = element_text(hjust = 0.5, size = 11)) +
      theme(axis.title = element_text(size = 9)) 

heart_plot <- ggplot(data = mRNA_expression,
       aes(x = heart_small,
           y = heart_1_total)) +
       ggtitle("Heart") +
       xlab("smallRNA-seq Expression level") +
       ylab("RNA-seq Expression level") +
       geom_point(alpha=0.4, col = "darkgrey") +
       geom_abline(slope = coef(heart.lm)[["mRNA_expression$heart_small"]], 
                   intercept = coef(heart.lm)[["(Intercept)"]],
                   col = 'darkred') +
       annotate(geom = 'text', x = 13, y = 10, size = 3, label = paste("R^2 ==", round(summary(heart.lm)$adj.r.squared, 2)), parse = TRUE) +
       theme_light()  +
       theme(plot.title = element_text(hjust = 0.5, size = 11)) +
       theme(axis.title = element_text(size = 9)) 

kidney_plot <- ggplot(data = mRNA_expression,
       aes(x = kidney_small,
           y = kidney_total)) +
       ggtitle("Kidney") +
       xlab("smallRNA-seq Expression level") +
       ylab("RNA-seq Expression level") +
       geom_point(alpha=0.4, col = "darkgrey") +
       geom_abline(slope = coef(kidney.lm)[["mRNA_expression$kidney_small"]], 
                   intercept = coef(kidney.lm)[["(Intercept)"]],
                   col = 'darkred') +
       annotate(geom = 'text', x = 15, y = 10, size = 3, label = paste("R^2 ==", round(summary(kidney.lm)$adj.r.squared, 2)), parse = TRUE) +
       theme_light() +
       theme(plot.title = element_text(hjust = 0.5, size = 11)) +
       theme(axis.title = element_text(size = 9)) 

(brain_plot | cerebellum_plot)/(heart_plot | kidney_plot) + plot_annotation(tag_levels = "A")
ggsave((brain_plot | cerebellum_plot)/(heart_plot | kidney_plot) + plot_annotation(tag_levels = "A"), file = "plots/Figure_1.png")


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

(brain_CDS_plot | cerebellum_CDS_plot)/(heart_CDS_plot | kidney_CDS_plot) + plot_annotation(tag_levels = "A")
ggsave((brain_CDS_plot | cerebellum_CDS_plot)/(heart_CDS_plot | kidney_CDS_plot) + plot_annotation(tag_levels = "A"), file = "plots/Figure_2.png")

q()
