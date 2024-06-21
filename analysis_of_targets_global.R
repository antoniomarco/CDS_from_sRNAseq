library(ggplot2)
library(patchwork)
# Bioconductor
library(edgeR)


# TARGET ANALYSIS

# MirTarBase
targets_MTB <- read.table("datasets/mirtarbase_hsa_9.tab.gz")
# Small RNA expression mapped to all genes
mRNA_expression <- read.table("datasets/Meunier_mRNA.tab.gz", header = TRUE, row.names = 1)
mRNA_expression <- mRNA_expression[,-c(1:5)]
# Small RNAs mapped to microRNAs
microRNA_expression <- read.table("datasets/Meunier_microRNA.tab.gz", header = TRUE, row.names = 1)
microRNA_expression <- microRNA_expression[,-c(1:5)]
# Remove genes/miRNAs with zero reads in all samples
mRNA_expression <- mRNA_expression[which(rowSums(mRNA_expression)!=0),]
microRNA_expression <- microRNA_expression[which(rowSums(microRNA_expression)!=0),]
# mRNAs and microRNA with target information, subset
# remove those inteaction reported for human annotated microRNAs
microRNAs_MTB <- unique(targets_MTB[,1][grepl("^hsa-",targets_MTB[,1])])
mRNA_MTB <- unique(targets_MTB[,2][grepl("^hsa-",targets_MTB[,1])])
# Only mRNs/microRNAs in all sets
microRNA_expression <- microRNA_expression[names(which(table(c(row.names(microRNA_expression),microRNAs_MTB))==2)),]
mRNA_expression <- mRNA_expression[names(which(table(c(row.names(mRNA_expression),mRNA_MTB))==2)),]

# Normalization, expression levels. VOOM TRANSFORATION
mRNA_expression_voom <- voom(mRNA_expression)$E
microRNA_expression_voom <- voom(microRNA_expression)$E

# # Uncomment to run from scratch, otherwise, using pre-computed correlations
# # Correlations table
# correlations_table <- matrix(nrow=nrow(mRNA_expression_voom), ncol = 0)
# total <- length(rownames(microRNA_expression_voom))
# for(i in 1:total)
# {
#   print(paste("MicroRNA ", i, " out of ", total, sep = ""))
#   mycor <- function(input){return(cor(as.numeric(input),as.numeric(microRNA_expression_voom[i,]), method = "spearman"))}
#   correlations_table <- cbind(correlations_table, apply(mRNA_expression_voom,1,mycor))
# }
# correlations.df <- as.data.frame(correlations_table)
# colnames(correlations.df) <- rownames(microRNA_expression_voom)
# save(correlations.df, file = "datasets/correlations.Rdata")
load("datasets/correlations.Rdata")

# Filter targets for only those in the expression dataset
targets_MTB_clean <- targets_MTB[targets_MTB$V1 %in% colnames(correlations.df),]
targets_MTB_clean <- targets_MTB_clean[targets_MTB_clean$V2 %in% rownames(correlations.df),]

# Correlation between target and non-target pairs
corr_targets <- vector()
correlations.df_notargets <- correlations.df
for(i in 1:nrow(targets_MTB_clean)){
	if((i %% 1000) == 0){
		cat(paste(i, " out of ", nrow(targets_MTB_clean), " (", round((i/nrow(targets_MTB_clean)*100), digits = 0), "%)\r", sep = ""))
	}
	corr_targets <- c(corr_targets, correlations.df[targets_MTB_clean[i,2],targets_MTB_clean[i,1]])
	correlations.df_notargets[targets_MTB_clean[i,2],targets_MTB_clean[i,1]] <- 1000 # 1000 can't result from cor function
}

correlations.df_notargets.unlisted <- unlist(correlations.df_notargets)
corr_non_targets <- as.vector(correlations.df_notargets.unlisted)
corr_non_targets <- corr_non_targets[corr_non_targets < 1000] # Remove values from target interactions
corr_comparison <- data.frame(corr = c(corr_targets, corr_non_targets), targ = c(rep("target", length(corr_targets)), rep("non-target", length(corr_non_targets))))

# MW TEST
wilcox.test(corr_targets, corr_non_targets, alternative = "less")
# W = 2.1545e+12, p-value < 2.2e-16

# remove NAs
corr_targets <- corr_targets[!is.na(corr_targets)]
corr_non_targets <- corr_non_targets[!is.na(corr_non_targets)]

# Smoothing function
smoothing_exp_correl <- function (correl = seq(-1, 1, 0.1), w = 0.1, s = 0.1) {
  size <- length(correl)
  mids <- seq(-1+(w/2), 1-(w/2), s)
  counts <- rep(0,length(mids))
  for(i in 1:length(correl)){
    if(correl[i] != 1){
      which_bin <- which((correl[i] >= mids-(w/2)) & (correl[i] < mids+(w/2)))
      counts[which_bin] <- counts[which_bin] + 1
    } else {
      counts[length(counts)] <- counts[length(counts)] + 1
    }
  }
  return(counts/size)
}

# Smooth
corr_targets_smooth_s02_w02 <- smoothing_exp_correl(corr_targets, s = 0.02, w = 0.2)
corr_non_targets_smooth_s02_w02 <- smoothing_exp_correl(corr_non_targets, s = 0.02, w = 0.2)
MTB_corrs_s02_w02 <- log2(corr_targets_smooth_s02_w02/corr_non_targets_smooth_s02_w02)
data.corrs.s02.w02 <- data.frame(mid = seq(-0.9, 0.9, 0.02), corr = MTB_corrs_s02_w02)

# Plot
plot_MTB <- ggplot(data = data.corrs.s02.w02, aes(x = mid, y = corr, col = 'darkred'))+ 
  geom_line() +
  geom_hline(yintercept=0, linetype="solid", color = "black") +
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  xlab("Correlation (Pearson's coefficient)") +
  ylab("log2(target/non-target)") +
  theme_bw() +
  theme(legend.position="none")

plot_MTB
ggsave(plot_MTB , file = "plots/Figure_3.png")


q()
