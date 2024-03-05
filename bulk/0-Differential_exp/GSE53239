setwd("/home/leandro/Documents/IC/GSE53239_DLPFC/")

library(GEOquery)
library(tidyverse)
library(dplyr)
library(DESeq2)

gset <- getGEO("GSE53239")

meta1 <- gset$`GSE53239-GPL15433_series_matrix.txt.gz`@phenoData@data

meta2 <- gset$`GSE53239-GPL9115_series_matrix.txt.gz`@phenoData@data

metadata <- rbind(meta1, meta2)

exp_set1 <- read.delim("GSE53239_HGB1-10_geneLevel_HTSeqCountData.txt")
exp_set2 <- read.delim("GSE53239_HGB11-22_geneLevel_HTSeqCountData.txt")
exp_set2 <- exp_set2[,-c(1)]
exp_set <- cbind(exp_set1, exp_set2)

rm(exp_set1, exp_set2, meta1, meta2, gset)

###################### Collapsing genes by  sum ###################### 
# Checking for duplicated genes
n_occur <- data.frame(table(exp_set$GeneName))
n_occur[n_occur$Freq > 1,]

# Sum duplicated genes to normalize the expression set
exp_set = aggregate(. ~ GeneName, data=exp_set, FUN=sum)

# define genes as row names and delete the gene_symbol column
row.names(exp_set) = exp_set$GeneName; exp_set$GeneName = NULL

# Clean environment
rm(ensembl_genes, ensembl, gene_symbols, n_occur)

##################### BD vs CTRL #################
control_label <- "disease state: Control"
classes <- "disease state: BD"

# Subset metadata to include only samples of the chosen class and controls
metadata <- metadata %>%
  filter(characteristics_ch1.1 %in% c("disease state: Control",
                                      "disease state: BD"))

# Transform class column into a factor with the levels in the correct order:
# FIRST CONTROL, SECOND TREATED.
metadata$characteristics_ch1.1 <- factor(metadata$characteristics_ch1.1, levels = 
                                       c(control_label, classes))

# Subset count matrix to keep samples of the chosen class and controls
colnames <- str_split(metadata$title, pattern=" ",simplify = T)
rownames(metadata) <- colnames[,2]
metadata <- metadata[!rownames(metadata) == "HGB_5",]

exp_set <- exp_set[, rownames(metadata)]


cnt <- as.matrix(exp_set)

meta <- metadata

# create dds object
meta$platform_id <- as.factor(meta$platform_id)
dds <- DESeqDataSetFromMatrix(countData = cnt,
                                   colData = meta,
                                   design = ~ platform_id + characteristics_ch1.1,
                                   tidy = F)


# Plot PCA based on diagnostic
vsd <- vst(dds, blind=FALSE)

norm <- as.data.frame(assay(vsd))


###################### Using MDP function on data normalized with DESeq ####
library(mdp)

metatable <- metadata %>% rownames_to_column(var= "Sample")
rownames(metatable) <- metatable$Sample

samplesinfo <- metatable[, c("Sample","characteristics_ch1.1")]
colnames(samplesinfo) <- c("Sample", 
                           "Class")

norm <- norm[, rownames(samplesinfo)]


mdp_deseq <- mdp(data = norm,
                 pdata = samplesinfo, 
                 control_lab = "disease state: Control")

###################### Plot boxplot for all the samples ####
library(ggplot2)

boxplot_deseq <- mdp_deseq$sample_scores$allgenes
boxplot_deseq <- boxplot_deseq[boxplot_deseq$Class == "disease state: Control",]
boxplot_deseq %>% 
  ggplot(aes(x = Class, y = zscore_class, fill = Class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "black") +
  labs(title = "All genes - GSE53239") +
  theme_bw() +
  theme(legend.position = "none")
boxplot_deseq <- mdp_deseq$sample_scores$perturbedgenes
boxplot_deseq <- boxplot_deseq[boxplot_deseq$Class == "disease state: Control",]
boxplot_deseq %>% 
  ggplot(aes(x = Class, y = zscore_class, fill = Class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "black") +
  labs(title = "Perturbed genes - GSE53239") +
  theme_bw() + 
  theme(legend.position = "none")
## Get outlier samples

outliers_deseq <- c()

for(i in which(mdp_deseq$sample_scores$allgenes$outlier == 1)){
  outliers_deseq <- c(outliers_deseq, mdp_deseq$sample_scores$allgenes$Sample[i])
}
print(paste0("Number of outlier samples considering all genes + DESeq normalization: ", length(outliers_deseq)))

exp_set <- exp_set[,!colnames(exp_set) %in% outliers_deseq]
norm <- norm[,colnames(exp_set)]
metatable <- metatable[colnames(exp_set),]

outliers_deseq <- c()

for(i in which(mdp_deseq$sample_scores$perturbedgenes$outlier == 1)){
  outliers_deseq <- c(outliers_deseq, mdp_deseq$sample_scores$perturbedgenes$Sample[i])
}
print(paste0("Number of outlier samples considering perturbed genes + DESeq normalization: ", length(outliers_deseq)))

exp_set <- exp_set[,!colnames(exp_set) %in% outliers_deseq]
norm <- norm[,colnames(exp_set)]
metatable <- metatable[colnames(exp_set),]

write.table(norm,
            file = 'deconvolution/norm.tsv',
            sep = '\t',
            row.names = TRUE)
write.table(metatable,
            file = 'metatable.tsv',
            sep = '\t',
            row.names = TRUE)


# 2 in total

########################

# Define the targets/conditions (all classes but control/healthy)
classes <- "disease state: BD"

# Define the control group label
control_label <- "disease state: Control"

# Subset metadata to include only samples of the chosen class and controls
meta <- metatable %>%
  filter(characteristics_ch1.1 %in% c("disease state: BD",
                                      "disease state: Control"))

# Transform class column into a factor with the levels in the correct order:
# FIRST CONTROL, SECOND TREATED.
meta$class <- factor(meta$characteristics_ch1.1, levels = c(control_label, classes))

colnames <- str_split(meta$title, pattern=" ",simplify = T)
rownames(meta) <- colnames[,2]
meta <- meta[!rownames(meta) == "HGB_5",]
# Subset count matrix to keep samples of the chosen class and controls
cnt <- as.matrix(exp_set[,colnames(exp_set) %in% rownames(meta)])

# create dds object. There is no BD patient from cohort B
meta$platform_id <- as.factor(meta$platform_id)
dds <- DESeqDataSetFromMatrix(countData = cnt,
                                   colData = meta,
                                   design = ~ platform_id + characteristics_ch1.1,
                                   tidy = F)
cpm <- fpm(dds,robust = F)

write.table(cpm,
            file = 'deconvolution/cpm.tsv',
            sep = '\t',
            row.names = TRUE)

# Remove genes that have zero counts in all samples (i.e. not expressed genes)
keep <- rowSums(counts(dds)) > 1
dds_trat <- dds[keep,]

# Run DE analysis with default DESeq2 settings
dds_trat <- DESeq(dds)

# extract results defining the order of comparison
# must be: constrast=c("CLASS COLUMN IN META","DISEASE/TREATMENT GROUP","CONTROL)
# This will give log2FoldChange values that correspond to disease vs control
res <- as.data.frame(results(dds_trat, contrast = c("characteristics_ch1.1",
                                                    classes,
                                                    control_label),cooksCutoff = F))
res <- res[(res$log2FoldChange < -1 | res$log2FoldChange > 1) & res$padj < 0.05
           & !is.na(res$padj),]
write.table(res,
            file = 'DEGs_BD_vs_Control.tsv',
            sep = '\t')
