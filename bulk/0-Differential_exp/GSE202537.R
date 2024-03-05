setwd("/home/leandro/Documents/IC/GSE202537_basal/")

library(GEOquery)
library(tidyverse)
library(qdapRegex)

gset <- getGEO("GSE202537")

exp_set <- read.delim("GSE202537_Read_counts_NAc_Caudate_Putamen_psychosis_mcontrols.csv",
                      sep = ",")
exp_set <- exp_set %>% remove_rownames() %>% column_to_rownames(var="X")

metatable <- gset$GSE202537_series_matrix.txt.gz@phenoData@data

rm(gset); gc()

samples <- ex_between(metatable$title, "[", "]")

metatable$sample <- paste0("X", samples)
rownames(metatable) <- metatable$sample

###################### Get BioMart ###################### 
library(biomaRt)

# Load Database
ensembl = useMart(biomart="ensembl",
                  dataset="hsapiens_gene_ensembl")

# Ensembl gene IDs
ensembl_genes = rownames(exp_set)
ensembl_genes = gsub("\\..*","", ensembl_genes)
exp_set$gene_id = ensembl_genes

# Dataframe of attributes if you want to check
# at = ensembl@attributes

# Get gene_symbols from Ensembl(BioMart)
gene_symbols = getBM(attributes = c('hgnc_symbol',
                                    'ensembl_gene_id'),
                     filters = 'ensembl_gene_id',
                     values = ensembl_genes,
                     mart = ensembl)

# Changing genes names
library(dplyr)
colnames(gene_symbols)[2] = 'gene_id'
exp_set = left_join(gene_symbols, exp_set, by='gene_id')

# Removing ensembl_ids
exp_set$gene_id = NULL

# Removing empty rows
exp_set <- exp_set[!exp_set$hgnc_symbol == "",]

# Removing .counts from the column names
names(exp_set) = gsub("\\..*","", names(exp_set))


###################### Collapsing genes by  sum ###################### 
# Checking for duplicated genes
n_occur <- data.frame(table(exp_set$hgnc_symbol))
n_occur[n_occur$Freq > 1,]

# Sum duplicated genes to normalize the expression set
exp_set = aggregate(. ~ hgnc_symbol, data=exp_set, FUN=sum)

# define genes as row names and delete the gene_symbol column
row.names(exp_set) = exp_set$hgnc_symbol; exp_set$hgnc_symbol = NULL

# Clean environment
rm(ensembl_genes, ensembl, gene_symbols, n_occur)

write.table(exp_set,
            file = 'exp_set.tsv',
            sep = '\t',
            row.names = T)

write.table(metatable,
            file = 'samplesinfo.tsv',
            sep = '\t',
            row.names = T)


# Normalization with DESeq
library(DESeq2)

metatable <- metatable[metatable$characteristics_ch1.7 %in% c("disease state: match control",
                                                              "disease state: psychosis_bipolar"),]

meta <- metatable
meta$`age:ch1` <- as.numeric(meta$`age:ch1`)
meta$age <- (meta$`age:ch1` - mean(meta$`age:ch1`) ) / sd(meta$`age:ch1`)
meta$`rin:ch1` <- as.numeric(meta$`rin:ch1`)
meta$rin <- (meta$`rin:ch1` - mean(meta$`rin:ch1`) ) / sd(meta$`rin:ch1`)
meta$`pmi:ch1` <- as.numeric(meta$`pmi:ch1`)
meta$pmi <- (meta$`pmi:ch1` - mean(meta$`pmi:ch1`) ) / sd(meta$`pmi:ch1`)
meta$`ph:ch1` <- as.numeric(meta$`ph:ch1`)
meta $ph <- (meta$`ph:ch1` - mean(meta$`ph:ch1`) ) / sd(meta$`ph:ch1`)
meta$manner <- as.factor(meta$`manner of death:ch1`)
meta$tissue <- as.factor(meta$`tissue:ch1`)

exp_set <- exp_set[, meta$sample]

# Subset count matrix to keep samples of the chosen class and controls
cnt <- as.matrix(exp_set[,colnames(exp_set) %in% meta$sample])

cl <- "disease state: psychosis_bipolar"

# Define the control group label
control_label <- "disease state: match control"

meta$class <- factor(meta$characteristics_ch1.7, levels = c(control_label, cl))


dds <- DESeqDataSetFromMatrix(countData = cnt,
                              colData = meta,
                              design = ~ ph + pmi + rin + age + tissue + manner + class,
                              tidy = F)

# Plot PCA based on diagnostic
vsd <- vst(dds, blind=FALSE)

norm <- as.data.frame(assay(vsd))


###################### Using MDP function on data normalized with DESeq ####
library(mdp)

metatable <- metatable %>% rownames_to_column(var= "Sample")
rownames(metatable) <- metatable$Sample

samplesinfo <- metatable[, c("Sample","characteristics_ch1.7")]
colnames(samplesinfo) <- c("Sample", 
                           "Class")

norm <- norm[, rownames(samplesinfo)]

mdp_deseq <- mdp(data = norm,
                 pdata = samplesinfo, 
                 control_lab = "disease state: match control")

###################### Plot boxplot for all the samples ####
library(ggplot2)

boxplot_deseq <- mdp_deseq$sample_scores$allgenes
boxplot_deseq <- boxplot_deseq[boxplot_deseq$Class == "disease state: match control",]
boxplot_deseq %>% 
  ggplot(aes(x = Class, y = zscore_class, fill = Class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "black") +
  labs(title = "All genes - GSE202537") + 
  theme_bw() +
  theme(legend.position = "none")
boxplot_deseq <- mdp_deseq$sample_scores$perturbedgenes
boxplot_deseq <- boxplot_deseq[boxplot_deseq$Class == "disease state: match control",]

boxplot_deseq %>% 
  ggplot(aes(x = Class, y = zscore_class, fill = Class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "black") +
  labs(title = "Perturbed genes - GSE202537") +
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


#####################################################
cl <- "disease state: psychosis_bipolar"

# Define the control group label
control_label <- "disease state: match control"

meta <- metatable %>%
  filter(characteristics_ch1.7 %in% c(cl, control_label))# FIRST CONTROL, SECOND TREATED.
meta$class <- factor(meta$characteristics_ch1.7, levels = c(control_label, cl))

# Subset count matrix to keep samples of the chosen class and controls
cnt <- as.matrix(exp_set[,colnames(exp_set) %in% meta$sample])

# create dds object
meta$`age:ch1` <- as.numeric(meta$`age:ch1`)
meta$age <- (meta$`age:ch1` - mean(meta$`age:ch1`) ) / sd(meta$`age:ch1`)
meta$`rin:ch1` <- as.numeric(meta$`rin:ch1`)
meta$rin <- (meta$`rin:ch1` - mean(meta$`rin:ch1`) ) / sd(meta$`rin:ch1`)
meta$`pmi:ch1` <- as.numeric(meta$`pmi:ch1`)
meta$pmi <- (meta$`pmi:ch1` - mean(meta$`pmi:ch1`) ) / sd(meta$`pmi:ch1`)
meta$`ph:ch1` <- as.numeric(meta$`ph:ch1`)
meta $ph <- (meta$`ph:ch1` - mean(meta$`ph:ch1`) ) / sd(meta$`ph:ch1`)
meta$manner <- as.factor(meta$`manner of death:ch1`)
meta$tissue <- as.factor(meta$`tissue:ch1`)
meta$gender <- as.factor(meta$`gender:ch1`)
dds <- DESeqDataSetFromMatrix(countData = cnt,
                              colData = meta,
                              design = ~ age + gender + rin + pmi + ph + manner + tissue + class,
                              tidy = F)

# Remove genes that have zero counts in all samples (i.e. not expressed genes)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

# Run DE analysis with default DESeq2 settings
dds <- DESeq(dds)
ddsClean <- dds[which(mcols(dds)$betaConv),]

# extract results defining the order of comparison
# must be: constrast=c("CLASS COLUMN IN META","DISEASE/TREATMENT GROUP","CONTROL)
# This will give log2FoldChange values that correspond to disease vs control
res <- as.data.frame(results(ddsClean, contrast = c("class", cl, control_label)))
res <- res[res$padj < 0.05,]
res <- res[res$log2FoldChange > 1 | res$log2FoldChange < -1,]
res <- res[!is.na(res$padj),]

write.table(res,
            "DEGs_BDvsCtrl_filtered.tsv",
            row.names = T)
