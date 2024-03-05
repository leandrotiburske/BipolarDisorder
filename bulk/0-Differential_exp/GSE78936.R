setwd("~/Documents/IC/GSE78936_FC/")

library(tidyverse)

metatable <- read.csv("diagnosis.csv")
exp_set <- read.delim("GSE78936_countMatrix.txt.gz")

###################### Get BioMart ###################### 
library(biomaRt)

# Load Database
ensembl = useMart(biomart="ensembl",
                  dataset="hsapiens_gene_ensembl")

# Ensembl gene IDs
ensembl_genes = rownames(exp_set)
#ensembl_genes = gsub("\\..*","", ensembl_genes)
exp_set$gene_id = ensembl_genes

# Dataframe of attributes if you want to check
# at = ensembl@attributes

# Get gene_symbols from Ensembl(BioMart)
gene_symbols = getBM(attributes = c('hgnc_symbol',
                                    'entrezgene_id'),
                     filters = 'entrezgene_id',
                     values = ensembl_genes,
                     mart = ensembl)

# Changing genes names
library(dplyr)
exp_set$gene_id <- rownames(exp_set)
colnames(gene_symbols)[2] = 'gene_id'
class(gene_symbols$gene_id) <- "character"
exp_set <- exp_set[exp_set$gene_id %in% gene_symbols$gene_id,]
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

# Normalization with DESeq
library(DESeq2)

metatable <- metatable[metatable$Diagnosis != "SZ",]
exp_set <- exp_set[, colnames(exp_set) %in% metatable$Sample]

dds <- DESeqDataSetFromMatrix(countData = exp_set,
                              colData = metatable, 
                              design = ~ Region + Diagnosis)


vsd <- vst(dds, blind=FALSE)
norm <- as.data.frame(assay(vsd))




###################### Using MDP function on data normalized with DESeq ####
library(mdp)
library(tidyverse)

metadata <- metatable[,colnames(metatable) %in% c("Sample", "Diagnosis")]

colnames(metadata) <- c("Sample", 
                           "Class")


mdp_deseq <- mdp(data = norm,
                 pdata = metadata, 
                 control_lab = "Ctrl")

###################### Plot boxplot for all the samples ####
library(ggplot2)

boxplot_deseq <- mdp_deseq$sample_scores$allgenes
boxplot_deseq <- boxplot_deseq[boxplot_deseq$Class == "dorsal striatum, control",]
boxplot_deseq %>% 
  ggplot(aes(x = Class, y = zscore_class,fill=Class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "black") +
  labs(title = "All genes - GSE80336") +
  theme_bw() +
  theme(legend.position = "none")
boxplot_deseq <- mdp_deseq$sample_scores$perturbedgenes
boxplot_deseq <- boxplot_deseq[boxplot_deseq$Class == "dorsal striatum, control",]
boxplot_deseq %>% 
  ggplot(aes(x = Class, y = zscore_class,fill=Class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "black") +
  labs(title = "Perturbed genes - GSE80336")+
  theme_bw() +
  theme(legend.position = "none")

## Get outlier samples

outliers_deseq <- c()

for(i in which(mdp_deseq$sample_scores$allgenes$outlier == 1)){
  outliers_deseq <- c(outliers_deseq, mdp_deseq$sample_scores$allgenes$Sample[i])
}
print(paste0("Number of outlier samples considering all genes + DESeq normalization: ", length(outliers_deseq)))

exp_set <- exp_set[, !colnames(exp_set) %in% outliers_deseq]
norm <- norm[, !colnames(norm) %in% outliers_deseq]
metatable <- metatable[metatable$Sample %in% colnames(exp_set),]

outliers_deseq <- c()

for(i in which(mdp_deseq$sample_scores$perturbedgenes$outlier == 1)){
  outliers_deseq <- c(outliers_deseq, mdp_deseq$sample_scores$perturbedgenes$Sample[i])
}
print(paste0("Number of outlier samples considering perturbed genes + DESeq normalization: ", length(outliers_deseq)))

exp_set <- exp_set[, !colnames(exp_set) %in% outliers_deseq]
norm <- norm[, !colnames(norm) %in% outliers_deseq]
metatable <- metatable[metatable$Sample %in% colnames(exp_set),]


dds_trat <- DESeqDataSetFromMatrix(countData = exp_set,
                                   colData = metatable,
                                   design = ~ Region + Diagnosis,
                                   tidy = F)

cpm <- fpm(dds_trat, robust = F)
vsd <- vst(dds, blind=FALSE)
norm <- as.data.frame(assay(vsd))

write.table(norm,
            file = 'deconvolution/norm.tsv',
            sep = '\t',
            row.names = TRUE)

write.table(cpm,
            file = 'deconvolution/cpm.tsv',
            sep = '\t',
            row.names = TRUE)

# Remove genes that have zero counts in all samples (i.e. not expressed genes)
keep <- rowSums(counts(dds_trat)) > 1
dds_trat <- dds_trat[keep,]

# Run DE analysis with default DESeq2 settings
dds_trat <- DESeq(dds_trat)

# extract results defining the order of comparison
# must be: constrast=c("CLASS COLUMN IN META","DISEASE/TREATMENT GROUP","CONTROL)
# This will give log2FoldChange values that correspond to disease vs control
classes <- "BD"
control_label <- "Ctrl"

res <- as.data.frame(results(dds_trat, contrast = c("Diagnosis",
                                                    classes,
                                                    control_label)))
res <- res[res$padj < 0.05 & (res$log2FoldChange < -1 | res$log2FoldChange > 1),]
res <- res[!is.na(res$padj),]

write.table(x = res,
            file = "DEGs_BDvsCtrl.tsv",
            sep="\t")
