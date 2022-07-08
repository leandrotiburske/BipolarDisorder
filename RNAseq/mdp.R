###################### Changing directory ####

setwd('~/Documents/IC/bpd_2020_krebs_phycomed/bpd_2020_krebs_phycomed')

###################### Loading ###################### 
library(data.table)
exp_set = data.table::fread('data/GSE124326/GSE124326_count_matrix.txt')

###################### Get BioMart ###################### 
library(biomaRt)

# Load Database
ensembl = useMart(biomart="ensembl",
                  dataset="hsapiens_gene_ensembl")

# Ensembl gene IDs
ensembl_genes = exp_set$gene
ensembl_genes = gsub("\\..*","", ensembl_genes)
exp_set$gene = ensembl_genes

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
colnames(gene_symbols)[2] = 'gene'
exp_set = left_join(gene_symbols, exp_set, by='gene')

# Removing ensembl_ids
exp_set$gene = NULL

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

###################### Setup Meta and Expression ###################### 
samplesinfo <- read.table('data/GSE124326/table_phenodata.tsv', header = T)

# Identify the sample names in metadata
head(colnames(exp_set))
head(table(samplesinfo$title))

# Filter samples info
names(samplesinfo)
samplesinfo <- samplesinfo[,
                           c("title",
                             "geo_accession",
                             "source_name_ch1",
                             "age.ch1",
                             "b.cells.memory.ch1",
                             "b.cells.naive.ch1",
                             "bipolar.disorder.diagnosis.ch1",
                             "dendritic.cells.activated.ch1",
                             "dendritic.cells.resting.ch1",
                             "eosinophils.ch1",
                             "included.in.final.analysis.ch1",
                             "lithium.use..non.user.0..user...1..ch1",
                             "macrophages.m0.ch1",
                             "macrophages.m1.ch1",
                             "macrophages.m2.ch1",
                             "mast.cells.activated.ch1",
                             "mast.cells.resting.ch1",
                             "monocytes.ch1",
                             "neutrophils.ch1",
                             "nk.cells.activated.ch1",
                             "nk.cells.resting.ch1",
                             "plasma.cells.ch1",
                             "rin.ch1",
                             "Sex.ch1",
                             "t.cells.cd4.memory.activated.ch1",
                             "t.cells.cd4.memory.resting.ch1",
                             "t.cells.cd4.naive.ch1",
                             "t.cells.cd8.ch1",
                             "t.cells.follicular.helper.ch1",
                             "t.cells.gamma.delta.ch1",
                             "t.cells.regulatory.ch1",
                             "tissue.cell.type.ch1",
                             "tobacco.use.ch1"
                           )]

colnames(samplesinfo)[which(names(samplesinfo) == 'geo_accession')] <- 'Sample'
colnames(exp_set) <- samplesinfo$Sample

# Subset samplesinfo object to keep only samples that are 
# present in the expression matrix
print(paste("Total number of samples:",
            nrow(samplesinfo[samplesinfo$Sample %in% colnames(exp_set),])))
samplesinfo = samplesinfo[samplesinfo$Sample %in% colnames(exp_set),]

# reorder expression matrix according to the order of the samplesinfo object
# VERY IMPORTANT!! ALWAYS DO THIS!!
exp_set <- exp_set[, samplesinfo$Sample]

# Check if all are in the same order
all(samplesinfo$Sample %in% colnames(exp_set))

# Set the classes you need to compare
samplesinfo$Class = samplesinfo$bipolar.disorder.diagnosis.ch1

# Change characters in sex column to numbers so that removeBatchEffect() can analyse it
samplesinfo$Sex.ch1[samplesinfo$Sex.ch1 == "F"] <- 0
samplesinfo$Sex.ch1[samplesinfo$Sex.ch1 == "M"] <- 1

library(tidyverse)

drop_indexes <- which(is.na(samplesinfo$Sex.ch1))
samplesinfo <- samplesinfo %>% drop_na(Sex.ch1)
exp_set <- subset(exp_set, select = -c(drop_indexes))

metatable <- data.frame(sample_name=samplesinfo$Sample, # Name of the samples
                        class=samplesinfo$Class,
                        tobacco_use=samplesinfo$tobacco.use.ch1,
                        lithium_use=samplesinfo$lithium.use..non.user.0..user...1..ch1,
                        sex=samplesinfo$Sex.ch1,
                        rin=samplesinfo$rin.ch1) # Classes
metatable <- metatable %>% remove_rownames %>% column_to_rownames(var="sample_name")
metatable <- data.frame(metatable,
                        sample_name=samplesinfo$Sample)


###################### Normalization with DESeq ####
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = exp_set, colData = metatable, design = ~ class)
dds <- estimateSizeFactors(dds)
normalized_deseq <- counts(dds, normalized=TRUE)
vsd <- vst(dds, blind=FALSE)

library(limma)

assay(vsd) <- limma::removeBatchEffect(assay(vsd),
                                       batch = vsd$rin,
                                       batch2 = vsd$sex,
                                       batch3 = vsd$tobacco_use,
                                       batch4 = vsd$lithium_use)
vsd <- assay(vsd)
normalized_deseq <- as.data.frame(vsd)

###################### Normalization with edgeR ####
library(edgeR)

dds <- DESeqDataSetFromMatrix(countData = exp_set, colData = samplesinfo, design = ~ Class)
normalized_tmm <- calcNormFactors(dds, method = "TMM")
tmm <- cpm(normalized_tmm)
tmm <- as.data.frame(tmm)

###################### Using MDP function on data normalized with DESeq ####
library(mdp)

mdp_deseq <- mdp(normalized_deseq, samplesinfo, control_lab = "Control")

###################### Plot boxplot for all the samples ####
library(ggplot2)

boxplot_deseq <- mdp_deseq$sample_scores$allgenes
boxplot_deseq %>% 
  ggplot(aes(x = Class, y = zscore_class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "orange") +
  labs(title = "Considering all genes and DESeq normalization")
boxplot_deseq <- mdp_deseq$sample_scores$perturbedgenes
boxplot_deseq %>% 
  ggplot(aes(x = Class, y = zscore_class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "orange") +
  labs(title = "Considering perturbed genes and DESeq normalization")

## Get outlier samples

outliers_deseq <- c()

for(i in which(mdp_deseq$sample_scores$allgenes$outlier == 1)){
  outliers_deseq <- c(outliers_deseq, mdp_deseq$sample_scores$allgenes$Sample[i])
}
print(paste0("Number of outlier samples considering all genes + DESeq normalization: ", length(outliers_deseq)))

outliers_deseq <- c()

for(i in which(mdp_deseq$sample_scores$perturbedgenes$outlier == 1)){
  outliers_deseq <- c(outliers_deseq, mdp_deseq$sample_scores$perturbedgenes$Sample[i])
}
print(paste0("Number of outlier samples considering perturbed genes + DESeq normalization: ", length(outliers_deseq)))

###################### Use MDP function on data normalized with edgeR ####

mdp_tmm <- mdp(tmm, samplesinfo, control_lab = "Control")

###################### Plot boxplot for all samples ####

boxplot_tmm <- mdp_tmm$sample_scores$allgenes
boxplot_tmm %>% 
  ggplot(aes(x = Class, y = zscore_class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "orange") +
  labs(title = "Considering all genes and TMM normalization")
boxplot_tmm <- mdp_tmm$sample_scores$perturbedgenes
boxplot_tmm %>% 
  ggplot(aes(x = Class, y = zscore_class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "orange") +
  labs(title = "Considering perturbed genes and TMM normalization")

###################### Get outlier samples ####

outliers_tmm <- c()

for(i in which(mdp_tmm$sample_scores$perturbedgenes$outlier == 1)){
  outliers_tmm <- c(outliers_tmm, mdp_tmm$sample_scores$perturbedgenes$Sample[i])
}
print(paste0("Number of outlier samples considering perturbed genes + TMM normalization: ", length(outliers_tmm)))

outliers_tmm <- c()

for(i in which(mdp_tmm$sample_scores$allgenes$outlier == 1)){
  outliers_tmm <- c(outliers_tmm, mdp_tmm$sample_scores$allgenes$Sample[i])
}

print(paste0("Number of outlier samples considering all genes + TMM normalization: ", length(outliers_tmm)))