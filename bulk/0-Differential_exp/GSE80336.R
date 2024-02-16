setwd('~/Documents/IC/Davis_2016_MolPsychiatry')

###################### Loading ###################### 
library(data.table)
exp_set = data.table::fread('~/Documents/IC/Davis_2016_MolPsychiatry/data/GSE80336/GSE80336_Counts.txt')
exp_set <- exp_set[,-c(1,3,4)]
###################### Collapsing genes by  sum ###################### 
# Checking for duplicated genes
n_occur <- data.frame(table(exp_set$GeneSymbol))
n_occur[n_occur$Freq > 1,]

# Sum duplicated genes to normalize the expression set
exp_set = aggregate(. ~ GeneSymbol, data=exp_set, FUN=sum)

# define genes as row names and delete the gene_symbol column
row.names(exp_set) = exp_set$GeneSymbol; exp_set$GeneSymbol = NULL

# Clean environment
rm(ensembl_genes, ensembl, gene_symbols, n_occur)

###################### Setup Meta and Expression ###################### 
samplesinfo <- read.csv('~/Documents/IC/Davis_2016_MolPsychiatry/samplesinfo.csv')
# Identify the sample names in metadata
head(colnames(exp_set))
head(table(samplesinfo$title))

# Filter samples info
names(samplesinfo)

# Subset samplesinfo object to keep only samples that are 
# present in the expression matrix
samplesinfo$Sample = colnames(exp_set)
print(paste("Total number of samples:",
            nrow(samplesinfo[samplesinfo$Sample %in% colnames(exp_set),])))
samplesinfo = samplesinfo[samplesinfo$Sample %in% colnames(exp_set),]

# Check if all are in the same order
all(samplesinfo$Sample %in% colnames(exp_set))

# Set the classes you need to compare
samplesinfo$Class = samplesinfo$source_name_ch1


write.table(exp_set,
            file = 'intermediate/exp_set.tsv',
            sep = '\t',
            row.names = T)

write.table(samplesinfo,
            file = 'intermediate/samplesinfo.tsv',
            sep = '\t',
            row.names = T)



############################################################################### 

setwd('~/Documents/IC/Davis_2016_MolPsychiatry')
library(dplyr)
library(tibble)

exp_set = data.table::fread('intermediate/exp_set.tsv')
exp_set <- exp_set %>% remove_rownames %>% column_to_rownames(var="V1")

samplesinfo = data.table::fread('intermediate/samplesinfo.tsv')
samplesinfo <- samplesinfo %>% remove_rownames %>% column_to_rownames(var="V1")
rownames(samplesinfo) <-samplesinfo$Sample


# Define sample classes in a metatable
metatable <- data.frame(sample_name=samplesinfo$Sample, # Name of the samples
                        class=samplesinfo$Status,
                        sex=samplesinfo$Sex,
                        interval=samplesinfo$PMI,
                        age=samplesinfo$Age,
                        rin=samplesinfo$RIN)
metatable <- metatable %>% remove_rownames %>% column_to_rownames(var="sample_name")
metatable <- data.frame(metatable,
                        sample_name=samplesinfo$Sample)


# Normalization with DESeq
library(DESeq2)

metatable$interval <- as.numeric(metatable$interval)
metatable$rin <- as.numeric(metatable$rin)

metatable$age <- (metatable$age - mean(metatable$age) ) / sd(metatable$age)
metatable$interval <- (metatable$interval - mean(metatable$interval) ) / sd(metatable$interval)
metatable$rin <- (metatable$rin - mean(metatable$rin) ) / sd(metatable$rin)
dds <- DESeqDataSetFromMatrix(countData = exp_set,
                              colData = metatable, 
                              design = ~ age + sex + interval + rin + class)


vsd <- vst(dds, blind=FALSE)
norm <- as.data.frame(assay(vsd))




###################### Using MDP function on data normalized with DESeq ####
library(mdp)
library(tidyverse)

samplesinfo <- samplesinfo[,colnames(samplesinfo) %in% c("Sample", "Status")]

colnames(samplesinfo) <- c("Sample", 
                           "Class")

mdp_deseq <- mdp(data = norm,
                 pdata = samplesinfo, 
                 control_lab = "control")

###################### Plot boxplot for all the samples ####
library(ggplot2)

boxplot_deseq <- mdp_deseq$sample_scores$allgenes
boxplot_deseq <- boxplot_deseq[boxplot_deseq$Class == "control",]
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

outliers_deseq <- c()

for(i in which(mdp_deseq$sample_scores$perturbedgenes$outlier == 1)){
  outliers_deseq <- c(outliers_deseq, mdp_deseq$sample_scores$perturbedgenes$Sample[i])
}
print(paste0("Number of outlier samples considering perturbed genes + DESeq normalization: ", length(outliers_deseq)))

exp_set <- exp_set[, colnames(exp_set) != "C_28"]
norm <- norm[, colnames(norm) != "C_28"]
metatable <- metatable[colnames(exp_set),]

write.table(norm,
            file = 'deconvolution/norm.tsv',
            sep = '\t',
            row.names = TRUE)

###################### DEGs with DESeq2
###################### BP vs Ctrl ####
samplesinfo = data.table::fread('samplesinfo.csv')

# Define sample classes in a metatable
# Define sample classes in a metatable
metatable <- data.frame(sample_name=samplesinfo$Sample, # Name of the samples
                        class=samplesinfo$Status,
                        sex=samplesinfo$Sex,
                        interval=samplesinfo$PMI,
                        age=samplesinfo$Age,
                        rin=samplesinfo$RIN)
metatable <- metatable %>% remove_rownames %>% column_to_rownames(var="sample_name")
metatable <- data.frame(metatable,
                        sample_name=samplesinfo$Sample)
metatable <- metatable[colnames(exp_set),]

############### DEGs BD vs Ctrl #################
classes <- "bipolar"
control_label <- "control"
  
meta <- metatable

# Subset count matrix to keep samples of the chosen class and controls
cnt <- as.matrix(exp_set[,colnames(exp_set) %in% meta$sample_name])
meta$age <- (meta$age - mean(meta$age) ) / sd(meta$age)
meta$interval <- (meta$interval - mean(meta$interval) ) / sd(meta$interval)
meta$rin <- (meta$rin - mean(meta$rin) ) / sd(meta$rin)
meta$class <- as.factor(meta$class)
# create dds object. There is no BD patient from cohort B
dds_trat <- DESeqDataSetFromMatrix(countData = cnt,
                                   colData = meta,
                                   design = ~ rin + interval + class,
                                   tidy = F)

# Remove genes that have zero counts in all samples (i.e. not expressed genes)
keep <- rowSums(counts(dds_trat)) > 1
dds_trat <- dds_trat[keep,]

# Run DE analysis with default DESeq2 settings
dds_trat <- DESeq(dds_trat)

# extract results defining the order of comparison
# must be: constrast=c("CLASS COLUMN IN META","DISEASE/TREATMENT GROUP","CONTROL)
# This will give log2FoldChange values that correspond to disease vs control
res <- as.data.frame(results(dds_trat, contrast = c("class",
                                                    classes,
                                                    control_label),cooksCutoff = F))
res <- res[res$padj < 0.05 & (res$log2FoldChange < -1 | res$log2FoldChange > 1),]
res <- res[!is.na(res$padj),]

write.table(x = res,
            file = "DEGs_BDvsCtrl.tsv",
            sep="\t")
