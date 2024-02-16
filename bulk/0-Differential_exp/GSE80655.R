setwd('~/Documents/IC/Ramaker_2017_GenomeMed')

###################### Loading ###################### 
library(data.table)
exp_set = data.table::fread('~/Documents/IC/Ramaker_2017_GenomeMed/GSE80655/GSE80655_GeneExpressionData_Updated_3-26-2018.txt')

###################### Get BioMart ###################### 
library(biomaRt)

# Load Database
ensembl = useMart(biomart="ensembl",
                  dataset="hsapiens_gene_ensembl")

# Ensembl gene IDs
ensembl_genes = exp_set$gene_id
ensembl_genes = gsub("\\..*","", ensembl_genes)
exp_set$gene_id = ensembl_genes

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

###################### Setup Meta and Expression ###################### 
samplesinfo <- read.table('~/Documents/IC/Ramaker_2017_GenomeMed/data/GSE80655/table_phenodata.tsv', header = T)
# Identify the sample names in metadata
head(colnames(exp_set))
head(table(samplesinfo$title))

# Filter samples info
names(samplesinfo)
samplesinfo <- samplesinfo[,
                           c("description",
                             "source_name_ch1",
                             "brain.region.ch1",
                             "clinical.diagnosis.ch1",
                             "age.at.death.ch1",
                             "ethnicity.ch1",
                             "post.mortem.interval.ch1",
                             "gender.ch1",
                             "brain.ph.ch1"
                           )]

colnames(samplesinfo)[which(names(samplesinfo) == 'description')] <- 'Sample'

# Subset samplesinfo object to keep only samples that are 
# present in the expression matrix
print(paste("Total number of samples:",
            nrow(samplesinfo[samplesinfo$Sample %in% colnames(exp_set),])))
samplesinfo = samplesinfo[samplesinfo$Sample %in% colnames(exp_set),]

# Check if all are in the same order
all(samplesinfo$Sample %in% colnames(exp_set))

# Set the classes you need to compare
samplesinfo$Class = samplesinfo$source_name_ch1

# Exclude Schizophrenia and Major Depression subjects
samplesinfo <- samplesinfo[!samplesinfo$clinical.diagnosis.ch1 %in% c("Major Depression", "Schizophrenia"),]
exp_set <- exp_set[, samplesinfo$Sample]
all(samplesinfo$Sample %in% colnames(exp_set))

write.table(exp_set,
            file = 'intermediate/exp_set.tsv',
            sep = '\t',
            row.names = T)

write.table(samplesinfo,
            file = 'intermediate/samplesinfo.tsv',
            sep = '\t',
            row.names = T)

###############################################################################

setwd('~/Documents/IC/Ramaker_2017_GenomeMed')
library(dplyr)
library(tibble)

exp_set = data.table::fread('intermediate/exp_set.tsv')
exp_set <- exp_set %>% remove_rownames %>% column_to_rownames(var="V1")

samplesinfo = data.table::fread('intermediate/samplesinfo.tsv')
samplesinfo <- samplesinfo %>% remove_rownames %>% column_to_rownames(var="V1")

# Define sample classes in a metatable
metatable <- data.frame(sample_name=samplesinfo$Sample, # Name of the samples
                        class=samplesinfo$Class,
                        sex=samplesinfo$gender.ch1,
                        region=samplesinfo$brain.region.ch1,
                        diagnosis=samplesinfo$clinical.diagnosis.ch1,
                        ethnicity=samplesinfo$ethnicity.ch1,
                        interval=samplesinfo$post.mortem.interval.ch1,
                        age=samplesinfo$age.at.death.ch1,
                        ph=samplesinfo$brain.ph.ch1)
metatable <- metatable %>% remove_rownames %>% column_to_rownames(var="sample_name")
metatable <- data.frame(metatable,
                        sample_name=samplesinfo$Sample)


# Normalization with DESeq
library(DESeq2)

# Very important step! Get standard deviation for each continuous variable
meta <- metatable
meta$age <- meta$age / sd(meta$age)
meta$interval <- meta$interval / sd(meta$interval)
meta$ph <- meta$ph / sd(meta$ph)
meta$diagnosis <- factor(meta$diagnosis)
meta$ethnicity <- factor(meta$ethnicity)
dds <- DESeqDataSetFromMatrix(countData = exp_set,
                              colData = meta, 
                              design = ~ age + sex + region + ethnicity + interval + diagnosis)


vsd <- vst(dds, blind=FALSE)
norm <- assay(vsd)


###################### Using MDP function on data normalized with DESeq ####
library(mdp)

metatable <- metatable %>% rownames_to_column(var= "Sample")
rownames(metatable) <- metatable$Sample

samplesinfo <- metatable[, c("Sample","diagnosis")]
colnames(samplesinfo) <- c("Sample", 
                           "Class")

norm <- norm[, rownames(samplesinfo)]
norm <- as.data.frame(norm)

mdp_deseq <- mdp(data = norm,
                 pdata = samplesinfo, 
                 control_lab = "Control")

###################### Plot boxplot for all the samples ####
library(ggplot2)

boxplot_deseq <- mdp_deseq$sample_scores$allgenes
boxplot_deseq <- boxplot_deseq[boxplot_deseq$Class == "Control",]

boxplot_deseq %>% 
  ggplot(aes(x = Class, y = zscore_class,fill=Class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "black") +
  labs(title = "All genes - GSE80655") +
  theme_bw() +
  theme(legend.position = "none")
boxplot_deseq <- mdp_deseq$sample_scores$perturbedgenes
boxplot_deseq <- boxplot_deseq[boxplot_deseq$Class == "Control",]

boxplot_deseq %>% 
  ggplot(aes(x = Class, y = zscore_class,fill=Class)) + 
  geom_boxplot(width = .2, outlier.colour = NA, coef = 1000) + 
  geom_jitter(width = 0.05, alpha = 0.4, color = "black") +
  labs(title = "Perturbed genes - GSE80655")+
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

write.table(exp_set,
            file = 'exp_set_after_mdp.tsv',
            sep = '\t',
            row.names = TRUE)


metatable <- metatable[metatable$region != "nAcc",]
exp_set <- exp_set[,colnames(exp_set) %in% rownames(metatable)]


meta <- metatable
meta$age <- meta$age / sd(meta$age)
meta$interval <- meta$interval / sd(meta$interval)
meta$ph <- meta$ph / sd(meta$ph)
meta$diagnosis <- factor(meta$diagnosis)
meta$ethnicity <- factor(meta$ethnicity)
dds <- DESeqDataSetFromMatrix(countData = exp_set,
                              colData = meta, 
                              design = ~ age + sex + region + ethnicity + interval + diagnosis)

vsd <- vst(dds, blind=FALSE)
norm <- assay(vsd)
cpm <- fpm(dds, robust = FALSE)

write.table(norm,
            file = 'deconvolution/DLPFC/cortex_norm.tsv',
            sep = '\t',
            row.names = TRUE)

write.table(cpm,
            file = 'deconvolution/DLPFC/cpm_cortex.tsv',
            sep = '\t',
            row.names = TRUE)


####################################

classes <- "Bipolar Disorder"

control_label <- "Control"

exp_set <- as.data.frame(exp_set)

# Subset metadata to include only samples of the chosen class and controls
meta <- metatable %>%
  filter(diagnosis == control_label | diagnosis == classes)

# Transform class column into a factor with the levels in the correct order:
# FIRST CONTROL, SECOND TREATED.
meta$class <- factor(meta$class, levels = c(control_label, classes))

# Subset count matrix to keep samples of the chosen class and controls
cnt <- as.matrix(exp_set[,colnames(exp_set) %in% meta$sample_name])

# create dds object. 
library(DESeq2)
meta$age <- as.numeric(meta$age)
meta$age <- (meta$age - mean(meta$age) ) / sd(meta$age)
meta$interval <- as.numeric(meta$interval)
meta$interval <- (meta$interval - mean(meta$interval) ) / sd(meta$interval)
meta$ph <- as.numeric(meta$ph)
meta$ph <- (meta$ph - mean(meta$ph) ) / sd(meta$ph)
meta$region <- as.factor(meta$region)
meta$ethnicity <- as.factor(meta$ethnicity)
meta$sex <- as.factor(meta$sex)
dds <- DESeqDataSetFromMatrix(countData = cnt,
                              colData = meta,
                              design = ~ age + interval + sex + region + diagnosis,
                              tidy = F)

# Remove genes that have zero counts in all samples (i.e. not expressed genes)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]


# Run DE analysis with default DESeq2 settings
dds_trat <- DESeq(dds)
ddsClean <- dds_trat[which(mcols(dds_trat)$betaConv),]

# extract results defining the order of comparison
# must be: constrast=c("CLASS COLUMN IN META","DISEASE/TREATMENT GROUP","CONTROL)
# This will give log2FoldChange values that correspond to disease vs control
res <- as.data.frame(results(ddsClean, contrast = c("diagnosis",
                                                    classes,
                                                    control_label)))

res <- res[res$padj < 0.05,]
res <- res[res$log2FoldChange < -1 | res$log2FoldChange > 1,]
res <- res[!is.na(res$padj),]

write.table(res,
            file = 'intermediate/DEGs_BD_vs_Ctrl.tsv',
            sep = '\t',
            row.names = T,
            quote = F)
