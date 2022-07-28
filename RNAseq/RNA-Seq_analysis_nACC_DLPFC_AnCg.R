#### Transcriptome Study of Bipolar Disorder Sexual Dimorphism ----

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

# Change characters in sex column to numbers so that removeBatchEffect() can analyse it
samplesinfo$gender.ch1[samplesinfo$Sex.ch1 == "F"] <- 0
samplesinfo$gender.ch1[samplesinfo$Sex.ch1 == "M"] <- 1


write.table(exp_set,
            file = 'intermediate/exp_set.tsv',
            sep = '\t',
            row.names = T)

write.table(samplesinfo,
            file = 'intermediate/samplesinfo.tsv',
            sep = '\t',
            row.names = T)






###################### Plot PCAs ##################################### 

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

metatable$age <- (metatable$age - mean(metatable$age) ) / sd(metatable$age)
metatable$interval <- (metatable$interval - mean(metatable$interval) ) / sd(metatable$interval)
metatable$ph <- (metatable$ph - mean(metatable$ph) ) / sd(metatable$ph)
metatable$class <- factor(metatable$class)
metatable$ethnicity <- factor(metatable$ethnicity)
dds <- DESeqDataSetFromMatrix(countData = exp_set,
                              colData = metatable, 
                              design = ~ age + ethnicity + interval + ph + class)

library(FactoMineR)
library(limma)

# Plot PCA based on diagnostic
vsd <- vst(dds, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd),
                                       batch = vsd$age,
                                       batch1 = vsd$interval,
                                       batch2 = vsd$sex,
                                       batch3 = vsd$rin)

PCA <- plotPCA(vsd, intgroup = 'class')
print(PCA)

# Plot PCA using fast.prcomp() for comparison
library(ggplot2)
library(gmodels)

bpctrl <- assay(vsd)
bpctrl <- as.data.frame(bpctrl)
PCA <- fast.prcomp(t(bpctrl))
pca_plot <- data.frame(Sample = row.names(PCA$x),
                       PC1 = PCA$x[, 1], # selecting PC1
                       PC2 = PCA$x[, 2], # selecting PC2
                       metatable)

pca_var <- PCA$sdev ^ 2
pca_var_per <- round(pca_var / sum(pca_var) * 100, 1) # calculando a porcentagem de variação de cada PC
percentage <- paste0(colnames(PCA$x), " ", as.character(pca_var_per), "%") # juntando os termos PC e suas respectivas porcentagens

## Based on diagnostic
ggplot(pca_plot,
       aes(x = PC1, 
           y = PC2, 
           color = class)) +
  geom_point() +
  xlab(percentage[1]) +
  ylab(percentage[2])

rm(bpctrl, PCA, pca_plot, vsd, pca_var, pca_var_per, percentage)


###################### DEGs with DESeq2
###################### BP vs Ctrl ####

# Define the control group label
control_label <- "dorsal striatum, control"

# Create sexes variable to iterate through
sexes <- unique(metatable$sex)

# Define the targets/conditions (all classes but control/healthy)
classes <- "dorsal striatum, bipolar"

for(i in 1:length(sexes)){
  # Define class to be compared in DE analysis
  sx <- sexes[i]
  
  # Subset metadata to include only samples of the chosen class and controls
  meta <- metatable %>%
    filter(sex == sx)
  
  # Transform class column into a factor with the levels in the correct order:
  # FIRST CONTROL, SECOND TREATED.
  meta$class <- factor(meta$class, levels = c(control_label, classes))
  
  # Subset count matrix to keep samples of the chosen class and controls
  cnt <- as.matrix(exp_set[,colnames(exp_set) %in% meta$sample_name])
  
  # create dds object. There is no BD patient from cohort B
  dds_trat <- DESeqDataSetFromMatrix(countData = cnt,
                                     colData = meta,
                                     design = ~ age + interval + rin + class,
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
                                                      control_label)))
  res_annot <- res %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(group = sx)
  
  if(i==1){
    DEGs_all_DESeq2 <- res_annot
  }else{
    DEGs_all_DESeq2 <- rbind(DEGs_all_DESeq2, res_annot)
  }
  
}
rm(res_annot, res, meta, cnt, classes, control_label)

write.table(DEGs_all_DESeq2,
            file = 'intermediate/DEGs_BP_vs_Ctrl_OWM.tsv',
            sep = '\t',
            row.names = F,
            quote = F)

# Subset DEGs to keep only data from women and filter based on padj and log2FC
female_filtered <- subset.data.frame(DEGs_all_DESeq2, group == 0)
female_filtered <- female_filtered[female_filtered$log2FoldChange < -1.5 | 
                                     female_filtered$log2FoldChange > 1.5,]
female_filtered <- female_filtered[female_filtered$padj < 0.05,]

write.table(female_filtered,
            file = 'intermediate/female_filtered.tsv',
            sep = '\t',
            row.names = F,
            quote = F)

# Subset DEGs to keep only data from men and filter based on padj and log2FC
male_filtered <- subset.data.frame(DEGs_all_DESeq2, group == 1)
male_filtered <- male_filtered[male_filtered$log2FoldChange < -1.5 | 
                                 male_filtered$log2FoldChange > 1.5,]
male_filtered <- male_filtered[male_filtered$padj < 0.05,]
male_filtered <- male_filtered[!is.na(male_filtered$padj),]

write.table(male_filtered,
            file = 'intermediate/male_filtered.tsv',
            sep = '\t',
            row.names = F,
            quote = F)
