#### Transcriptome Study of Bipolar Disorder Type I and II ----

setwd('/home/leandro/Documents/IC/bpd_2020_krebs_phycomed/bpd_2020_krebs_phycomed')

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
                             "tobacco.use.ch1",
                             "assessment.group.ch1"
                           )]

colnames(samplesinfo)[which(names(samplesinfo) == 'title')] <- 'Sample'

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

# Filter rows with NA for specified columns so that batch effect can be removed
print(paste0("Number of rows with NA in sex column: ", sum(is.na(samplesinfo$Sex.ch1))))
library(tidyverse)
drop_indexes <- which(is.na(samplesinfo$Sex.ch1))
samplesinfo <- samplesinfo %>% drop_na(Sex.ch1)
exp_set <- subset(exp_set, select = -c(drop_indexes))
# Only keep samples from assessment group A
drop_indexes <- which(samplesinfo$assessment.group.ch1 == "B")
samplesinfo <- samplesinfo[samplesinfo$assessment.group.ch1 == "A",]
exp_set <- subset(exp_set, select = -c(drop_indexes))

# Too many samples have NA in tobacco use column, therefore they were not removed
print(paste0("Number of rows with NA in tobacco column: ", sum(is.na(samplesinfo$tobacco.use.ch1))))
# RIN, class and lithium columns do not have missing values
print(paste0("Number of rows with NA in RIN column: ", sum(is.na(samplesinfo$rin.ch1))))
print(paste0("Number of rows with NA in class column: ", sum(is.na(samplesinfo$Class))))
print(paste0("Number of rows with NA in lithium column: ", sum(is.na(samplesinfo$lithium.use..non.user.0..user...1..ch1))))

write.table(exp_set,
            file = 'intermediate/exp_set.tsv',
            sep = '\t',
            row.names = T)

write.table(samplesinfo,
            file = 'intermediate/samplesinfo.tsv',
            sep = '\t',
            row.names = T)






###################### Plot PCAs ##################################### 






setwd('/home/leandro/Documents/IC/bpd_2020_krebs_phycomed/bpd_2020_krebs_phycomed')
library(dplyr)
library(tibble)

exp_set = data.table::fread('intermediate/exp_set.tsv')
exp_set <- exp_set %>% remove_rownames %>% column_to_rownames(var="V1")

samplesinfo = data.table::fread('intermediate/samplesinfo.tsv')
samplesinfo <- samplesinfo %>% remove_rownames %>% column_to_rownames(var="V1")

# Define sample classes in a metatable
metatable <- data.frame(sample_name=samplesinfo$Sample, # Name of the samples
                        class=samplesinfo$Class,
                        tobacco_use=samplesinfo$tobacco.use.ch1,
                        lithium_use=samplesinfo$lithium.use..non.user.0..user...1..ch1,
                        sex=samplesinfo$Sex.ch1,
                        rin=samplesinfo$rin.ch1,
                        assessment = samplesinfo$assessment.group.ch1)
metatable <- metatable %>% remove_rownames %>% column_to_rownames(var="sample_name")
metatable <- data.frame(metatable,
                        sample_name=samplesinfo$Sample)


# Normalization with DESeq
library(DESeq2)
metatable$sex <- factor(metatable$sex)
metatable$lithium_use <- factor(metatable$lithium_use)
dds <- DESeqDataSetFromMatrix(countData = exp_set, 
                              colData = metatable, 
                              design = ~ lithium_use + sex + class)

library(FactoMineR)
library(limma)

# Plot PCA based on diagnostic
vsd <- vst(dds, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd),
                                       batch = vsd$rin,
                                       batch2 = vsd$sex,
                                       batch3 = vsd$tobacco_use,
                                       batch4 = vsd$lithium_use)

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
###################### BP1 and BP2 vs Ctrl ####

classes <- unique(samplesinfo$Class)[-c(3)]

# Define the control group label
control_label <- "Control"

# Run DESeq for all classes (In this case, BP1 or BP2 vs Control)
for(i in 1:length(classes)){
  
  # Define class to be compared in DE analysis
  cl <- classes[i]
  
  # Subset metadata to include only samples of the chosen class and controls
  meta <- metatable %>%
    filter(class %in% c(cl, control_label))
  
  # Transform class column into a factor with the levels in the correct order:
  # FIRST CONTROL, SECOND TREATED.
  meta$class <- factor(meta$class, levels = c(control_label, cl))
  
  # Subset count matrix to keep samples of the chosen class and controls
  cnt <- as.matrix(exp_set[,colnames(exp_set) %in% meta$sample_name])
  
  # create dds object
  meta$sex <- factor(meta$sex)
  dds <- DESeqDataSetFromMatrix(countData = cnt,
                                colData = meta,
                                design = ~ sex + class,
                                tidy = F)
  
  # Remove genes that have zero counts in all samples (i.e. not expressed genes)
  keep <- rowSums(counts(dds)) > 1
  dds <- dds[keep,]
  
  # Run DE analysis with default DESeq2 settings
  dds <- DESeq(dds)
  
  
  # extract results defining the order of comparison
  # must be: constrast=c("CLASS COLUMN IN META","DISEASE/TREATMENT GROUP","CONTROL)
  # This will give log2FoldChange values that correspond to disease vs control
  res <- as.data.frame(results(dds, contrast = c("class",
                                                 cl,
                                                 control_label)))
  
  res_annot <- res %>%
    tibble::rownames_to_column("Gene") %>%
    # left_join(metatable$class, by = "Gene") %>%
    mutate(group = cl)
  
  # Run lfcShrink to reduce the amplitude of the log2FoldChange values and remove
  # the log2FC to pvalue dependence/bias. We decided to keep this chunk commented
  # to avoid confusions in the future.
  # res_shrink <- lfcShrink(dds = dds,coef = "class_covid_moderate_vs_healthy",
  #                         type = "apeglm")
  # 
  # res_shrink_annot <- as.data.frame(res_shrink) %>%
  #   rownames_to_column("Gene") %>%
  #   left_join(id.dictionary,by="Gene") %>%
  #   mutate(group=cl)
  
  if(i==1){
    DEGs_all_DESeq2 <- res_annot
  }else{
    DEGs_all_DESeq2 <- rbind(DEGs_all_DESeq2, res_annot)
  }
  
}
rm(res_annot, res, meta, cnt, cl, classes, control_label)





###################### BP2 vs BP1 #### 

# Run DESeq for TEST classes (In this case, BP2 vs BP1)

# Redefine sample classes in a metatable
metatable <- metatable[metatable$class != 'Control', ]

# Define the targets/conditions (all classes but control/healthy)
classes <- "BP2"

# Define the control group label
control_label <- "BP1"

# Subset metadata to include only samples of the chosen class and controls
meta <- metatable %>%
  filter(class %in% c(classes, control_label))

# Transform class column into a factor with the levels in the correct order:
# FIRST CONTROL, SECOND TREATED.
meta$class <- factor(meta$class, levels = c(control_label, classes))

# Subset count matrix to keep samples of the chosen class and controls
cnt <- as.matrix(exp_set[,colnames(exp_set) %in% meta$sample_name])

# create dds object. There is no BD patient from cohort B
meta$sex <- factor(meta$sex)
meta$lithium_use <- factor(meta$lithium_use)
dds_trat <- DESeqDataSetFromMatrix(countData = cnt,
                              colData = meta,
                              design = ~ sex + lithium_use + class,
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

# DEGs of BP2 vs BP1
DEGs_test_groups = res

# Save DEGs as data frames 
dir.create('intermediate')
write.table(DEGs_all_DESeq2,
            file = 'intermediate/DEGs_BP1-2_vs_Control.tsv',
            sep = '\t',
            row.names = F,
            quote = F)

DEGs_test_groups$Gene <- rownames(DEGs_test_groups)
write.table(DEGs_test_groups,
            file = 'intermediate/DEGs_BP2_vs_BP1.tsv',
            sep = '\t',
            row.names = F,
            quote = F)




################################ BP2 vs BP1 (only women/men) ####
# Run DESeq for TEST classes (In this case, BP2 vs BP1)

# Redefine sample classes in a metatable
metatable <- metatable[metatable$class != 'Control', ]
sexes <- unique(metatable$sex)

# Define the targets/conditions (all classes but control/healthy)
classes <- "BP2"

# Define the control group label
control_label <- "BP1"

for(i in 1:length(sexes)){
  # Define class to be compared in DE analysis
  sx <- sexes[i]
  
  # Subset metadata to include only samples of the chosen class and controls
  meta <- metatable %>%
    filter(sex == sx)
  
  # Transform class column into a factor with the levels in the correct order:
  # FIRST CONTROL, SECOND TREATED.
  meta$class <- factor(meta$class, levels = c(control_label, classes))
  

  # Subset metadata to include only samples of the chosen class and controls
  meta <- metatable %>%
    filter(class %in% c(classes, control_label))

  # Transform class column into a factor with the levels in the correct order:
  # FIRST CONTROL, SECOND TREATED.
  meta$class <- factor(meta$class, levels = c(control_label, classes))

  # Subset count matrix to keep samples of the chosen class and controls
  cnt <- as.matrix(exp_set[,colnames(exp_set) %in% meta$sample_name])

  # create dds object. There is no BD patient from cohort B
  meta$lithium_use <- factor(meta$lithium_use)
  dds_trat <- DESeqDataSetFromMatrix(countData = cnt,
                                     colData = meta,
                                     design = ~ lithium_use + class,
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

# DEGs of BP2 vs BP1
DEGs_test_groups = DEGs_all_DESeq2

DEGs_test_groups$Gene <- rownames(DEGs_test_groups)
write.table(DEGs_test_groups,
            file = 'intermediate/DEGs_BP2_vs_BP1_OWM.tsv',
            sep = '\t',
            row.names = F,
            quote = F)




############################## BP vs Ctrl (women and men) ####
# Run DESeq for TEST classes (In this case, BP vs Control)

# Redefine sample classes in a metatable
metatable <- data.frame(sample_name=samplesinfo$Sample, # Name of the samples
                                                 class=samplesinfo$Class,
                                                 tobacco_use=samplesinfo$tobacco.use.ch1,
                                                 lithium_use=samplesinfo$lithium.use..non.user.0..user...1..ch1,
                                                 sex=samplesinfo$Sex.ch1,
                                                 rin=samplesinfo$rin.ch1,
                                                 assessment = samplesinfo$assessment.group.ch1)
metatable <- metatable %>% remove_rownames %>% column_to_rownames(var="sample_name")
metatable <- data.frame(metatable,
                        sample_name=samplesinfo$Sample)

# Define the targets/conditions (all classes but control/healthy)
metatable$class[metatable$class == "BP1" | metatable$class == "BP2"] <- "BP"
cl <- "BP"

# Define the control group label
control_label <- "Control"

sexes <- unique(metatable$sex)

# Run DESeq for all classes (In this case, BP1 or BP2 vs Control)
for(i in 1:length(sexes)){
  
  # Define class to be compared in DE analysis
  sx <- sexes[i]
  
  # Subset metadata to include only samples of the chosen class and controls
  meta <- metatable %>%
    filter(sex == sx)
  
  # Transform class column into a factor with the levels in the correct order:
  # FIRST CONTROL, SECOND TREATED.
  meta$class <- factor(meta$class, levels = c(control_label, cl))
  
  # Subset count matrix to keep samples of the chosen class and controls
  cnt <- as.matrix(exp_set[,colnames(exp_set) %in% meta$sample_name])
  
  # create dds object
  dds <- DESeqDataSetFromMatrix(countData = cnt,
                                colData = meta,
                                design = ~ class,
                                tidy = F)
  
  # Remove genes that have zero counts in all samples (i.e. not expressed genes)
  keep <- rowSums(counts(dds)) > 1
  dds <- dds[keep,]
  
  # Run DE analysis with default DESeq2 settings
  dds <- DESeq(dds)
  
  
  # extract results defining the order of comparison
  # must be: constrast=c("CLASS COLUMN IN META","DISEASE/TREATMENT GROUP","CONTROL)
  # This will give log2FoldChange values that correspond to disease vs control
  res <- as.data.frame(results(dds, contrast = c("class", cl, control_label)))
  
  res_annot <- res %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(group = sx)
  
  if(i==1){
    DEGs_all_DESeq2 <- res_annot
  }else{
    DEGs_all_DESeq2 <- rbind(DEGs_all_DESeq2, res_annot)
  }
  
}
rm(res_annot, res, meta, cnt, cl, classes, control_label)

DEGs_all_DESeq2$group[DEGs_all_DESeq2$group == 0] <- "F"
DEGs_all_DESeq2$group[DEGs_all_DESeq2$group == 1] <- "M"

write.table(DEGs_all_DESeq2,
            file = 'intermediate/DEGs_BP_vs_Ctrl.tsv',
            sep = '\t',
            row.names = F,
            quote = F)