#### Transcriptome Study of Bipolar Disorder Sexual Dimorphism ----
## Author: Leandro Tiburske
## Date: 23/07/2022

setwd('~/Documents/IC/Kohen_2014_TranslPsychiatry')

###################### Loading ###################### 
library(data.table)
exp_set = data.table::fread('~/Documents/IC/Kohen_2014_TranslPsychiatry/GSE42546/GSE42546_CleanedRawCounts.txt')

###################### Collapsing genes by  sum ###################### 
# Checking for duplicated genes
n_occur <- data.frame(table(exp_set$gene))
n_occur[n_occur$Freq > 1,]

# Sum duplicated genes to normalize the expression set
exp_set = aggregate(. ~ GeneSymbol, data=exp_set, FUN=sum)

# define genes as row names and delete the gene_symbol column
row.names(exp_set) = exp_set$GeneSymbol; exp_set$GeneSymbol = NULL

# Clean environment
rm(ensembl_genes, ensembl, gene_symbols, n_occur)

###################### Setup Meta and Expression ###################### 
samplesinfo <- read.table('~/Documents/IC/Kohen_2014_TranslPsychiatry/data/GSE42546/table_phenodata.tsv', header = T)
# Identify the sample names in metadata
head(colnames(exp_set))
head(table(samplesinfo$title))

# Filter samples info
names(samplesinfo)
samplesinfo <- samplesinfo[,
                           c("title",
                             "age.at.death.ch1",
                             "gender.ch1",
                             "diagnosis.ch1"
                           )]

colnames(samplesinfo)[which(names(samplesinfo) == 'title')] <- 'Sample'

# Subset samplesinfo object to keep only samples that are 
# present in the expression matrix
print(paste("Total number of samples:",
            nrow(samplesinfo[samplesinfo$Sample %in% colnames(exp_set),])))
samplesinfo = samplesinfo[samplesinfo$Sample %in% colnames(exp_set),]

# Check if all are in the same order
exp_set <- exp_set[, samplesinfo$Sample]
all(samplesinfo$Sample %in% colnames(exp_set))

# Set the classes you need to compare
samplesinfo$Class = samplesinfo$diagnosis.ch1
samplesinfo <- subset.data.frame(samplesinfo,
                                       !Class %in% c("depression", "schizophrenia"))
exp_set <- exp_set[, samplesinfo$Sample]

write.table(exp_set,
            file = 'intermediate/exp_set.tsv',
            sep = '\t',
            row.names = T)

write.table(samplesinfo,
            file = 'intermediate/samplesinfo.tsv',
            sep = '\t',
            row.names = T)


###################### DEG analysis ##################################### 

setwd('~/Documents/IC/Kohen_2014_TranslPsychiatry')
library(dplyr)
library(tibble)

samplesinfo = data.table::fread('intermediate/samplesinfo.tsv')
samplesinfo <- samplesinfo %>% remove_rownames %>% column_to_rownames(var="V1")

# Define sample classes in a metatable
metatable <- data.frame(sample_name=samplesinfo$Sample, # Name of the samples
                        class=samplesinfo$Class,
                        sex=samplesinfo$gender.ch1,
                        age=samplesinfo$age.at.death.ch1)
metatable <- metatable %>% remove_rownames %>% column_to_rownames(var="sample_name")
metatable <- data.frame(metatable,
                        sample_name=samplesinfo$Sample)


# Differentially expressed genes with edgeR
library(edgeR)

x <- read.delim("intermediate/exp_set.tsv")

## Obtain DEGs between bipolar and control women ####
### Subset data frames to keep data only from women

meta_women <- metatable[metatable$sex == "F",]
select <- meta_women$sample_name
exp_women <- x[select]
exp_women <- exp_women[, meta_women$sample_name]

### Create design variable to analyse gene expression based on diagnostic
### while considering possible batch effects due to age differences

design <- model.matrix(~ meta_women$age + meta_women$class)
DGE <- DGEList(counts = exp_women,
               samples = meta_women,
               group = meta_women$class)

### Filter genes based on counts, group sample sizes and library sizes
keep <- filterByExpr(DGE, design)
DGE_list <- DGE[keep,,keep.lib.sizes = FALSE]

### Normalize data
DGE <- calcNormFactors(DGE_list)

### Explore multidimensional scaling sample relations
plotMDS(DGE, method="bcv", col=as.numeric(DGE$samples$group))

### Estimate dispersion
y <- estimateDisp(DGE)

### Differential expression
et <- exactTest(y,
                pair = c("control", "bipolar"))
topTags(et, n = 10)


de <- decideTestsDGE(et,
                     adjust.method = "BH",
                     p.value = 0.05,
                     lfc = 2)
summary(de)
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags = detags)
abline(h = c(-2, 2), col = "blue")

DEGs <- et$table[as.logical(de),]

DGE_counts <- DGE_list$counts[rownames(DGE_list$counts) %in% rownames(DEGs), ]
DGE_counts <- cpm(DGE_counts, log = TRUE)

library(pheatmap)
pheatmap(DGE_counts,
         color = greenred(75),
         cluster_cols = FALSE)

### Write file
write.table(DEGs %>% rownames_to_column("Genes"),
            file = 'intermediate/DEGs_BP_vs_Ctrl_OW.tsv',
            sep = '\t',
            row.names = F)

## Obtain DEGs between bipolar and control men ####
### Subset data frames to keep data only from men

meta_men <- metatable[metatable$sex == "M",]
select <- meta_men$sample_name
exp_men <- x[select]

### Create design variable to analyse gene expression based on diagnostic
### while considering possible batch effects due to age differences

design <- model.matrix(~ meta_men$age + meta_men$class)
DGE <- DGEList(counts = exp_men,
               samples = meta_men,
               group = meta_men$class)

### Filter genes based on counts, group sample sizes and library sizes
keep <- filterByExpr(DGE, design)
DGE_list <- DGE[keep,,keep.lib.sizes = FALSE]

### Normalize data
DGE <- calcNormFactors(DGE_list)

### Estimate dispersion
y <- estimateDisp(DGE)


### Differential expression
et <- exactTest(y,
                pair = c("control", "bipolar"))
topTags(et, n = 10)

de <- decideTestsDGE(et,
                     adjust.method = "BH",
                     p.value = 0.05,
                     lfc = 2)
summary(de)
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags = detags)
abline(h = c(-2, 2), col = "blue")

DEGs <- et$table[as.logical(de),]
DGE_counts <- DGE$counts[rownames(DGE$counts) %in% rownames(DEGs), ]
DGE_counts <- cpm(DGE_counts, log = TRUE)

library(pheatmap)
pheatmap(DGE_counts,
         color = greenred(75),
         cluster_cols = FALSE)

### Write file
write.table(DEGs %>% rownames_to_column("Genes"),
            file = 'intermediate/DEGs_BP_vs_Ctrl_OM.tsv',
            sep = '\t',
            row.names = F)
