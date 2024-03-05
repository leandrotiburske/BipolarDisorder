library(Seurat) # Seurat v5


################################### Batch 1 #####################################

D19_4295 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-4295/outs/filtered_feature_bc_matrix/")

D19_4303 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-4303_S3/umi_count/",
                    gene.column = 1) 

colnames(D19_4303) <- paste0(colnames(D19_4303), "-1")

gc()

joint.bcs <- intersect(colnames(D19_4295), colnames(D19_4303))

D19_4295 <- D19_4295[, joint.bcs]
D19_4303 <- as.matrix(D19_4303[, joint.bcs])

# Setup Seurat object
D19_4295 <- CreateSeuratObject(counts = D19_4295)
gc()

# Normalize RNA data with log normalization
D19_4295 <- NormalizeData(D19_4295)
# Find and scale variable features
D19_4295 <- FindVariableFeatures(D19_4295, selection.method = "mean.var.plot")
D19_4295 <- ScaleData(D19_4295, features = VariableFeatures(D19_4295))

gc()

# Add HTO data as a new assay independent from RNA
D19_4295[["HTO"]] <- CreateAssayObject(counts = D19_4303)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_4295 <- NormalizeData(D19_4295, assay = "HTO", normalization.method = "CLR")

D19_4295 <- HTODemux(D19_4295, assay = "HTO", positive.quantile = 0.99)

table(D19_4295$HTO_classification.global)

Idents(D19_4295) <- "HTO_classification.global"

D19_4295 <- subset(D19_4295, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_4303, joint.bcs); gc()

Idents(D19_4295) <- "hash.ID"

D19_4295 <- subset(D19_4295, 
                   idents = c("SZ2-ATGATGAA", "SZ1-GGTAGATG", "SZ3-CTCGAACG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_4295$HTO_classification)



D19_4296 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-4296/outs/filtered_feature_bc_matrix/")

D19_4304 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-4304_S4/umi_count/",
                    gene.column = 1) 

colnames(D19_4304) <- paste0(colnames(D19_4304), "-1")

gc()

joint.bcs <- intersect(colnames(D19_4296), colnames(D19_4304))

D19_4296 <- D19_4296[, joint.bcs]
D19_4304 <- as.matrix(D19_4304[, joint.bcs])

# Setup Seurat object
D19_4296 <- CreateSeuratObject(counts = D19_4296)
gc()

# Normalize RNA data with log normalization
D19_4296 <- NormalizeData(D19_4296)
# Find and scale variable features
D19_4296 <- FindVariableFeatures(D19_4296, selection.method = "mean.var.plot")
D19_4296 <- ScaleData(D19_4296, features = VariableFeatures(D19_4296))

gc()

# Add HTO data as a new assay independent from RNA
D19_4296[["HTO"]] <- CreateAssayObject(counts = D19_4304)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_4296 <- NormalizeData(D19_4296, assay = "HTO", normalization.method = "CLR")

D19_4296 <- HTODemux(D19_4296, assay = "HTO", positive.quantile = 0.99)

table(D19_4296$HTO_classification.global)

Idents(D19_4296) <- "HTO_classification.global"

D19_4296 <- subset(D19_4296, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_4304, joint.bcs); gc()

table(D19_4296$HTO_classification)

Idents(D19_4296) <- "hash.ID"

D19_4296 <- subset(D19_4296, 
                   idents = c("SZ2-ATGATGAA", "SZ1-GGTAGATG", "SZ3-CTCGAACG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_4296$HTO_classification)



D19_4297 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-4297/outs/filtered_feature_bc_matrix/")

D19_4305 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-4305_S3/umi_count/",
                    gene.column = 1) 

colnames(D19_4305) <- paste0(colnames(D19_4305), "-1")

gc()

joint.bcs <- intersect(colnames(D19_4297), colnames(D19_4305))

D19_4297 <- D19_4297[, joint.bcs]
D19_4305 <- as.matrix(D19_4305[, joint.bcs])

# Setup Seurat object
D19_4297 <- CreateSeuratObject(counts = D19_4297)
gc()

# Normalize RNA data with log normalization
D19_4297 <- NormalizeData(D19_4297)
# Find and scale variable features
D19_4297 <- FindVariableFeatures(D19_4297, selection.method = "mean.var.plot")
D19_4297 <- ScaleData(D19_4297, features = VariableFeatures(D19_4297))

gc()

# Add HTO data as a new assay independent from RNA
D19_4297[["HTO"]] <- CreateAssayObject(counts = D19_4305)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_4297 <- NormalizeData(D19_4297, assay = "HTO", normalization.method = "CLR")

D19_4297 <- HTODemux(D19_4297, assay = "HTO", positive.quantile = 0.99)

table(D19_4297$HTO_classification.global)

Idents(D19_4297) <- "HTO_classification.global"

D19_4297 <- subset(D19_4297, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_4305, joint.bcs); gc()

table(D19_4297$HTO_classification)

Idents(D19_4297) <- "hash.ID"

D19_4297 <- subset(D19_4297, 
                   idents = c("SZ2-ATGATGAA", "SZ1-GGTAGATG", "SZ3-CTCGAACG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_4297$HTO_classification)



D19_4298 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-4298/outs/filtered_feature_bc_matrix/")

D19_4306 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-4306_S4/umi_count/",
                    gene.column = 1) 

colnames(D19_4306) <- paste0(colnames(D19_4306), "-1")

gc()

joint.bcs <- intersect(colnames(D19_4298), colnames(D19_4306))

D19_4298 <- D19_4298[, joint.bcs]
D19_4306 <- as.matrix(D19_4306[, joint.bcs])

# Setup Seurat object
D19_4298 <- CreateSeuratObject(counts = D19_4298)
gc()

# Normalize RNA data with log normalization
D19_4298 <- NormalizeData(D19_4298)
# Find and scale variable features
D19_4298 <- FindVariableFeatures(D19_4298, selection.method = "mean.var.plot")
D19_4298 <- ScaleData(D19_4298, features = VariableFeatures(D19_4298))

gc()

# Add HTO data as a new assay independent from RNA
D19_4298[["HTO"]] <- CreateAssayObject(counts = D19_4306)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_4298 <- NormalizeData(D19_4298, assay = "HTO", normalization.method = "CLR")

D19_4298 <- HTODemux(D19_4298, assay = "HTO", positive.quantile = 0.99)

table(D19_4298$HTO_classification.global)

Idents(D19_4298) <- "HTO_classification.global"

D19_4298 <- subset(D19_4298, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_4306, joint.bcs); gc()

table(D19_4298$HTO_classification)

Idents(D19_4298) <- "hash.ID"

D19_4298 <- subset(D19_4298, 
                   idents = c("SZ2-ATGATGAA", "SZ1-GGTAGATG", "SZ3-CTCGAACG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_4298$HTO_classification)




D19_4299 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-4299/outs/filtered_feature_bc_matrix/")

D19_4307 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-4307_S3/umi_count/",
                    gene.column = 1) 

colnames(D19_4307) <- paste0(colnames(D19_4307), "-1")

gc()

joint.bcs <- intersect(colnames(D19_4299), colnames(D19_4307))

D19_4299 <- D19_4299[, joint.bcs]
D19_4307 <- as.matrix(D19_4307[, joint.bcs])

# Setup Seurat object
D19_4299 <- CreateSeuratObject(counts = D19_4299)
gc()

# Normalize RNA data with log normalization
D19_4299 <- NormalizeData(D19_4299)
# Find and scale variable features
D19_4299 <- FindVariableFeatures(D19_4299, selection.method = "mean.var.plot")
D19_4299 <- ScaleData(D19_4299, features = VariableFeatures(D19_4299))

gc()

# Add HTO data as a new assay independent from RNA
D19_4299[["HTO"]] <- CreateAssayObject(counts = D19_4307)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_4299 <- NormalizeData(D19_4299, assay = "HTO", normalization.method = "CLR")

D19_4299 <- HTODemux(D19_4299, assay = "HTO", positive.quantile = 0.99)

table(D19_4299$HTO_classification.global)

Idents(D19_4299) <- "HTO_classification.global"

D19_4299 <- subset(D19_4299, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_4307, joint.bcs); gc()

table(D19_4299$HTO_classification)

Idents(D19_4299) <- "hash.ID"

D19_4299 <- subset(D19_4299, 
                   idents = c("SZ2-ATGATGAA", "SZ1-GGTAGATG", "SZ3-CTCGAACG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_4299$HTO_classification)




D19_4300 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-4300/outs/filtered_feature_bc_matrix/")

D19_4308 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-4308_S4/umi_count/",
                    gene.column = 1) 

colnames(D19_4308) <- paste0(colnames(D19_4308), "-1")

gc()

joint.bcs <- intersect(colnames(D19_4300), colnames(D19_4308))

D19_4300 <- D19_4300[, joint.bcs]
D19_4308 <- as.matrix(D19_4308[, joint.bcs])

# Setup Seurat object
D19_4300 <- CreateSeuratObject(counts = D19_4300)
gc()

# Normalize RNA data with log normalization
D19_4300 <- NormalizeData(D19_4300)
# Find and scale variable features
D19_4300 <- FindVariableFeatures(D19_4300, selection.method = "mean.var.plot")
D19_4300 <- ScaleData(D19_4300, features = VariableFeatures(D19_4300))

gc()

# Add HTO data as a new assay independent from RNA
D19_4300[["HTO"]] <- CreateAssayObject(counts = D19_4308)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_4300 <- NormalizeData(D19_4300, assay = "HTO", normalization.method = "CLR")

D19_4300 <- HTODemux(D19_4300, assay = "HTO", positive.quantile = 0.99)

table(D19_4300$HTO_classification.global)

Idents(D19_4300) <- "HTO_classification.global"

D19_4300 <- subset(D19_4300, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_4308, joint.bcs); gc()

table(D19_4300$HTO_classification)

Idents(D19_4300) <- "hash.ID"

D19_4300 <- subset(D19_4300, 
                   idents = c("SZ2-ATGATGAA", "SZ1-GGTAGATG", "SZ3-CTCGAACG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_4300$HTO_classification)



D19_4301 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-4301/outs/filtered_feature_bc_matrix/")

D19_4309 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-4309_S3/umi_count/",
                    gene.column = 1) 

colnames(D19_4309) <- paste0(colnames(D19_4309), "-1")

gc()

joint.bcs <- intersect(colnames(D19_4301), colnames(D19_4309))

D19_4301 <- D19_4301[, joint.bcs]
D19_4309 <- as.matrix(D19_4309[, joint.bcs])

# Setup Seurat object
D19_4301 <- CreateSeuratObject(counts = D19_4301)
gc()

# Normalize RNA data with log normalization
D19_4301 <- NormalizeData(D19_4301)
# Find and scale variable features
D19_4301 <- FindVariableFeatures(D19_4301, selection.method = "mean.var.plot")
D19_4301 <- ScaleData(D19_4301, features = VariableFeatures(D19_4301))

gc()

# Add HTO data as a new assay independent from RNA
D19_4301[["HTO"]] <- CreateAssayObject(counts = D19_4309)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_4301 <- NormalizeData(D19_4301, assay = "HTO", normalization.method = "CLR")

D19_4301 <- HTODemux(D19_4301, assay = "HTO", positive.quantile = 0.99)

table(D19_4301$HTO_classification.global)

Idents(D19_4301) <- "HTO_classification.global"

D19_4301 <- subset(D19_4301, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_4309, joint.bcs); gc()

table(D19_4301$HTO_classification)

Idents(D19_4301) <- "hash.ID"

D19_4301 <- subset(D19_4301, 
                   idents = c("SZ2-ATGATGAA", "SZ1-GGTAGATG", "SZ3-CTCGAACG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_4301$HTO_classification)



D19_4302 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-4302/outs/filtered_feature_bc_matrix/")

D19_4310 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-4310_S4/umi_count/",
                    gene.column = 1) 

colnames(D19_4310) <- paste0(colnames(D19_4310), "-1")

gc()

joint.bcs <- intersect(colnames(D19_4302), colnames(D19_4310))

D19_4302 <- D19_4302[, joint.bcs]
D19_4310 <- as.matrix(D19_4310[, joint.bcs])

# Setup Seurat object
D19_4302 <- CreateSeuratObject(counts = D19_4302)
gc()

# Normalize RNA data with log normalization
D19_4302 <- NormalizeData(D19_4302)
# Find and scale variable features
D19_4302 <- FindVariableFeatures(D19_4302, selection.method = "mean.var.plot")
D19_4302 <- ScaleData(D19_4302, features = VariableFeatures(D19_4302))

gc()

# Add HTO data as a new assay independent from RNA
D19_4302[["HTO"]] <- CreateAssayObject(counts = D19_4310)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_4302 <- NormalizeData(D19_4302, assay = "HTO", normalization.method = "CLR")

D19_4302 <- HTODemux(D19_4302, assay = "HTO", positive.quantile = 0.99)

table(D19_4302$HTO_classification.global)

Idents(D19_4302) <- "HTO_classification.global"

D19_4302 <- subset(D19_4302, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_4310, joint.bcs); gc()

table(D19_4302$HTO_classification)

Idents(D19_4302) <- "hash.ID"

D19_4302 <- subset(D19_4302, 
                   idents = c("SZ2-ATGATGAA", "SZ1-GGTAGATG", "SZ3-CTCGAACG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_4302$HTO_classification)



# merge all batch 1 seurat objects


batch1 <- merge(D19_4295, y = c(D19_4296, D19_4297, D19_4298,
                                D19_4299, D19_4300, D19_4301, D19_4302),
                project = "Batch1")

rm(D19_4295, D19_4296, D19_4297, D19_4298,
   D19_4299, D19_4300, D19_4301, D19_4302); gc()

batch1[["percent.mt"]] <- PercentageFeatureSet(batch1, pattern = "^MT-"); gc()
# Visualize QC metrics as a violin plot
VlnPlot(batch1,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,group.by = "orig.ident")

batch1 <- subset(batch1, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 15); gc()


saveRDS(batch1, "/media/leandro/SAMSUNG/FASTQs/batch1.rds")


################################### Batch 2 #####################################

D19_5859 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-5859/outs/filtered_feature_bc_matrix/")

D19_5867 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-5867_S3/umi_count/",
                    gene.column = 1) 

colnames(D19_5867) <- paste0(colnames(D19_5867), "-1")

gc()

joint.bcs <- intersect(colnames(D19_5859), colnames(D19_5867))

D19_5859 <- D19_5859[, joint.bcs]
D19_5867 <- as.matrix(D19_5867[, joint.bcs])

# Setup Seurat object
D19_5859 <- CreateSeuratObject(counts = D19_5859)
gc()

# Normalize RNA data with log normalization
D19_5859 <- NormalizeData(D19_5859)
# Find and scale variable features
D19_5859 <- FindVariableFeatures(D19_5859, selection.method = "mean.var.plot")
D19_5859 <- ScaleData(D19_5859, features = VariableFeatures(D19_5859))

gc()

# Add HTO data as a new assay independent from RNA
D19_5859[["HTO"]] <- CreateAssayObject(counts = D19_5867)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_5859 <- NormalizeData(D19_5859, assay = "HTO", normalization.method = "CLR")

D19_5859 <- HTODemux(D19_5859, assay = "HTO", positive.quantile = 0.99)

table(D19_5859$HTO_classification.global)

Idents(D19_5859) <- "HTO_classification.global"

D19_5859 <- subset(D19_5859, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_5867, joint.bcs); gc()

table(D19_5859$HTO_classification)

Idents(D19_5859) <- "hash.ID"

D19_5859 <- subset(D19_5859, 
                   idents = c("SZ4-TTCCTGCC", "SZ5-CCGTACCT", "SZ6-TGACGCCG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_5859$HTO_classification)



D19_5860 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-5860/outs/filtered_feature_bc_matrix/")

D19_5868 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-5868_S4/umi_count/",
                    gene.column = 1) 

colnames(D19_5868) <- paste0(colnames(D19_5868), "-1")

gc()

joint.bcs <- intersect(colnames(D19_5860), colnames(D19_5868))

D19_5860 <- D19_5860[, joint.bcs]
D19_5868 <- as.matrix(D19_5868[, joint.bcs])

# Setup Seurat object
D19_5860 <- CreateSeuratObject(counts = D19_5860)
gc()

# Normalize RNA data with log normalization
D19_5860 <- NormalizeData(D19_5860)
# Find and scale variable features
D19_5860 <- FindVariableFeatures(D19_5860, selection.method = "mean.var.plot")
D19_5860 <- ScaleData(D19_5860, features = VariableFeatures(D19_5860))

gc()

# Add HTO data as a new assay independent from RNA
D19_5860[["HTO"]] <- CreateAssayObject(counts = D19_5868)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_5860 <- NormalizeData(D19_5860, assay = "HTO", normalization.method = "CLR")

D19_5860 <- HTODemux(D19_5860, assay = "HTO", positive.quantile = 0.99)

table(D19_5860$HTO_classification.global)

Idents(D19_5860) <- "HTO_classification.global"

D19_5860 <- subset(D19_5860, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_5868, joint.bcs); gc()

table(D19_5860$HTO_classification)

Idents(D19_5860) <- "hash.ID"

D19_5860 <- subset(D19_5860, 
                   idents = c("SZ4-TTCCTGCC", "SZ5-CCGTACCT", "SZ6-TGACGCCG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_5860$HTO_classification)




D19_5861 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-5861/outs/filtered_feature_bc_matrix/")

D19_5869 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-5869_S3/umi_count/",
                    gene.column = 1) 

colnames(D19_5869) <- paste0(colnames(D19_5869), "-1")

gc()

joint.bcs <- intersect(colnames(D19_5861), colnames(D19_5869))

D19_5861 <- D19_5861[, joint.bcs]
D19_5869 <- as.matrix(D19_5869[, joint.bcs])

# Setup Seurat object
D19_5861 <- CreateSeuratObject(counts = D19_5861)
gc()

# Normalize RNA data with log normalization
D19_5861 <- NormalizeData(D19_5861)
# Find and scale variable features
D19_5861 <- FindVariableFeatures(D19_5861, selection.method = "mean.var.plot")
D19_5861 <- ScaleData(D19_5861, features = VariableFeatures(D19_5861))

gc()

# Add HTO data as a new assay independent from RNA
D19_5861[["HTO"]] <- CreateAssayObject(counts = D19_5869)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_5861 <- NormalizeData(D19_5861, assay = "HTO", normalization.method = "CLR")

D19_5861 <- HTODemux(D19_5861, assay = "HTO", positive.quantile = 0.99)

table(D19_5861$HTO_classification.global)

Idents(D19_5861) <- "HTO_classification.global"

D19_5861 <- subset(D19_5861, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_5869, joint.bcs); gc()

table(D19_5861$HTO_classification)

Idents(D19_5861) <- "hash.ID"

D19_5861 <- subset(D19_5861, 
                   idents = c("SZ4-TTCCTGCC", "SZ5-CCGTACCT", "SZ6-TGACGCCG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_5861$HTO_classification)




D19_5862 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-5862/outs/filtered_feature_bc_matrix/")

D19_5870 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-5870_S4/umi_count/",
                    gene.column = 1) 

colnames(D19_5870) <- paste0(colnames(D19_5870), "-1")

gc()

joint.bcs <- intersect(colnames(D19_5862), colnames(D19_5870))

D19_5862 <- D19_5862[, joint.bcs]
D19_5870 <- as.matrix(D19_5870[, joint.bcs])

# Setup Seurat object
D19_5862 <- CreateSeuratObject(counts = D19_5862)
gc()

# Normalize RNA data with log normalization
D19_5862 <- NormalizeData(D19_5862)
# Find and scale variable features
D19_5862 <- FindVariableFeatures(D19_5862, selection.method = "mean.var.plot")
D19_5862 <- ScaleData(D19_5862, features = VariableFeatures(D19_5862))

gc()

# Add HTO data as a new assay independent from RNA
D19_5862[["HTO"]] <- CreateAssayObject(counts = D19_5870)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_5862 <- NormalizeData(D19_5862, assay = "HTO", normalization.method = "CLR")

D19_5862 <- HTODemux(D19_5862, assay = "HTO", positive.quantile = 0.99)

table(D19_5862$HTO_classification.global)

Idents(D19_5862) <- "HTO_classification.global"

D19_5862 <- subset(D19_5862, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_5870, joint.bcs); gc()

table(D19_5862$HTO_classification)

Idents(D19_5862) <- "hash.ID"

D19_5862 <- subset(D19_5862, 
                   idents = c("SZ4-TTCCTGCC", "SZ5-CCGTACCT", "SZ6-TGACGCCG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_5862$HTO_classification)



D19_5863 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-5863/outs/filtered_feature_bc_matrix/")

D19_5871 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-5871_S3/umi_count/",
                    gene.column = 1) 

colnames(D19_5871) <- paste0(colnames(D19_5871), "-1")

gc()

joint.bcs <- intersect(colnames(D19_5863), colnames(D19_5871))

D19_5863 <- D19_5863[, joint.bcs]
D19_5871 <- as.matrix(D19_5871[, joint.bcs])

# Setup Seurat object
D19_5863 <- CreateSeuratObject(counts = D19_5863)
gc()

# Normalize RNA data with log normalization
D19_5863 <- NormalizeData(D19_5863)
# Find and scale variable features
D19_5863 <- FindVariableFeatures(D19_5863, selection.method = "mean.var.plot")
D19_5863 <- ScaleData(D19_5863, features = VariableFeatures(D19_5863))

gc()

# Add HTO data as a new assay independent from RNA
D19_5863[["HTO"]] <- CreateAssayObject(counts = D19_5871)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_5863 <- NormalizeData(D19_5863, assay = "HTO", normalization.method = "CLR")

D19_5863 <- HTODemux(D19_5863, assay = "HTO", positive.quantile = 0.99)

table(D19_5863$HTO_classification.global)

Idents(D19_5863) <- "HTO_classification.global"

D19_5863 <- subset(D19_5863, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_5871, joint.bcs); gc()

table(D19_5863$HTO_classification)

Idents(D19_5863) <- "hash.ID"

D19_5863 <- subset(D19_5863, 
                   idents = c("SZ4-TTCCTGCC", "SZ5-CCGTACCT", "SZ6-TGACGCCG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_5863$HTO_classification)



D19_5864 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-5864/outs/filtered_feature_bc_matrix/")

D19_5872 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-5872_S4/umi_count/",
                    gene.column = 1) 

colnames(D19_5872) <- paste0(colnames(D19_5872), "-1")

gc()

joint.bcs <- intersect(colnames(D19_5864), colnames(D19_5872))

D19_5864 <- D19_5864[, joint.bcs]
D19_5872 <- as.matrix(D19_5872[, joint.bcs])

# Setup Seurat object
D19_5864 <- CreateSeuratObject(counts = D19_5864)
gc()

# Normalize RNA data with log normalization
D19_5864 <- NormalizeData(D19_5864)
# Find and scale variable features
D19_5864 <- FindVariableFeatures(D19_5864, selection.method = "mean.var.plot")
D19_5864 <- ScaleData(D19_5864, features = VariableFeatures(D19_5864))

gc()

# Add HTO data as a new assay independent from RNA
D19_5864[["HTO"]] <- CreateAssayObject(counts = D19_5872)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_5864 <- NormalizeData(D19_5864, assay = "HTO", normalization.method = "CLR")

D19_5864 <- HTODemux(D19_5864, assay = "HTO", positive.quantile = 0.99)

table(D19_5864$HTO_classification.global)

Idents(D19_5864) <- "HTO_classification.global"

D19_5864 <- subset(D19_5864, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_5872, joint.bcs); gc()

table(D19_5864$HTO_classification)

Idents(D19_5864) <- "hash.ID"

D19_5864 <- subset(D19_5864, 
                   idents = c("SZ4-TTCCTGCC", "SZ5-CCGTACCT", "SZ6-TGACGCCG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_5864$HTO_classification)



D19_5865 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-5865/outs/filtered_feature_bc_matrix/")

D19_6155 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-6155_S4/umi_count/",
                    gene.column = 1) 

colnames(D19_6155) <- paste0(colnames(D19_6155), "-1")

gc()

joint.bcs <- intersect(colnames(D19_5865), colnames(D19_6155))

D19_5865 <- D19_5865[, joint.bcs]
D19_6155 <- as.matrix(D19_6155[, joint.bcs])

# Setup Seurat object
D19_5865 <- CreateSeuratObject(counts = D19_5865)
gc()

# Normalize RNA data with log normalization
D19_5865 <- NormalizeData(D19_5865)
# Find and scale variable features
D19_5865 <- FindVariableFeatures(D19_5865, selection.method = "mean.var.plot")
D19_5865 <- ScaleData(D19_5865, features = VariableFeatures(D19_5865))

gc()

# Add HTO data as a new assay independent from RNA
D19_5865[["HTO"]] <- CreateAssayObject(counts = D19_6155)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_5865 <- NormalizeData(D19_5865, assay = "HTO", normalization.method = "CLR")

D19_5865 <- HTODemux(D19_5865, assay = "HTO", positive.quantile = 0.99)

table(D19_5865$HTO_classification.global)

Idents(D19_5865) <- "HTO_classification.global"

D19_5865 <- subset(D19_5865, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6155, joint.bcs); gc()

table(D19_5865$HTO_classification)

Idents(D19_5865) <- "hash.ID"

D19_5865 <- subset(D19_5865, 
                   idents = c("SZ4-TTCCTGCC", "SZ5-CCGTACCT", "SZ6-TGACGCCG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_5865$HTO_classification)



D19_5866 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-5866/outs/filtered_feature_bc_matrix/")

D19_5874 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/D19-5874_S3/umi_count/",
                    gene.column = 1) 

colnames(D19_5874) <- paste0(colnames(D19_5874), "-1")

gc()

joint.bcs <- intersect(colnames(D19_5866), colnames(D19_5874))

D19_5866 <- D19_5866[, joint.bcs]
D19_5874 <- as.matrix(D19_5874[, joint.bcs])

# Setup Seurat object
D19_5866 <- CreateSeuratObject(counts = D19_5866)
gc()

# Normalize RNA data with log normalization
D19_5866 <- NormalizeData(D19_5866)
# Find and scale variable features
D19_5866 <- FindVariableFeatures(D19_5866, selection.method = "mean.var.plot")
D19_5866 <- ScaleData(D19_5866, features = VariableFeatures(D19_5866))

gc()

# Add HTO data as a new assay independent from RNA
D19_5866[["HTO"]] <- CreateAssayObject(counts = D19_5874)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_5866 <- NormalizeData(D19_5866, assay = "HTO", normalization.method = "CLR")

D19_5866 <- HTODemux(D19_5866, assay = "HTO", positive.quantile = 0.99)

table(D19_5866$HTO_classification.global)

Idents(D19_5866) <- "HTO_classification.global"

D19_5866 <- subset(D19_5866, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_5874, joint.bcs); gc()

table(D19_5866$HTO_classification)

Idents(D19_5866) <- "hash.ID"

D19_5866 <- subset(D19_5866, 
                   idents = c("SZ4-TTCCTGCC", "SZ5-CCGTACCT", "SZ6-TGACGCCG",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_5866$HTO_classification)


# merge all batch 2 seurat objects


batch2 <- merge(D19_5859, y = c(D19_5860, D19_5861, D19_5862,
                                D19_5863, D19_5864, D19_5865, D19_5866),
                project = "Batch2")

rm(D19_5859, D19_5860, D19_5861, D19_5862,
   D19_5863, D19_5864, D19_5865, D19_5866); gc()

batch2[["percent.mt"]] <- PercentageFeatureSet(batch2, pattern = "^MT-"); gc()
# Visualize QC metrics as a violin plot
VlnPlot(batch2,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,group.by = "orig.ident")

batch2 <- subset(batch2, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 15); gc()


saveRDS(batch2, "/media/leandro/SAMSUNG/FASTQs/batch2.rds")


################################### Batch 3 #####################################

D19_6771 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6771/outs/filtered_feature_bc_matrix/")

D19_6787 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6787/umi_count/",
                    gene.column = 1) 

colnames(D19_6787) <- paste0(colnames(D19_6787), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6771), colnames(D19_6787))

D19_6771 <- D19_6771[, joint.bcs]
D19_6787 <- as.matrix(D19_6787[, joint.bcs])

# Setup Seurat object
D19_6771 <- CreateSeuratObject(counts = D19_6771)
gc()

# Normalize RNA data with log normalization
D19_6771 <- NormalizeData(D19_6771)
# Find and scale variable features
D19_6771 <- FindVariableFeatures(D19_6771, selection.method = "mean.var.plot")
D19_6771 <- ScaleData(D19_6771, features = VariableFeatures(D19_6771))

gc()

# Add HTO data as a new assay independent from RNA
D19_6771[["HTO"]] <- CreateAssayObject(counts = D19_6787)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6771 <- NormalizeData(D19_6771, assay = "HTO", normalization.method = "CLR")

D19_6771 <- HTODemux(D19_6771, assay = "HTO", positive.quantile = 0.99)

table(D19_6771$HTO_classification.global)

Idents(D19_6771) <- "HTO_classification.global"

D19_6771 <- subset(D19_6771, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6787, joint.bcs); gc()

Idents(D19_6771) <- "hash.ID"

D19_6771 <- subset(D19_6771, 
                   idents = c("SZ7-TTCCTGCC", "SZ8-CCGTACCT", "SZ9-GCCTAGTA"),
                   invert = TRUE); gc()


D19_6772 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6772/outs/filtered_feature_bc_matrix/")

D19_6788 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6788/umi_count/",
                    gene.column = 1) 

colnames(D19_6788) <- paste0(colnames(D19_6788), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6772), colnames(D19_6788))

D19_6772 <- D19_6772[, joint.bcs]
D19_6788 <- as.matrix(D19_6788[, joint.bcs])

# Setup Seurat object
D19_6772 <- CreateSeuratObject(counts = D19_6772)
gc()

# Normalize RNA data with log normalization
D19_6772 <- NormalizeData(D19_6772)
# Find and scale variable features
D19_6772 <- FindVariableFeatures(D19_6772, selection.method = "mean.var.plot")
D19_6772 <- ScaleData(D19_6772, features = VariableFeatures(D19_6772))

gc()

# Add HTO data as a new assay independent from RNA
D19_6772[["HTO"]] <- CreateAssayObject(counts = D19_6788)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6772 <- NormalizeData(D19_6772, assay = "HTO", normalization.method = "CLR")

D19_6772 <- HTODemux(D19_6772, assay = "HTO", positive.quantile = 0.99)

table(D19_6772$HTO_classification.global)

Idents(D19_6772) <- "HTO_classification.global"

D19_6772 <- subset(D19_6772, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6788, joint.bcs); gc()

Idents(D19_6772) <- "hash.ID"

D19_6772 <- subset(D19_6772, 
                   idents = c("SZ7-TTCCTGCC", "SZ8-CCGTACCT", "SZ9-GCCTAGTA"),
                   invert = TRUE); gc()


D19_6773 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6773/outs/filtered_feature_bc_matrix/")

D19_6789 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6789/umi_count/",
                    gene.column = 1) 

colnames(D19_6789) <- paste0(colnames(D19_6789), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6773), colnames(D19_6789))

D19_6773 <- D19_6773[, joint.bcs]
D19_6789 <- as.matrix(D19_6789[, joint.bcs])

# Setup Seurat object
D19_6773 <- CreateSeuratObject(counts = D19_6773)
gc()

# Normalize RNA data with log normalization
D19_6773 <- NormalizeData(D19_6773)
# Find and scale variable features
D19_6773 <- FindVariableFeatures(D19_6773, selection.method = "mean.var.plot")
D19_6773 <- ScaleData(D19_6773, features = VariableFeatures(D19_6773))

gc()

# Add HTO data as a new assay independent from RNA
D19_6773[["HTO"]] <- CreateAssayObject(counts = D19_6789)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6773 <- NormalizeData(D19_6773, assay = "HTO", normalization.method = "CLR")

D19_6773 <- HTODemux(D19_6773, assay = "HTO", positive.quantile = 0.99)

table(D19_6773$HTO_classification.global)

Idents(D19_6773) <- "HTO_classification.global"

D19_6773 <- subset(D19_6773, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6789, joint.bcs); gc()

Idents(D19_6773) <- "hash.ID"

D19_6773 <- subset(D19_6773, 
                   idents = c("SZ7-TTCCTGCC", "SZ8-CCGTACCT", "SZ9-GCCTAGTA"),
                   invert = TRUE); gc()


D19_6774 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6774/outs/filtered_feature_bc_matrix/")

D19_6790 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6790/umi_count/",
                    gene.column = 1) 

colnames(D19_6790) <- paste0(colnames(D19_6790), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6774), colnames(D19_6790))

D19_6774 <- D19_6774[, joint.bcs]
D19_6790 <- as.matrix(D19_6790[, joint.bcs])

# Setup Seurat object
D19_6774 <- CreateSeuratObject(counts = D19_6774)
gc()

# Normalize RNA data with log normalization
D19_6774 <- NormalizeData(D19_6774)
# Find and scale variable features
D19_6774 <- FindVariableFeatures(D19_6774, selection.method = "mean.var.plot")
D19_6774 <- ScaleData(D19_6774, features = VariableFeatures(D19_6774))

gc()

# Add HTO data as a new assay independent from RNA
D19_6774[["HTO"]] <- CreateAssayObject(counts = D19_6790)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6774 <- NormalizeData(D19_6774, assay = "HTO", normalization.method = "CLR")

D19_6774 <- HTODemux(D19_6774, assay = "HTO", positive.quantile = 0.99)

table(D19_6774$HTO_classification.global)

Idents(D19_6774) <- "HTO_classification.global"

D19_6774 <- subset(D19_6774, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6790, joint.bcs); gc()

Idents(D19_6774) <- "hash.ID"

D19_6774 <- subset(D19_6774, 
                   idents = c("SZ7-TTCCTGCC", "SZ8-CCGTACCT", "SZ9-GCCTAGTA"),
                   invert = TRUE); gc()


D19_6775 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6775/outs/filtered_feature_bc_matrix/")

D19_6791 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6791/umi_count/",
                    gene.column = 1) 

colnames(D19_6791) <- paste0(colnames(D19_6791), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6775), colnames(D19_6791))

D19_6775 <- D19_6775[, joint.bcs]
D19_6791 <- as.matrix(D19_6791[, joint.bcs])

# Setup Seurat object
D19_6775 <- CreateSeuratObject(counts = D19_6775)
gc()

# Normalize RNA data with log normalization
D19_6775 <- NormalizeData(D19_6775)
# Find and scale variable features
D19_6775 <- FindVariableFeatures(D19_6775, selection.method = "mean.var.plot")
D19_6775 <- ScaleData(D19_6775, features = VariableFeatures(D19_6775))

gc()

# Add HTO data as a new assay independent from RNA
D19_6775[["HTO"]] <- CreateAssayObject(counts = D19_6791)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6775 <- NormalizeData(D19_6775, assay = "HTO", normalization.method = "CLR")

D19_6775 <- HTODemux(D19_6775, assay = "HTO", positive.quantile = 0.99)

table(D19_6775$HTO_classification.global)

Idents(D19_6775) <- "HTO_classification.global"

D19_6775 <- subset(D19_6775, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6791, joint.bcs); gc()

Idents(D19_6775) <- "hash.ID"

D19_6775 <- subset(D19_6775, 
                   idents = c("SZ7-TTCCTGCC", "SZ8-CCGTACCT", "SZ9-GCCTAGTA"),
                   invert = TRUE); gc()


D19_6776 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6776/outs/filtered_feature_bc_matrix/")

D19_6792 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6792/umi_count/",
                    gene.column = 1) 

colnames(D19_6792) <- paste0(colnames(D19_6792), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6776), colnames(D19_6792))

D19_6776 <- D19_6776[, joint.bcs]
D19_6792 <- as.matrix(D19_6792[, joint.bcs])

# Setup Seurat object
D19_6776 <- CreateSeuratObject(counts = D19_6776)
gc()

# Normalize RNA data with log normalization
D19_6776 <- NormalizeData(D19_6776)
# Find and scale variable features
D19_6776 <- FindVariableFeatures(D19_6776, selection.method = "mean.var.plot")
D19_6776 <- ScaleData(D19_6776, features = VariableFeatures(D19_6776))

gc()

# Add HTO data as a new assay independent from RNA
D19_6776[["HTO"]] <- CreateAssayObject(counts = D19_6792)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6776 <- NormalizeData(D19_6776, assay = "HTO", normalization.method = "CLR")

D19_6776 <- HTODemux(D19_6776, assay = "HTO", positive.quantile = 0.99)

table(D19_6776$HTO_classification.global)

Idents(D19_6776) <- "HTO_classification.global"

D19_6776 <- subset(D19_6776, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6792, joint.bcs); gc()

Idents(D19_6776) <- "hash.ID"

D19_6776 <- subset(D19_6776, 
                   idents = c("SZ7-TTCCTGCC", "SZ8-CCGTACCT", "SZ9-GCCTAGTA"),
                   invert = TRUE); gc()


D19_6777 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6777/outs/filtered_feature_bc_matrix/")

D19_6793 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6793/umi_count/",
                    gene.column = 1) 

colnames(D19_6793) <- paste0(colnames(D19_6793), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6777), colnames(D19_6793))

D19_6777 <- D19_6777[, joint.bcs]
D19_6793 <- as.matrix(D19_6793[, joint.bcs])

# Setup Seurat object
D19_6777 <- CreateSeuratObject(counts = D19_6777)
gc()

# Normalize RNA data with log normalization
D19_6777 <- NormalizeData(D19_6777)
# Find and scale variable features
D19_6777 <- FindVariableFeatures(D19_6777, selection.method = "mean.var.plot")
D19_6777 <- ScaleData(D19_6777, features = VariableFeatures(D19_6777))

gc()

# Add HTO data as a new assay independent from RNA
D19_6777[["HTO"]] <- CreateAssayObject(counts = D19_6793)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6777 <- NormalizeData(D19_6777, assay = "HTO", normalization.method = "CLR")

D19_6777 <- HTODemux(D19_6777, assay = "HTO", positive.quantile = 0.99)

table(D19_6777$HTO_classification.global)

Idents(D19_6777) <- "HTO_classification.global"

D19_6777 <- subset(D19_6777, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6793, joint.bcs); gc()

Idents(D19_6777) <- "hash.ID"

D19_6777 <- subset(D19_6777, 
                   idents = c("SZ7-TTCCTGCC", "SZ8-CCGTACCT", "SZ9-GCCTAGTA"),
                   invert = TRUE); gc()



D19_6778 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6778/outs/filtered_feature_bc_matrix/")

D19_6794 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6794/umi_count/",
                    gene.column = 1) 

colnames(D19_6794) <- paste0(colnames(D19_6794), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6778), colnames(D19_6794))

D19_6778 <- D19_6778[, joint.bcs]
D19_6794 <- as.matrix(D19_6794[, joint.bcs])

# Setup Seurat object
D19_6778 <- CreateSeuratObject(counts = D19_6778)
gc()

# Normalize RNA data with log normalization
D19_6778 <- NormalizeData(D19_6778)
# Find and scale variable features
D19_6778 <- FindVariableFeatures(D19_6778, selection.method = "mean.var.plot")
D19_6778 <- ScaleData(D19_6778, features = VariableFeatures(D19_6778))

gc()

# Add HTO data as a new assay independent from RNA
D19_6778[["HTO"]] <- CreateAssayObject(counts = D19_6794)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6778 <- NormalizeData(D19_6778, assay = "HTO", normalization.method = "CLR")

D19_6778 <- HTODemux(D19_6778, assay = "HTO", positive.quantile = 0.99)

table(D19_6778$HTO_classification.global)

Idents(D19_6778) <- "HTO_classification.global"

D19_6778 <- subset(D19_6778, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6794, joint.bcs); gc()

Idents(D19_6778) <- "hash.ID"

D19_6778 <- subset(D19_6778, 
                   idents = c("SZ7-TTCCTGCC", "SZ8-CCGTACCT", "SZ9-GCCTAGTA"),
                   invert = TRUE); gc()


# merge all batch 3 seurat objects


batch3 <- merge(D19_6771, y = c(D19_6772, D19_6773, D19_6774,
                D19_6775, D19_6776, D19_6777, D19_6778),
                project = "Batch3")

rm(D19_6771, D19_6772, D19_6773, D19_6774,
   D19_6775, D19_6776, D19_6777, D19_6778); gc()

batch3[["percent.mt"]] <- PercentageFeatureSet(batch3, pattern = "^MT-"); gc()
# Visualize QC metrics as a violin plot
VlnPlot(batch3,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,group.by = "orig.ident")

batch3 <- subset(batch3, 
               subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 15); gc()


saveRDS(batch3, "/media/leandro/SAMSUNG/FASTQs/batch3.rds")



################################### Batch 4 #####################################

D19_6779 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6779/outs/filtered_feature_bc_matrix/")

D19_6795 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6795/umi_count/",
                    gene.column = 1) 

colnames(D19_6795) <- paste0(colnames(D19_6795), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6779), colnames(D19_6795))

D19_6779 <- D19_6779[, joint.bcs]
D19_6795 <- as.matrix(D19_6795[, joint.bcs])

# Setup Seurat object
D19_6779 <- CreateSeuratObject(counts = D19_6779)
gc()

# Normalize RNA data with log normalization
D19_6779 <- NormalizeData(D19_6779)
# Find and scale variable features
D19_6779 <- FindVariableFeatures(D19_6779, selection.method = "mean.var.plot")
D19_6779 <- ScaleData(D19_6779, features = VariableFeatures(D19_6779))

gc()

# Add HTO data as a new assay independent from RNA
D19_6779[["HTO"]] <- CreateAssayObject(counts = D19_6795)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6779 <- NormalizeData(D19_6779, assay = "HTO", normalization.method = "CLR")

D19_6779 <- HTODemux(D19_6779, assay = "HTO", positive.quantile = 0.99)

table(D19_6779$HTO_classification.global)

Idents(D19_6779) <- "HTO_classification.global"

D19_6779 <- subset(D19_6779, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6795, joint.bcs); gc()

Idents(D19_6779) <- "hash.ID"

D19_6779 <- subset(D19_6779, 
                   idents = c("SZ12-TGACGCCG", "SZ11-CTTATCAC", "SZ10-CCGTACCT"),
                   invert = TRUE); gc()


D19_6780 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6780/outs/filtered_feature_bc_matrix/")

D19_6796 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6796/umi_count/",
                    gene.column = 1) 

colnames(D19_6796) <- paste0(colnames(D19_6796), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6780), colnames(D19_6796))

D19_6780 <- D19_6780[, joint.bcs]
D19_6796 <- as.matrix(D19_6796[, joint.bcs])

# Setup Seurat object
D19_6780 <- CreateSeuratObject(counts = D19_6780)
gc()

# Normalize RNA data with log normalization
D19_6780 <- NormalizeData(D19_6780)
# Find and scale variable features
D19_6780 <- FindVariableFeatures(D19_6780, selection.method = "mean.var.plot")
D19_6780 <- ScaleData(D19_6780, features = VariableFeatures(D19_6780))

gc()

# Add HTO data as a new assay independent from RNA
D19_6780[["HTO"]] <- CreateAssayObject(counts = D19_6796)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6780 <- NormalizeData(D19_6780, assay = "HTO", normalization.method = "CLR")

D19_6780 <- HTODemux(D19_6780, assay = "HTO", positive.quantile = 0.99)

table(D19_6780$HTO_classification.global)

Idents(D19_6780) <- "HTO_classification.global"

D19_6780 <- subset(D19_6780, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6796, joint.bcs); gc()

Idents(D19_6780) <- "hash.ID"

D19_6780 <- subset(D19_6780, 
                   idents = c("SZ12-TGACGCCG", "SZ11-CTTATCAC", "SZ10-CCGTACCT"),
                   invert = TRUE); gc()


D19_6781 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6781/outs/filtered_feature_bc_matrix/")

D19_6797 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6797/umi_count/",
                    gene.column = 1) 

colnames(D19_6797) <- paste0(colnames(D19_6797), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6781), colnames(D19_6797))

D19_6781 <- D19_6781[, joint.bcs]
D19_6797 <- as.matrix(D19_6797[, joint.bcs])

# Setup Seurat object
D19_6781 <- CreateSeuratObject(counts = D19_6781)
gc()

# Normalize RNA data with log normalization
D19_6781 <- NormalizeData(D19_6781)
# Find and scale variable features
D19_6781 <- FindVariableFeatures(D19_6781, selection.method = "mean.var.plot")
D19_6781 <- ScaleData(D19_6781, features = VariableFeatures(D19_6781))

gc()

# Add HTO data as a new assay independent from RNA
D19_6781[["HTO"]] <- CreateAssayObject(counts = D19_6797)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6781 <- NormalizeData(D19_6781, assay = "HTO", normalization.method = "CLR")

D19_6781 <- HTODemux(D19_6781, assay = "HTO", positive.quantile = 0.99)

table(D19_6781$HTO_classification.global)

Idents(D19_6781) <- "HTO_classification.global"

D19_6781 <- subset(D19_6781, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6797, joint.bcs); gc()

Idents(D19_6781) <- "hash.ID"

D19_6781 <- subset(D19_6781, 
                   idents = c("SZ12-TGACGCCG", "SZ11-CTTATCAC", "SZ10-CCGTACCT"),
                   invert = TRUE); gc()


D19_6782 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6782/outs/filtered_feature_bc_matrix/")

D19_6798 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6798/umi_count/",
                    gene.column = 1) 

colnames(D19_6798) <- paste0(colnames(D19_6798), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6782), colnames(D19_6798))

D19_6782 <- D19_6782[, joint.bcs]
D19_6798 <- as.matrix(D19_6798[, joint.bcs])

# Setup Seurat object
D19_6782 <- CreateSeuratObject(counts = D19_6782)
gc()

# Normalize RNA data with log normalization
D19_6782 <- NormalizeData(D19_6782)
# Find and scale variable features
D19_6782 <- FindVariableFeatures(D19_6782, selection.method = "mean.var.plot")
D19_6782 <- ScaleData(D19_6782, features = VariableFeatures(D19_6782))

gc()

# Add HTO data as a new assay independent from RNA
D19_6782[["HTO"]] <- CreateAssayObject(counts = D19_6798)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6782 <- NormalizeData(D19_6782, assay = "HTO", normalization.method = "CLR")

D19_6782 <- HTODemux(D19_6782, assay = "HTO", positive.quantile = 0.99)

table(D19_6782$HTO_classification.global)

Idents(D19_6782) <- "HTO_classification.global"

D19_6782 <- subset(D19_6782, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6798, joint.bcs); gc()

Idents(D19_6782) <- "hash.ID"

D19_6782 <- subset(D19_6782, 
                   idents = c("SZ12-TGACGCCG", "SZ11-CTTATCAC", "SZ10-CCGTACCT"),
                   invert = TRUE); gc()



D19_6783 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6783/outs/filtered_feature_bc_matrix/")

D19_6799 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6799/umi_count/",
                    gene.column = 1) 

colnames(D19_6799) <- paste0(colnames(D19_6799), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6783), colnames(D19_6799))

D19_6783 <- D19_6783[, joint.bcs]
D19_6799 <- as.matrix(D19_6799[, joint.bcs])

# Setup Seurat object
D19_6783 <- CreateSeuratObject(counts = D19_6783)
gc()

# Normalize RNA data with log normalization
D19_6783 <- NormalizeData(D19_6783)
# Find and scale variable features
D19_6783 <- FindVariableFeatures(D19_6783, selection.method = "mean.var.plot")
D19_6783 <- ScaleData(D19_6783, features = VariableFeatures(D19_6783))

gc()

# Add HTO data as a new assay independent from RNA
D19_6783[["HTO"]] <- CreateAssayObject(counts = D19_6799)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6783 <- NormalizeData(D19_6783, assay = "HTO", normalization.method = "CLR")

D19_6783 <- HTODemux(D19_6783, assay = "HTO", positive.quantile = 0.99)

table(D19_6783$HTO_classification.global)

Idents(D19_6783) <- "HTO_classification.global"

D19_6783 <- subset(D19_6783, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6799, joint.bcs); gc()

Idents(D19_6783) <- "hash.ID"

D19_6783 <- subset(D19_6783, 
                   idents = c("SZ12-TGACGCCG", "SZ11-CTTATCAC", "SZ10-CCGTACCT"),
                   invert = TRUE); gc()



D19_6784 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6784/outs/filtered_feature_bc_matrix/")

D19_6800 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6800/umi_count/",
                    gene.column = 1) 

colnames(D19_6800) <- paste0(colnames(D19_6800), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6784), colnames(D19_6800))

D19_6784 <- D19_6784[, joint.bcs]
D19_6800 <- as.matrix(D19_6800[, joint.bcs])

# Setup Seurat object
D19_6784 <- CreateSeuratObject(counts = D19_6784)
gc()

# Normalize RNA data with log normalization
D19_6784 <- NormalizeData(D19_6784)
# Find and scale variable features
D19_6784 <- FindVariableFeatures(D19_6784, selection.method = "mean.var.plot")
D19_6784 <- ScaleData(D19_6784, features = VariableFeatures(D19_6784))

gc()

# Add HTO data as a new assay independent from RNA
D19_6784[["HTO"]] <- CreateAssayObject(counts = D19_6800)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6784 <- NormalizeData(D19_6784, assay = "HTO", normalization.method = "CLR")

D19_6784 <- HTODemux(D19_6784, assay = "HTO", positive.quantile = 0.99)

table(D19_6784$HTO_classification.global)

Idents(D19_6784) <- "HTO_classification.global"

D19_6784 <- subset(D19_6784, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6800, joint.bcs); gc()

Idents(D19_6784) <- "hash.ID"

D19_6784 <- subset(D19_6784, 
                   idents = c("SZ12-TGACGCCG", "SZ11-CTTATCAC", "SZ10-CCGTACCT"),
                   invert = TRUE); gc()


D19_6785 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6785/outs/filtered_feature_bc_matrix/")

D19_6801 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6801/umi_count/",
                    gene.column = 1) 

colnames(D19_6801) <- paste0(colnames(D19_6801), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6785), colnames(D19_6801))

D19_6785 <- D19_6785[, joint.bcs]
D19_6801 <- as.matrix(D19_6801[, joint.bcs])

# Setup Seurat object
D19_6785 <- CreateSeuratObject(counts = D19_6785)
gc()

# Normalize RNA data with log normalization
D19_6785 <- NormalizeData(D19_6785)
# Find and scale variable features
D19_6785 <- FindVariableFeatures(D19_6785, selection.method = "mean.var.plot")
D19_6785 <- ScaleData(D19_6785, features = VariableFeatures(D19_6785))

gc()

# Add HTO data as a new assay independent from RNA
D19_6785[["HTO"]] <- CreateAssayObject(counts = D19_6801)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6785 <- NormalizeData(D19_6785, assay = "HTO", normalization.method = "CLR")

D19_6785 <- HTODemux(D19_6785, assay = "HTO", positive.quantile = 0.99)

table(D19_6785$HTO_classification.global)

Idents(D19_6785) <- "HTO_classification.global"

D19_6785 <- subset(D19_6785, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6801, joint.bcs); gc()

Idents(D19_6785) <- "hash.ID"

D19_6785 <- subset(D19_6785, 
                   idents = c("SZ12-TGACGCCG", "SZ11-CTTATCAC", "SZ10-CCGTACCT"),
                   invert = TRUE); gc()



D19_6786 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-6786/outs/filtered_feature_bc_matrix/")

D19_6802 <- Read10X("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/CITE-Seq/190620KelA_D19-6802/umi_count/",
                    gene.column = 1) 

colnames(D19_6802) <- paste0(colnames(D19_6802), "-1")

gc()

joint.bcs <- intersect(colnames(D19_6786), colnames(D19_6802))

D19_6786 <- D19_6786[, joint.bcs]
D19_6802 <- as.matrix(D19_6802[, joint.bcs])

# Setup Seurat object
D19_6786 <- CreateSeuratObject(counts = D19_6786)
gc()

# Normalize RNA data with log normalization
D19_6786 <- NormalizeData(D19_6786)
# Find and scale variable features
D19_6786 <- FindVariableFeatures(D19_6786, selection.method = "mean.var.plot")
D19_6786 <- ScaleData(D19_6786, features = VariableFeatures(D19_6786))

gc()

# Add HTO data as a new assay independent from RNA
D19_6786[["HTO"]] <- CreateAssayObject(counts = D19_6802)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_6786 <- NormalizeData(D19_6786, assay = "HTO", normalization.method = "CLR")

D19_6786 <- HTODemux(D19_6786, assay = "HTO", positive.quantile = 0.99)

table(D19_6786$HTO_classification.global)

Idents(D19_6786) <- "HTO_classification.global"

D19_6786 <- subset(D19_6786, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_6802, joint.bcs); gc()

Idents(D19_6786) <- "hash.ID"

D19_6786 <- subset(D19_6786, 
                   idents = c("SZ12-TGACGCCG", "SZ11-CTTATCAC", "SZ10-CCGTACCT"),
                   invert = TRUE); gc()

# merge all batch 4 seurat objects


batch4 <- merge(D19_6779, y = c(D19_6780, D19_6781, D19_6782,
                                D19_6783, D19_6784, D19_6785, D19_6786),
                project = "Batch4")

rm(D19_6779, D19_6780, D19_6781, D19_6782,
   D19_6783, D19_6784, D19_6785, D19_6786); gc()

batch4[["percent.mt"]] <- PercentageFeatureSet(batch4, pattern = "^MT-"); gc()
# Visualize QC metrics as a violin plot
VlnPlot(batch4,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,group.by = "orig.ident")

batch4 <- subset(batch4, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 15); gc()


saveRDS(batch4, "/media/leandro/SAMSUNG/FASTQs/batch4.rds")

################################### Batch 5 #####################################

D19_7380 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7380/outs/filtered_feature_bc_matrix/")

D19_7348 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7348_S1/umi_count/",
                    gene.column = 1) 

colnames(D19_7348) <- paste0(colnames(D19_7348), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7380), colnames(D19_7348))

D19_7380 <- D19_7380[, joint.bcs]
D19_7348 <- as.matrix(D19_7348[, joint.bcs])

# Setup Seurat object
D19_7380 <- CreateSeuratObject(counts = D19_7380)
gc()

# Normalize RNA data with log normalization
D19_7380 <- NormalizeData(D19_7380)
# Find and scale variable features
D19_7380 <- FindVariableFeatures(D19_7380, selection.method = "mean.var.plot")
D19_7380 <- ScaleData(D19_7380, features = VariableFeatures(D19_7380))

gc()

# Add HTO data as a new assay independent from RNA
D19_7380[["HTO"]] <- CreateAssayObject(counts = D19_7348)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7380 <- NormalizeData(D19_7380, assay = "HTO", normalization.method = "CLR")

D19_7380 <- HTODemux(D19_7380, assay = "HTO", positive.quantile = 0.99)

table(D19_7380$HTO_classification.global)

Idents(D19_7380) <- "HTO_classification.global"

D19_7380 <- subset(D19_7380, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7348, joint.bcs); gc()

Idents(D19_7380) <- "hash.ID"

D19_7380 <- subset(D19_7380, 
                   idents = c("SZ13-GGTAGATG", "SZ15-GCCTAGTA", "SZ14-CTCGAACG"),
                   invert = TRUE); gc()


D19_7381 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7381/outs/filtered_feature_bc_matrix/")

D19_7349 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7349_S2/umi_count/",
                    gene.column = 1) 

colnames(D19_7349) <- paste0(colnames(D19_7349), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7381), colnames(D19_7349))

D19_7381 <- D19_7381[, joint.bcs]
D19_7349 <- as.matrix(D19_7349[, joint.bcs])

# Setup Seurat object
D19_7381 <- CreateSeuratObject(counts = D19_7381)
gc()

# Normalize RNA data with log normalization
D19_7381 <- NormalizeData(D19_7381)
# Find and scale variable features
D19_7381 <- FindVariableFeatures(D19_7381, selection.method = "mean.var.plot")
D19_7381 <- ScaleData(D19_7381, features = VariableFeatures(D19_7381))

gc()

# Add HTO data as a new assay independent from RNA
D19_7381[["HTO"]] <- CreateAssayObject(counts = D19_7349)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7381 <- NormalizeData(D19_7381, assay = "HTO", normalization.method = "CLR")

D19_7381 <- HTODemux(D19_7381, assay = "HTO", positive.quantile = 0.99)

table(D19_7381$HTO_classification.global)

Idents(D19_7381) <- "HTO_classification.global"

D19_7381 <- subset(D19_7381, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7349, joint.bcs); gc()

Idents(D19_7381) <- "hash.ID"

D19_7381 <- subset(D19_7381, 
                   idents = c("SZ13-GGTAGATG", "SZ15-GCCTAGTA", "SZ14-CTCGAACG"),
                   invert = TRUE); gc()


D19_7382 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7382/outs/filtered_feature_bc_matrix/")

D19_7350 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7350_S3/umi_count/",
                    gene.column = 1) 

colnames(D19_7350) <- paste0(colnames(D19_7350), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7382), colnames(D19_7350))

D19_7382 <- D19_7382[, joint.bcs]
D19_7350 <- as.matrix(D19_7350[, joint.bcs])

# Setup Seurat object
D19_7382 <- CreateSeuratObject(counts = D19_7382)
gc()

# Normalize RNA data with log normalization
D19_7382 <- NormalizeData(D19_7382)
# Find and scale variable features
D19_7382 <- FindVariableFeatures(D19_7382, selection.method = "mean.var.plot")
D19_7382 <- ScaleData(D19_7382, features = VariableFeatures(D19_7382))

gc()

# Add HTO data as a new assay independent from RNA
D19_7382[["HTO"]] <- CreateAssayObject(counts = D19_7350)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7382 <- NormalizeData(D19_7382, assay = "HTO", normalization.method = "CLR")

D19_7382 <- HTODemux(D19_7382, assay = "HTO", positive.quantile = 0.99)

table(D19_7382$HTO_classification.global)

Idents(D19_7382) <- "HTO_classification.global"

D19_7382 <- subset(D19_7382, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7350, joint.bcs); gc()

Idents(D19_7382) <- "hash.ID"

D19_7382 <- subset(D19_7382, 
                   idents = c("SZ13-GGTAGATG", "SZ15-GCCTAGTA", "SZ14-CTCGAACG"),
                   invert = TRUE); gc()


D19_7383 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7383/outs/filtered_feature_bc_matrix/")

D19_7351 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7351_S4/umi_count/",
                    gene.column = 1) 

colnames(D19_7351) <- paste0(colnames(D19_7351), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7383), colnames(D19_7351))

D19_7383 <- D19_7383[, joint.bcs]
D19_7351 <- as.matrix(D19_7351[, joint.bcs])

# Setup Seurat object
D19_7383 <- CreateSeuratObject(counts = D19_7383)
gc()

# Normalize RNA data with log normalization
D19_7383 <- NormalizeData(D19_7383)
# Find and scale variable features
D19_7383 <- FindVariableFeatures(D19_7383, selection.method = "mean.var.plot")
D19_7383 <- ScaleData(D19_7383, features = VariableFeatures(D19_7383))

gc()

# Add HTO data as a new assay independent from RNA
D19_7383[["HTO"]] <- CreateAssayObject(counts = D19_7351)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7383 <- NormalizeData(D19_7383, assay = "HTO", normalization.method = "CLR")

D19_7383 <- HTODemux(D19_7383, assay = "HTO", positive.quantile = 0.99)

table(D19_7383$HTO_classification.global)

Idents(D19_7383) <- "HTO_classification.global"

D19_7383 <- subset(D19_7383, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7351, joint.bcs); gc()

Idents(D19_7383) <- "hash.ID"

D19_7383 <- subset(D19_7383, 
                   idents = c("SZ13-GGTAGATG", "SZ15-GCCTAGTA", "SZ14-CTCGAACG"),
                   invert = TRUE); gc()



D19_7384 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7384/outs/filtered_feature_bc_matrix/")

D19_7352 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7352_S5/umi_count/",
                    gene.column = 1) 

colnames(D19_7352) <- paste0(colnames(D19_7352), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7384), colnames(D19_7352))

D19_7384 <- D19_7384[, joint.bcs]
D19_7352 <- as.matrix(D19_7352[, joint.bcs])

# Setup Seurat object
D19_7384 <- CreateSeuratObject(counts = D19_7384)
gc()

# Normalize RNA data with log normalization
D19_7384 <- NormalizeData(D19_7384)
# Find and scale variable features
D19_7384 <- FindVariableFeatures(D19_7384, selection.method = "mean.var.plot")
D19_7384 <- ScaleData(D19_7384, features = VariableFeatures(D19_7384))

gc()

# Add HTO data as a new assay independent from RNA
D19_7384[["HTO"]] <- CreateAssayObject(counts = D19_7352)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7384 <- NormalizeData(D19_7384, assay = "HTO", normalization.method = "CLR")

D19_7384 <- HTODemux(D19_7384, assay = "HTO", positive.quantile = 0.99)

table(D19_7384$HTO_classification.global)

Idents(D19_7384) <- "HTO_classification.global"

D19_7384 <- subset(D19_7384, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7352, joint.bcs); gc()

Idents(D19_7384) <- "hash.ID"

D19_7384 <- subset(D19_7384, 
                   idents = c("SZ13-GGTAGATG", "SZ15-GCCTAGTA", "SZ14-CTCGAACG"),
                   invert = TRUE); gc()


D19_7385 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7385/outs/filtered_feature_bc_matrix/")

D19_7353 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7353_S6/umi_count/",
                    gene.column = 1) 

colnames(D19_7353) <- paste0(colnames(D19_7353), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7385), colnames(D19_7353))

D19_7385 <- D19_7385[, joint.bcs]
D19_7353 <- as.matrix(D19_7353[, joint.bcs])

# Setup Seurat object
D19_7385 <- CreateSeuratObject(counts = D19_7385)
gc()

# Normalize RNA data with log normalization
D19_7385 <- NormalizeData(D19_7385)
# Find and scale variable features
D19_7385 <- FindVariableFeatures(D19_7385, selection.method = "mean.var.plot")
D19_7385 <- ScaleData(D19_7385, features = VariableFeatures(D19_7385))

gc()

# Add HTO data as a new assay independent from RNA
D19_7385[["HTO"]] <- CreateAssayObject(counts = D19_7353)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7385 <- NormalizeData(D19_7385, assay = "HTO", normalization.method = "CLR")

D19_7385 <- HTODemux(D19_7385, assay = "HTO", positive.quantile = 0.99)

table(D19_7385$HTO_classification.global)

Idents(D19_7385) <- "HTO_classification.global"

D19_7385 <- subset(D19_7385, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7353, joint.bcs); gc()

Idents(D19_7385) <- "hash.ID"

D19_7385 <- subset(D19_7385, 
                   idents = c("SZ13-GGTAGATG", "SZ15-GCCTAGTA", "SZ14-CTCGAACG"),
                   invert = TRUE); gc()


D19_7386 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7386/outs/filtered_feature_bc_matrix/")

D19_7354 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7354_S7/umi_count/",
                    gene.column = 1) 

colnames(D19_7354) <- paste0(colnames(D19_7354), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7386), colnames(D19_7354))

D19_7386 <- D19_7386[, joint.bcs]
D19_7354 <- as.matrix(D19_7354[, joint.bcs])

# Setup Seurat object
D19_7386 <- CreateSeuratObject(counts = D19_7386)
gc()

# Normalize RNA data with log normalization
D19_7386 <- NormalizeData(D19_7386)
# Find and scale variable features
D19_7386 <- FindVariableFeatures(D19_7386, selection.method = "mean.var.plot")
D19_7386 <- ScaleData(D19_7386, features = VariableFeatures(D19_7386))

gc()

# Add HTO data as a new assay independent from RNA
D19_7386[["HTO"]] <- CreateAssayObject(counts = D19_7354)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7386 <- NormalizeData(D19_7386, assay = "HTO", normalization.method = "CLR")

D19_7386 <- HTODemux(D19_7386, assay = "HTO", positive.quantile = 0.99)

table(D19_7386$HTO_classification.global)

Idents(D19_7386) <- "HTO_classification.global"

D19_7386 <- subset(D19_7386, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7354, joint.bcs); gc()

Idents(D19_7386) <- "hash.ID"

D19_7386 <- subset(D19_7386, 
                   idents = c("SZ13-GGTAGATG", "SZ15-GCCTAGTA", "SZ14-CTCGAACG"),
                   invert = TRUE); gc()

table(D19_7386$HTO_classification)



D19_7387 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7387/outs/filtered_feature_bc_matrix/")

D19_7355 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7355_S8/umi_count/",
                    gene.column = 1) 

colnames(D19_7355) <- paste0(colnames(D19_7355), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7387), colnames(D19_7355))

D19_7387 <- D19_7387[, joint.bcs]
D19_7355 <- as.matrix(D19_7355[, joint.bcs])

# Setup Seurat object
D19_7387 <- CreateSeuratObject(counts = D19_7387)
gc()

# Normalize RNA data with log normalization
D19_7387 <- NormalizeData(D19_7387)
# Find and scale variable features
D19_7387 <- FindVariableFeatures(D19_7387, selection.method = "mean.var.plot")
D19_7387 <- ScaleData(D19_7387, features = VariableFeatures(D19_7387))

gc()

# Add HTO data as a new assay independent from RNA
D19_7387[["HTO"]] <- CreateAssayObject(counts = D19_7355)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7387 <- NormalizeData(D19_7387, assay = "HTO", normalization.method = "CLR")

D19_7387 <- HTODemux(D19_7387, assay = "HTO", positive.quantile = 0.99)

table(D19_7387$HTO_classification.global)

Idents(D19_7387) <- "HTO_classification.global"

D19_7387 <- subset(D19_7387, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7355, joint.bcs); gc()

Idents(D19_7387) <- "hash.ID"

D19_7387 <- subset(D19_7387, 
                   idents = c("SZ13-GGTAGATG", "SZ15-GCCTAGTA", "SZ14-CTCGAACG"),
                   invert = TRUE); gc()

table(D19_7387$HTO_classification)


# merge all batch 5 seurat objects


batch5 <- merge(D19_7380, y = c(D19_7381, D19_7382, D19_7383,
                                D19_7384, D19_7385, D19_7386, D19_7387),
                project = "Batch5")

rm(D19_7380, D19_7381, D19_7382, D19_7383,
   D19_7384, D19_7385, D19_7386, D19_7387); gc()

batch5[["percent.mt"]] <- PercentageFeatureSet(batch5, pattern = "^MT-"); gc()
# Visualize QC metrics as a violin plot
VlnPlot(batch5,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,group.by = "orig.ident")

batch5 <- subset(batch5, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 15); gc()


saveRDS(batch5, "/media/leandro/SAMSUNG/FASTQs/batch5.rds")


################################### Batch 6 #####################################


D19_7388 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7388/outs/filtered_feature_bc_matrix/")

D19_7356 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7356_S9/umi_count/",
                    gene.column = 1) 

colnames(D19_7356) <- paste0(colnames(D19_7356), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7388), colnames(D19_7356))

D19_7388 <- D19_7388[, joint.bcs]
D19_7356 <- as.matrix(D19_7356[, joint.bcs])

# Setup Seurat object
D19_7388 <- CreateSeuratObject(counts = D19_7388)
gc()

# Normalize RNA data with log normalization
D19_7388 <- NormalizeData(D19_7388)
# Find and scale variable features
D19_7388 <- FindVariableFeatures(D19_7388, selection.method = "mean.var.plot")
D19_7388 <- ScaleData(D19_7388, features = VariableFeatures(D19_7388))

gc()

# Add HTO data as a new assay independent from RNA
D19_7388[["HTO"]] <- CreateAssayObject(counts = D19_7356)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7388 <- NormalizeData(D19_7388, assay = "HTO", normalization.method = "CLR")

D19_7388 <- HTODemux(D19_7388, assay = "HTO", positive.quantile = 0.99)

table(D19_7388$HTO_classification.global)

Idents(D19_7388) <- "HTO_classification.global"

D19_7388 <- subset(D19_7388, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7356, joint.bcs); gc()

Idents(D19_7388) <- "hash.ID"

D19_7388 <- subset(D19_7388, 
                   idents = c("SZ16-ATGATGAA", "SZ17-CTCGAACG", "SZ18-CTTATCAC"),
                   invert = TRUE); gc()

table(D19_7388$HTO_classification)



D19_7389 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7389/outs/filtered_feature_bc_matrix/")

D19_7357 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7357_S10/umi_count/",
                    gene.column = 1) 

colnames(D19_7357) <- paste0(colnames(D19_7357), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7389), colnames(D19_7357))

D19_7389 <- D19_7389[, joint.bcs]
D19_7357 <- as.matrix(D19_7357[, joint.bcs])

# Setup Seurat object
D19_7389 <- CreateSeuratObject(counts = D19_7389)
gc()

# Normalize RNA data with log normalization
D19_7389 <- NormalizeData(D19_7389)
# Find and scale variable features
D19_7389 <- FindVariableFeatures(D19_7389, selection.method = "mean.var.plot")
D19_7389 <- ScaleData(D19_7389, features = VariableFeatures(D19_7389))

gc()

# Add HTO data as a new assay independent from RNA
D19_7389[["HTO"]] <- CreateAssayObject(counts = D19_7357)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7389 <- NormalizeData(D19_7389, assay = "HTO", normalization.method = "CLR")

D19_7389 <- HTODemux(D19_7389, assay = "HTO", positive.quantile = 0.99)

table(D19_7389$HTO_classification.global)

Idents(D19_7389) <- "HTO_classification.global"

D19_7389 <- subset(D19_7389, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7357, joint.bcs); gc()

Idents(D19_7389) <- "hash.ID"

D19_7389 <- subset(D19_7389, 
                   idents = c("SZ16-ATGATGAA", "SZ17-CTCGAACG", "SZ18-CTTATCAC",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7389$HTO_classification)



D19_7390 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7390/outs/filtered_feature_bc_matrix/")

D19_7358 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7358_S11/umi_count/",
                    gene.column = 1) 

colnames(D19_7358) <- paste0(colnames(D19_7358), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7390), colnames(D19_7358))

D19_7390 <- D19_7390[, joint.bcs]
D19_7358 <- as.matrix(D19_7358[, joint.bcs])

# Setup Seurat object
D19_7390 <- CreateSeuratObject(counts = D19_7390)
gc()

# Normalize RNA data with log normalization
D19_7390 <- NormalizeData(D19_7390)
# Find and scale variable features
D19_7390 <- FindVariableFeatures(D19_7390, selection.method = "mean.var.plot")
D19_7390 <- ScaleData(D19_7390, features = VariableFeatures(D19_7390))

gc()

# Add HTO data as a new assay independent from RNA
D19_7390[["HTO"]] <- CreateAssayObject(counts = D19_7358)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7390 <- NormalizeData(D19_7390, assay = "HTO", normalization.method = "CLR")

D19_7390 <- HTODemux(D19_7390, assay = "HTO", positive.quantile = 0.99)

table(D19_7390$HTO_classification.global)

Idents(D19_7390) <- "HTO_classification.global"

D19_7390 <- subset(D19_7390, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7358, joint.bcs); gc()

Idents(D19_7390) <- "hash.ID"

D19_7390 <- subset(D19_7390, 
                   idents = c("SZ16-ATGATGAA", "SZ17-CTCGAACG", "SZ18-CTTATCAC"),
                   invert = TRUE); gc()

table(D19_7390$HTO_classification)


D19_7391 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7391/outs/filtered_feature_bc_matrix/")

D19_7359 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7359_S12/umi_count/",
                    gene.column = 1) 

colnames(D19_7359) <- paste0(colnames(D19_7359), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7391), colnames(D19_7359))

D19_7391 <- D19_7391[, joint.bcs]
D19_7359 <- as.matrix(D19_7359[, joint.bcs])

# Setup Seurat object
D19_7391 <- CreateSeuratObject(counts = D19_7391)
gc()

# Normalize RNA data with log normalization
D19_7391 <- NormalizeData(D19_7391)
# Find and scale variable features
D19_7391 <- FindVariableFeatures(D19_7391, selection.method = "mean.var.plot")
D19_7391 <- ScaleData(D19_7391, features = VariableFeatures(D19_7391))

gc()

# Add HTO data as a new assay independent from RNA
D19_7391[["HTO"]] <- CreateAssayObject(counts = D19_7359)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7391 <- NormalizeData(D19_7391, assay = "HTO", normalization.method = "CLR")

D19_7391 <- HTODemux(D19_7391, assay = "HTO", positive.quantile = 0.99)

table(D19_7391$HTO_classification.global)

Idents(D19_7391) <- "HTO_classification.global"

D19_7391 <- subset(D19_7391, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7359, joint.bcs); gc()

Idents(D19_7391) <- "hash.ID"

D19_7391 <- subset(D19_7391, 
                   idents = c("SZ16-ATGATGAA", "SZ17-CTCGAACG", "SZ18-CTTATCAC"),
                   invert = TRUE); gc()

table(D19_7391$HTO_classification)



D19_7392 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7392/outs/filtered_feature_bc_matrix/")

D19_7360 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7360_S13/umi_count/",
                    gene.column = 1) 

colnames(D19_7360) <- paste0(colnames(D19_7360), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7392), colnames(D19_7360))

D19_7392 <- D19_7392[, joint.bcs]
D19_7360 <- as.matrix(D19_7360[, joint.bcs])

# Setup Seurat object
D19_7392 <- CreateSeuratObject(counts = D19_7392)
gc()

# Normalize RNA data with log normalization
D19_7392 <- NormalizeData(D19_7392)
# Find and scale variable features
D19_7392 <- FindVariableFeatures(D19_7392, selection.method = "mean.var.plot")
D19_7392 <- ScaleData(D19_7392, features = VariableFeatures(D19_7392))

gc()

# Add HTO data as a new assay independent from RNA
D19_7392[["HTO"]] <- CreateAssayObject(counts = D19_7360)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7392 <- NormalizeData(D19_7392, assay = "HTO", normalization.method = "CLR")

D19_7392 <- HTODemux(D19_7392, assay = "HTO", positive.quantile = 0.99)

table(D19_7392$HTO_classification.global)

Idents(D19_7392) <- "HTO_classification.global"

D19_7392 <- subset(D19_7392, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7360, joint.bcs); gc()

Idents(D19_7392) <- "hash.ID"

D19_7392 <- subset(D19_7392, 
                   idents = c("SZ16-ATGATGAA", "SZ17-CTCGAACG", "SZ18-CTTATCAC"),
                   invert = TRUE); gc()

table(D19_7392$HTO_classification)



D19_7393 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7393/outs/filtered_feature_bc_matrix/")

D19_7361 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7361_S14/umi_count/",
                    gene.column = 1) 

colnames(D19_7361) <- paste0(colnames(D19_7361), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7393), colnames(D19_7361))

D19_7393 <- D19_7393[, joint.bcs]
D19_7361 <- as.matrix(D19_7361[, joint.bcs])

# Setup Seurat object
D19_7393 <- CreateSeuratObject(counts = D19_7393)
gc()

# Normalize RNA data with log normalization
D19_7393 <- NormalizeData(D19_7393)
# Find and scale variable features
D19_7393 <- FindVariableFeatures(D19_7393, selection.method = "mean.var.plot")
D19_7393 <- ScaleData(D19_7393, features = VariableFeatures(D19_7393))

gc()

# Add HTO data as a new assay independent from RNA
D19_7393[["HTO"]] <- CreateAssayObject(counts = D19_7361)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7393 <- NormalizeData(D19_7393, assay = "HTO", normalization.method = "CLR")

D19_7393 <- HTODemux(D19_7393, assay = "HTO", positive.quantile = 0.99)

table(D19_7393$HTO_classification.global)

Idents(D19_7393) <- "HTO_classification.global"

D19_7393 <- subset(D19_7393, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7361, joint.bcs); gc()

Idents(D19_7393) <- "hash.ID"

D19_7393 <- subset(D19_7393, 
                   idents = c("SZ16-ATGATGAA", "SZ17-CTCGAACG", "SZ18-CTTATCAC", 
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7393$HTO_classification)




D19_7394 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7394/outs/filtered_feature_bc_matrix/")

D19_7362 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7362_S15/umi_count/",
                    gene.column = 1) 

colnames(D19_7362) <- paste0(colnames(D19_7362), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7394), colnames(D19_7362))

D19_7394 <- D19_7394[, joint.bcs]
D19_7362 <- as.matrix(D19_7362[, joint.bcs])

# Setup Seurat object
D19_7394 <- CreateSeuratObject(counts = D19_7394)
gc()

# Normalize RNA data with log normalization
D19_7394 <- NormalizeData(D19_7394)
# Find and scale variable features
D19_7394 <- FindVariableFeatures(D19_7394, selection.method = "mean.var.plot")
D19_7394 <- ScaleData(D19_7394, features = VariableFeatures(D19_7394))

gc()

# Add HTO data as a new assay independent from RNA
D19_7394[["HTO"]] <- CreateAssayObject(counts = D19_7362)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7394 <- NormalizeData(D19_7394, assay = "HTO", normalization.method = "CLR")

D19_7394 <- HTODemux(D19_7394, assay = "HTO", positive.quantile = 0.99)

table(D19_7394$HTO_classification.global)

Idents(D19_7394) <- "HTO_classification.global"

D19_7394 <- subset(D19_7394, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7362, joint.bcs); gc()

Idents(D19_7394) <- "hash.ID"

D19_7394 <- subset(D19_7394, 
                   idents = c("SZ16-ATGATGAA", "SZ17-CTCGAACG", "SZ18-CTTATCAC", 
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7394$HTO_classification)



D19_7395 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7395/outs/filtered_feature_bc_matrix/")

D19_7363 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7363_S16/umi_count/",
                    gene.column = 1) 

colnames(D19_7363) <- paste0(colnames(D19_7363), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7395), colnames(D19_7363))

D19_7395 <- D19_7395[, joint.bcs]
D19_7363 <- as.matrix(D19_7363[, joint.bcs])

# Setup Seurat object
D19_7395 <- CreateSeuratObject(counts = D19_7395)
gc()

# Normalize RNA data with log normalization
D19_7395 <- NormalizeData(D19_7395)
# Find and scale variable features
D19_7395 <- FindVariableFeatures(D19_7395, selection.method = "mean.var.plot")
D19_7395 <- ScaleData(D19_7395, features = VariableFeatures(D19_7395))

gc()

# Add HTO data as a new assay independent from RNA
D19_7395[["HTO"]] <- CreateAssayObject(counts = D19_7363)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7395 <- NormalizeData(D19_7395, assay = "HTO", normalization.method = "CLR")

D19_7395 <- HTODemux(D19_7395, assay = "HTO", positive.quantile = 0.99)

table(D19_7395$HTO_classification.global)

Idents(D19_7395) <- "HTO_classification.global"

D19_7395 <- subset(D19_7395, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7363, joint.bcs); gc()

Idents(D19_7395) <- "hash.ID"

D19_7395 <- subset(D19_7395, 
                   idents = c("SZ16-ATGATGAA", "SZ17-CTCGAACG", "SZ18-CTTATCAC", 
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7395$HTO_classification)



# merge all batch 6 seurat objects


batch6 <- merge(D19_7388, y = c(D19_7389, D19_7390, D19_7391,
                                D19_7392, D19_7393, D19_7394, D19_7395),
                project = "Batch6")

rm(D19_7388, D19_7389, D19_7390, D19_7391,
   D19_7392, D19_7393, D19_7394, D19_7395); gc()

batch6[["percent.mt"]] <- PercentageFeatureSet(batch6, pattern = "^MT-"); gc()
# Visualize QC metrics as a violin plot
VlnPlot(batch6,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,group.by = "orig.ident")

batch6 <- subset(batch6, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 15); gc()


saveRDS(batch6, "/media/leandro/SAMSUNG/FASTQs/batch6.rds")


################################### Batch 7 #####################################

D19_7396 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7396/outs/filtered_feature_bc_matrix/")

D19_7364 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7364_S1/umi_count/",
                    gene.column = 1) 

colnames(D19_7364) <- paste0(colnames(D19_7364), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7396), colnames(D19_7364))

D19_7396 <- D19_7396[, joint.bcs]
D19_7364 <- as.matrix(D19_7364[, joint.bcs])

# Setup Seurat object
D19_7396 <- CreateSeuratObject(counts = D19_7396)
gc()

# Normalize RNA data with log normalization
D19_7396 <- NormalizeData(D19_7396)
# Find and scale variable features
D19_7396 <- FindVariableFeatures(D19_7396, selection.method = "mean.var.plot")
D19_7396 <- ScaleData(D19_7396, features = VariableFeatures(D19_7396))

gc()

# Add HTO data as a new assay independent from RNA
D19_7396[["HTO"]] <- CreateAssayObject(counts = D19_7364)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7396 <- NormalizeData(D19_7396, assay = "HTO", normalization.method = "CLR")

D19_7396 <- HTODemux(D19_7396, assay = "HTO", positive.quantile = 0.99)

table(D19_7396$HTO_classification.global)

Idents(D19_7396) <- "HTO_classification.global"

D19_7396 <- subset(D19_7396, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7364, joint.bcs); gc()

Idents(D19_7396) <- "hash.ID"

D19_7396 <- subset(D19_7396, 
                   idents = c("SZ20-TGACGCCG", "SZ19-TGGTGTCA", "SZ21-GCCTAGTA", 
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7396$HTO_classification)



D19_7397 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7397/outs/filtered_feature_bc_matrix/")

D19_7365 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7365_S2/umi_count/",
                    gene.column = 1) 

colnames(D19_7365) <- paste0(colnames(D19_7365), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7397), colnames(D19_7365))

D19_7397 <- D19_7397[, joint.bcs]
D19_7365 <- as.matrix(D19_7365[, joint.bcs])

# Setup Seurat object
D19_7397 <- CreateSeuratObject(counts = D19_7397)
gc()

# Normalize RNA data with log normalization
D19_7397 <- NormalizeData(D19_7397)
# Find and scale variable features
D19_7397 <- FindVariableFeatures(D19_7397, selection.method = "mean.var.plot")
D19_7397 <- ScaleData(D19_7397, features = VariableFeatures(D19_7397))

gc()

# Add HTO data as a new assay independent from RNA
D19_7397[["HTO"]] <- CreateAssayObject(counts = D19_7365)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7397 <- NormalizeData(D19_7397, assay = "HTO", normalization.method = "CLR")

D19_7397 <- HTODemux(D19_7397, assay = "HTO", positive.quantile = 0.99)

table(D19_7397$HTO_classification.global)

Idents(D19_7397) <- "HTO_classification.global"

D19_7397 <- subset(D19_7397, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7365, joint.bcs); gc()

Idents(D19_7397) <- "hash.ID"

D19_7397 <- subset(D19_7397, 
                   idents = c("SZ20-TGACGCCG", "SZ19-TGGTGTCA", "SZ21-GCCTAGTA"),
                   invert = TRUE); gc()

table(D19_7397$HTO_classification)



D19_7398 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7398/outs/filtered_feature_bc_matrix/")

D19_7366 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7366_S3/umi_count/",
                    gene.column = 1) 

colnames(D19_7366) <- paste0(colnames(D19_7366), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7398), colnames(D19_7366))

D19_7398 <- D19_7398[, joint.bcs]
D19_7366 <- as.matrix(D19_7366[, joint.bcs])

# Setup Seurat object
D19_7398 <- CreateSeuratObject(counts = D19_7398)
gc()

# Normalize RNA data with log normalization
D19_7398 <- NormalizeData(D19_7398)
# Find and scale variable features
D19_7398 <- FindVariableFeatures(D19_7398, selection.method = "mean.var.plot")
D19_7398 <- ScaleData(D19_7398, features = VariableFeatures(D19_7398))

gc()

# Add HTO data as a new assay independent from RNA
D19_7398[["HTO"]] <- CreateAssayObject(counts = D19_7366)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7398 <- NormalizeData(D19_7398, assay = "HTO", normalization.method = "CLR")

D19_7398 <- HTODemux(D19_7398, assay = "HTO", positive.quantile = 0.99)

table(D19_7398$HTO_classification.global)

Idents(D19_7398) <- "HTO_classification.global"

D19_7398 <- subset(D19_7398, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7366, joint.bcs); gc()

Idents(D19_7398) <- "hash.ID"

D19_7398 <- subset(D19_7398, 
                   idents = c("SZ20-TGACGCCG", "SZ19-TGGTGTCA", "SZ21-GCCTAGTA"),
                   invert = TRUE); gc()

table(D19_7398$HTO_classification)



D19_7399 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7399/outs/filtered_feature_bc_matrix/")

D19_7367 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7367_S4/umi_count/",
                    gene.column = 1) 

colnames(D19_7367) <- paste0(colnames(D19_7367), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7399), colnames(D19_7367))

D19_7399 <- D19_7399[, joint.bcs]
D19_7367 <- as.matrix(D19_7367[, joint.bcs])

# Setup Seurat object
D19_7399 <- CreateSeuratObject(counts = D19_7399)
gc()

# Normalize RNA data with log normalization
D19_7399 <- NormalizeData(D19_7399)
# Find and scale variable features
D19_7399 <- FindVariableFeatures(D19_7399, selection.method = "mean.var.plot")
D19_7399 <- ScaleData(D19_7399, features = VariableFeatures(D19_7399))

gc()

# Add HTO data as a new assay independent from RNA
D19_7399[["HTO"]] <- CreateAssayObject(counts = D19_7367)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7399 <- NormalizeData(D19_7399, assay = "HTO", normalization.method = "CLR")

D19_7399 <- HTODemux(D19_7399, assay = "HTO", positive.quantile = 0.99)

table(D19_7399$HTO_classification.global)

Idents(D19_7399) <- "HTO_classification.global"

D19_7399 <- subset(D19_7399, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7367, joint.bcs); gc()

Idents(D19_7399) <- "hash.ID"

D19_7399 <- subset(D19_7399, 
                   idents = c("SZ20-TGACGCCG", "SZ19-TGGTGTCA", "SZ21-GCCTAGTA"),
                   invert = TRUE); gc()

table(D19_7399$HTO_classification)



D19_7400 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7400/outs/filtered_feature_bc_matrix/")

D19_7368 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7368_S5/umi_count/",
                    gene.column = 1) 

colnames(D19_7368) <- paste0(colnames(D19_7368), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7400), colnames(D19_7368))

D19_7400 <- D19_7400[, joint.bcs]
D19_7368 <- as.matrix(D19_7368[, joint.bcs])

# Setup Seurat object
D19_7400 <- CreateSeuratObject(counts = D19_7400)
gc()

# Normalize RNA data with log normalization
D19_7400 <- NormalizeData(D19_7400)
# Find and scale variable features
D19_7400 <- FindVariableFeatures(D19_7400, selection.method = "mean.var.plot")
D19_7400 <- ScaleData(D19_7400, features = VariableFeatures(D19_7400))

gc()

# Add HTO data as a new assay independent from RNA
D19_7400[["HTO"]] <- CreateAssayObject(counts = D19_7368)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7400 <- NormalizeData(D19_7400, assay = "HTO", normalization.method = "CLR")

D19_7400 <- HTODemux(D19_7400, assay = "HTO", positive.quantile = 0.99)

table(D19_7400$HTO_classification.global)

Idents(D19_7400) <- "HTO_classification.global"

D19_7400 <- subset(D19_7400, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7368, joint.bcs); gc()

Idents(D19_7400) <- "hash.ID"

D19_7400 <- subset(D19_7400, 
                   idents = c("SZ20-TGACGCCG", "SZ19-TGGTGTCA", "SZ21-GCCTAGTA"),
                   invert = TRUE); gc()

table(D19_7400$HTO_classification)


D19_7401 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7401/outs/filtered_feature_bc_matrix/")

D19_7369 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7369_S6/umi_count/",
                    gene.column = 1) 

colnames(D19_7369) <- paste0(colnames(D19_7369), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7401), colnames(D19_7369))

D19_7401 <- D19_7401[, joint.bcs]
D19_7369 <- as.matrix(D19_7369[, joint.bcs])

# Setup Seurat object
D19_7401 <- CreateSeuratObject(counts = D19_7401)
gc()

# Normalize RNA data with log normalization
D19_7401 <- NormalizeData(D19_7401)
# Find and scale variable features
D19_7401 <- FindVariableFeatures(D19_7401, selection.method = "mean.var.plot")
D19_7401 <- ScaleData(D19_7401, features = VariableFeatures(D19_7401))

gc()

# Add HTO data as a new assay independent from RNA
D19_7401[["HTO"]] <- CreateAssayObject(counts = D19_7369)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7401 <- NormalizeData(D19_7401, assay = "HTO", normalization.method = "CLR")

D19_7401 <- HTODemux(D19_7401, assay = "HTO", positive.quantile = 0.99)

table(D19_7401$HTO_classification.global)

Idents(D19_7401) <- "HTO_classification.global"

D19_7401 <- subset(D19_7401, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7369, joint.bcs); gc()

Idents(D19_7401) <- "hash.ID"

D19_7401 <- subset(D19_7401, 
                   idents = c("SZ20-TGACGCCG", "SZ19-TGGTGTCA", "SZ21-GCCTAGTA", 
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7401$HTO_classification)



D19_7402 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7402/outs/filtered_feature_bc_matrix/")

D19_7370 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7370_S7/umi_count/",
                    gene.column = 1) 

colnames(D19_7370) <- paste0(colnames(D19_7370), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7402), colnames(D19_7370))

D19_7402 <- D19_7402[, joint.bcs]
D19_7370 <- as.matrix(D19_7370[, joint.bcs])

# Setup Seurat object
D19_7402 <- CreateSeuratObject(counts = D19_7402)
gc()

# Normalize RNA data with log normalization
D19_7402 <- NormalizeData(D19_7402)
# Find and scale variable features
D19_7402 <- FindVariableFeatures(D19_7402, selection.method = "mean.var.plot")
D19_7402 <- ScaleData(D19_7402, features = VariableFeatures(D19_7402))

gc()

# Add HTO data as a new assay independent from RNA
D19_7402[["HTO"]] <- CreateAssayObject(counts = D19_7370)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7402 <- NormalizeData(D19_7402, assay = "HTO", normalization.method = "CLR")

D19_7402 <- HTODemux(D19_7402, assay = "HTO", positive.quantile = 0.99)

table(D19_7402$HTO_classification.global)

Idents(D19_7402) <- "HTO_classification.global"

D19_7402 <- subset(D19_7402, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7370, joint.bcs); gc()

Idents(D19_7402) <- "hash.ID"

D19_7402 <- subset(D19_7402, 
                   idents = c("SZ20-TGACGCCG", "SZ19-TGGTGTCA", "SZ21-GCCTAGTA"),
                   invert = TRUE); gc()

table(D19_7402$HTO_classification)



D19_7403 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7403/outs/filtered_feature_bc_matrix/")

D19_7371 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7371_S8/umi_count/",
                    gene.column = 1) 

colnames(D19_7371) <- paste0(colnames(D19_7371), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7403), colnames(D19_7371))

D19_7403 <- D19_7403[, joint.bcs]
D19_7371 <- as.matrix(D19_7371[, joint.bcs])

# Setup Seurat object
D19_7403 <- CreateSeuratObject(counts = D19_7403)
gc()

# Normalize RNA data with log normalization
D19_7403 <- NormalizeData(D19_7403)
# Find and scale variable features
D19_7403 <- FindVariableFeatures(D19_7403, selection.method = "mean.var.plot")
D19_7403 <- ScaleData(D19_7403, features = VariableFeatures(D19_7403))

gc()

# Add HTO data as a new assay independent from RNA
D19_7403[["HTO"]] <- CreateAssayObject(counts = D19_7371)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7403 <- NormalizeData(D19_7403, assay = "HTO", normalization.method = "CLR")

D19_7403 <- HTODemux(D19_7403, assay = "HTO", positive.quantile = 0.99)

table(D19_7403$HTO_classification.global)

Idents(D19_7403) <- "HTO_classification.global"

D19_7403 <- subset(D19_7403, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7371, joint.bcs); gc()

Idents(D19_7403) <- "hash.ID"

D19_7403 <- subset(D19_7403, 
                   idents = c("SZ20-TGACGCCG", "SZ19-TGGTGTCA", "SZ21-GCCTAGTA"),
                   invert = TRUE); gc()

table(D19_7403$HTO_classification)


# merge all batch 7 seurat objects


batch7 <- merge(D19_7396, y = c(D19_7397, D19_7398, D19_7399,
                                D19_7400, D19_7401, D19_7402, D19_7403),
                project = "Batch7")

rm(D19_7396, D19_7397, D19_7398, D19_7399,
   D19_7400, D19_7401, D19_7402, D19_7403); gc()

batch7[["percent.mt"]] <- PercentageFeatureSet(batch7, pattern = "^MT-"); gc()
# Visualize QC metrics as a violin plot
VlnPlot(batch7,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,group.by = "orig.ident")

batch7 <- subset(batch7, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 15); gc()


saveRDS(batch7, "/media/leandro/SAMSUNG/FASTQs/batch7.rds")



################################### Batch 8 #####################################

D19_7404 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7404/outs/filtered_feature_bc_matrix/")

D19_7372 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7372_S9/umi_count/",
                    gene.column = 1) 

colnames(D19_7372) <- paste0(colnames(D19_7372), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7404), colnames(D19_7372))

D19_7404 <- D19_7404[, joint.bcs]
D19_7372 <- as.matrix(D19_7372[, joint.bcs])

# Setup Seurat object
D19_7404 <- CreateSeuratObject(counts = D19_7404,
                               min.features = 300)
gc()

# Normalize RNA data with log normalization
D19_7404 <- NormalizeData(D19_7404)
# Find and scale variable features
D19_7404 <- FindVariableFeatures(D19_7404, selection.method = "mean.var.plot")
D19_7404 <- ScaleData(D19_7404, features = VariableFeatures(D19_7404))

gc()

# Add HTO data as a new assay independent from RNA
D19_7404[["HTO"]] <- CreateAssayObject(counts = D19_7372[, colnames(D19_7404@assays$RNA$counts)])
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7404 <- NormalizeData(D19_7404, assay = "HTO", normalization.method = "CLR")

D19_7404 <- HTODemux(D19_7404, assay = "HTO", positive.quantile = 0.99)

table(D19_7404$HTO_classification.global)

Idents(D19_7404) <- "HTO_classification.global"

D19_7404 <- subset(D19_7404, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7372, joint.bcs); gc()

Idents(D19_7404) <- "hash.ID"
table(Idents(D19_7404))

D19_7404 <- subset(D19_7404, 
                   idents = c("SZ23-CTTATCAC", "SZ24-GCCTAGTA", "SZ22-CCGTACCT", 
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7404$HTO_classification)



D19_7405 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7405/outs/filtered_feature_bc_matrix/")

D19_7373 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7373_S10/umi_count/",
                    gene.column = 1) 

colnames(D19_7373) <- paste0(colnames(D19_7373), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7405), colnames(D19_7373))

D19_7405 <- D19_7405[, joint.bcs]
D19_7373 <- as.matrix(D19_7373[, joint.bcs])

# Setup Seurat object
D19_7405 <- CreateSeuratObject(counts = D19_7405)
gc()

# Normalize RNA data with log normalization
D19_7405 <- NormalizeData(D19_7405)
# Find and scale variable features
D19_7405 <- FindVariableFeatures(D19_7405, selection.method = "mean.var.plot")
D19_7405 <- ScaleData(D19_7405, features = VariableFeatures(D19_7405))

gc()

# Add HTO data as a new assay independent from RNA
D19_7405[["HTO"]] <- CreateAssayObject(counts = D19_7373)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7405 <- NormalizeData(D19_7405, assay = "HTO", normalization.method = "CLR")

D19_7405 <- HTODemux(D19_7405, assay = "HTO", positive.quantile = 0.99)

table(D19_7405$HTO_classification.global)

Idents(D19_7405) <- "HTO_classification.global"

D19_7405 <- subset(D19_7405, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7373, joint.bcs); gc()

Idents(D19_7405) <- "hash.ID"

D19_7405 <- subset(D19_7405, 
                   idents = c("SZ23-CTTATCAC", "SZ24-GCCTAGTA", "SZ22-CCGTACCT"),
                   invert = TRUE); gc()

table(D19_7405$HTO_classification)



D19_7406 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7406/outs/filtered_feature_bc_matrix/")

D19_7374 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7374_S11/umi_count/",
                    gene.column = 1) 

colnames(D19_7374) <- paste0(colnames(D19_7374), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7406), colnames(D19_7374))

D19_7406 <- D19_7406[, joint.bcs]
D19_7374 <- as.matrix(D19_7374[, joint.bcs])

# Setup Seurat object
D19_7406 <- CreateSeuratObject(counts = D19_7406)
gc()

# Normalize RNA data with log normalization
D19_7406 <- NormalizeData(D19_7406)
# Find and scale variable features
D19_7406 <- FindVariableFeatures(D19_7406, selection.method = "mean.var.plot")
D19_7406 <- ScaleData(D19_7406, features = VariableFeatures(D19_7406))

gc()

# Add HTO data as a new assay independent from RNA
D19_7406[["HTO"]] <- CreateAssayObject(counts = D19_7374)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7406 <- NormalizeData(D19_7406, assay = "HTO", normalization.method = "CLR")

D19_7406 <- HTODemux(D19_7406, assay = "HTO", positive.quantile = 0.99)

table(D19_7406$HTO_classification.global)

Idents(D19_7406) <- "HTO_classification.global"

D19_7406 <- subset(D19_7406, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7374, joint.bcs); gc()

Idents(D19_7406) <- "hash.ID"
table(Idents(D19_7406))

D19_7406 <- subset(D19_7406, 
                   idents = c("SZ23-CTTATCAC", "SZ24-GCCTAGTA", "SZ22-CCGTACCT",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7406$HTO_classification)



D19_7407 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7407/outs/filtered_feature_bc_matrix/")

D19_7375 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7375_S12/umi_count/",
                    gene.column = 1) 

colnames(D19_7375) <- paste0(colnames(D19_7375), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7407), colnames(D19_7375))

D19_7407 <- D19_7407[, joint.bcs]
D19_7375 <- as.matrix(D19_7375[, joint.bcs])

# Setup Seurat object
D19_7407 <- CreateSeuratObject(counts = D19_7407)
gc()

# Normalize RNA data with log normalization
D19_7407 <- NormalizeData(D19_7407)
# Find and scale variable features
D19_7407 <- FindVariableFeatures(D19_7407, selection.method = "mean.var.plot")
D19_7407 <- ScaleData(D19_7407, features = VariableFeatures(D19_7407))

gc()

# Add HTO data as a new assay independent from RNA
D19_7407[["HTO"]] <- CreateAssayObject(counts = D19_7375)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7407 <- NormalizeData(D19_7407, assay = "HTO", normalization.method = "CLR")

D19_7407 <- HTODemux(D19_7407, assay = "HTO", positive.quantile = 0.99)

table(D19_7407$HTO_classification.global)

Idents(D19_7407) <- "HTO_classification.global"

D19_7407 <- subset(D19_7407, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7375, joint.bcs); gc()

Idents(D19_7407) <- "hash.ID"
table(Idents(D19_7407))

D19_7407 <- subset(D19_7407, 
                   idents = c("SZ23-CTTATCAC", "SZ24-GCCTAGTA", "SZ22-CCGTACCT",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7407$HTO_classification)



D19_7408 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7408/outs/filtered_feature_bc_matrix/")

D19_7376 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7376_S13/umi_count/",
                    gene.column = 1) 

colnames(D19_7376) <- paste0(colnames(D19_7376), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7408), colnames(D19_7376))

D19_7408 <- D19_7408[, joint.bcs]
D19_7376 <- as.matrix(D19_7376[, joint.bcs])

# Setup Seurat object
D19_7408 <- CreateSeuratObject(counts = D19_7408)
gc()

# Normalize RNA data with log normalization
D19_7408 <- NormalizeData(D19_7408)
# Find and scale variable features
D19_7408 <- FindVariableFeatures(D19_7408, selection.method = "mean.var.plot")
D19_7408 <- ScaleData(D19_7408, features = VariableFeatures(D19_7408))

gc()

# Add HTO data as a new assay independent from RNA
D19_7408[["HTO"]] <- CreateAssayObject(counts = D19_7376)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7408 <- NormalizeData(D19_7408, assay = "HTO", normalization.method = "CLR")

D19_7408 <- HTODemux(D19_7408, assay = "HTO", positive.quantile = 0.99)

table(D19_7408$HTO_classification.global)

Idents(D19_7408) <- "HTO_classification.global"

D19_7408 <- subset(D19_7408, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7376, joint.bcs); gc()

Idents(D19_7408) <- "hash.ID"
table(Idents(D19_7408))

D19_7408 <- subset(D19_7408, 
                   idents = c("SZ23-CTTATCAC", "SZ24-GCCTAGTA", "SZ22-CCGTACCT",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7408$HTO_classification)



D19_7409 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7409/outs/filtered_feature_bc_matrix/")

D19_7377 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7377_S14/umi_count/",
                    gene.column = 1) 

colnames(D19_7377) <- paste0(colnames(D19_7377), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7409), colnames(D19_7377))

D19_7409 <- D19_7409[, joint.bcs]
D19_7377 <- as.matrix(D19_7377[, joint.bcs])

# Setup Seurat object
D19_7409 <- CreateSeuratObject(counts = D19_7409)
gc()

# Normalize RNA data with log normalization
D19_7409 <- NormalizeData(D19_7409)
# Find and scale variable features
D19_7409 <- FindVariableFeatures(D19_7409, selection.method = "mean.var.plot")
D19_7409 <- ScaleData(D19_7409, features = VariableFeatures(D19_7409))

gc()

# Add HTO data as a new assay independent from RNA
D19_7409[["HTO"]] <- CreateAssayObject(counts = D19_7377)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7409 <- NormalizeData(D19_7409, assay = "HTO", normalization.method = "CLR")

D19_7409 <- HTODemux(D19_7409, assay = "HTO", positive.quantile = 0.99)

table(D19_7409$HTO_classification.global)

Idents(D19_7409) <- "HTO_classification.global"

D19_7409 <- subset(D19_7409, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7377, joint.bcs); gc()

Idents(D19_7409) <- "hash.ID"
table(Idents(D19_7409))

D19_7409 <- subset(D19_7409, 
                   idents = c("SZ23-CTTATCAC", "SZ24-GCCTAGTA", "SZ22-CCGTACCT",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7409$HTO_classification)



D19_7410 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7410/outs/filtered_feature_bc_matrix/")

D19_7378 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7378_S15/umi_count/",
                    gene.column = 1) 

colnames(D19_7378) <- paste0(colnames(D19_7378), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7410), colnames(D19_7378))

D19_7410 <- D19_7410[, joint.bcs]
D19_7378 <- as.matrix(D19_7378[, joint.bcs])

# Setup Seurat object
D19_7410 <- CreateSeuratObject(counts = D19_7410)
gc()

# Normalize RNA data with log normalization
D19_7410 <- NormalizeData(D19_7410)
# Find and scale variable features
D19_7410 <- FindVariableFeatures(D19_7410, selection.method = "mean.var.plot")
D19_7410 <- ScaleData(D19_7410, features = VariableFeatures(D19_7410))

gc()

# Add HTO data as a new assay independent from RNA
D19_7410[["HTO"]] <- CreateAssayObject(counts = D19_7378)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7410 <- NormalizeData(D19_7410, assay = "HTO", normalization.method = "CLR")

D19_7410 <- HTODemux(D19_7410, assay = "HTO", positive.quantile = 0.99)

table(D19_7410$HTO_classification.global)

Idents(D19_7410) <- "HTO_classification.global"

D19_7410 <- subset(D19_7410, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7378, joint.bcs); gc()

Idents(D19_7410) <- "hash.ID"
table(Idents(D19_7410))

D19_7410 <- subset(D19_7410, 
                   idents = c("SZ23-CTTATCAC", "SZ24-GCCTAGTA", "SZ22-CCGTACCT",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7410$HTO_classification)



D19_7411 <- Read10X("/media/leandro/SAMSUNG/FASTQs/D19-7411/outs/filtered_feature_bc_matrix/")

D19_7379 <- Read10X("/media/leandro/SAMSUNG/FASTQs/CITE-seq/D19-7379_S16/umi_count/",
                    gene.column = 1) 

colnames(D19_7379) <- paste0(colnames(D19_7379), "-1")

gc()

joint.bcs <- intersect(colnames(D19_7411), colnames(D19_7379))

D19_7411 <- D19_7411[, joint.bcs]
D19_7379 <- as.matrix(D19_7379[, joint.bcs])

# Setup Seurat object
D19_7411 <- CreateSeuratObject(counts = D19_7411)
gc()

# Normalize RNA data with log normalization
D19_7411 <- NormalizeData(D19_7411)
# Find and scale variable features
D19_7411 <- FindVariableFeatures(D19_7411, selection.method = "mean.var.plot")
D19_7411 <- ScaleData(D19_7411, features = VariableFeatures(D19_7411))

gc()

# Add HTO data as a new assay independent from RNA
D19_7411[["HTO"]] <- CreateAssayObject(counts = D19_7379)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
D19_7411 <- NormalizeData(D19_7411, assay = "HTO", normalization.method = "CLR")

D19_7411 <- HTODemux(D19_7411, assay = "HTO", positive.quantile = 0.99)

table(D19_7411$HTO_classification.global)

Idents(D19_7411) <- "HTO_classification.global"

D19_7411 <- subset(D19_7411, 
                   idents = c("Negative", "Doublet"),
                   invert = TRUE); gc()
rm(D19_7379, joint.bcs); gc()

Idents(D19_7411) <- "hash.ID"
table(Idents(D19_7411))

D19_7411 <- subset(D19_7411, 
                   idents = c("SZ23-CTTATCAC", "SZ24-GCCTAGTA", "SZ22-CCGTACCT",
                              "unmapped"),
                   invert = TRUE); gc()

table(D19_7411$HTO_classification)



# merge all batch 8 seurat objects


batch8 <- merge(D19_7404, y = c(D19_7405, D19_7406, D19_7407,
                                D19_7408, D19_7409, D19_7410, D19_7411),
                project = "Batch8")

rm(D19_7404, D19_7405, D19_7406, D19_7407,
   D19_7408, D19_7409, D19_7410, D19_7411); gc()

batch8[["percent.mt"]] <- PercentageFeatureSet(batch8, pattern = "^MT-"); gc()
# Visualize QC metrics as a violin plot
VlnPlot(batch8,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,group.by = "orig.ident")

batch8 <- subset(batch8, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 15); gc()


saveRDS(batch8, "/media/leandro/SAMSUNG/FASTQs/batch8.rds")
