library(Seurat)

setwd("~/")

batches <- readRDS("~/integrated.rds")

######################### DEGs ###############################

batches$celltype_group <- paste(batches$celltype,
                                batches$Group, sep = "_")

Idents(batches) <- "celltype_group"

oligos <- FindMarkers(batches,
                      ident.1 = "Oligo_BD",
                      ident.2 = "Oligo_Ctrl", verbose = T); gc()
oligos$gene <- rownames(oligos)
write.csv(oligos,
          "oligos_DEGs.csv")

exc <- FindMarkers(batches,
                      ident.1 = "Exc_BD",
                      ident.2 = "Exc_Ctrl", verbose = T); gc()
exc$gene <- rownames(exc)
write.csv(exc,
          "exc_DEGs.csv")

vip <- FindMarkers(batches,
                   ident.1 = "VIP+ Inhib_BD",
                   ident.2 = "VIP+ Inhib_Ctrl", verbose = T); gc()
vip$gene <- rownames(vip)
write.table(vip,
            "vip_DEGs.csv")

pvalb <- FindMarkers(batches,
                   ident.1 = "PVALB+ Inhib_BD",
                   ident.2 = "PVALB+ Inhib_Ctrl", verbose = T); gc()
pvalb$gene <- rownames(pvalb)
write.csv(pvalb,
          "pvalb_DEGs.csv")

OPC <- FindMarkers(batches,
                     ident.1 = "OPC_BD",
                     ident.2 = "OPC_Ctrl", verbose = T); gc()
OPC$gene <- rownames(OPC)
write.csv(OPC,
          "OPC_DEGs.csv")

Astrocytes <- FindMarkers(batches,
                   ident.1 = "Astrocytes_BD",
                   ident.2 = "Astrocytes_Ctrl", verbose = T); gc()
Astrocytes$gene <- rownames(Astrocytes)
write.csv(Astrocytes,
          "astro_DEGs.csv")

npy <- FindMarkers(batches,
                          ident.1 = "NPY/SST+ Inhib_BD",
                          ident.2 = "NPY/SST+ Inhib_Ctrl", verbose = T); gc()
npy$gene <- rownames(npy)
write.csv(npy,
          "npy_DEGs.csv")

Endothelial <- FindMarkers(batches,
                   ident.1 = "Endothelial_BD",
                   ident.2 = "Endothelial_Ctrl", verbose = T); gc()
Endothelial$gene <- rownames(Endothelial)
write.csv(Endothelial,
          "endo_DEGs.csv")

lamp5 <- FindMarkers(batches,
                           ident.1 = "LAMP5+ Inhib_BD",
                           ident.2 = "LAMP5+ Inhib_Ctrl", verbose = T); gc()
lamp5$gene <- rownames(lamp5)
write.csv(lamp5,
          "lamp5_DEGs.csv")

Microglia <- FindMarkers(batches,
                     ident.1 = "Microglia_BD",
                     ident.2 = "Microglia_Ctrl", verbose = T); gc()
Microglia$gene <- rownames(Microglia)

write.csv(Microglia,
          "micro_DEGs.csv")

reln <- FindMarkers(batches,
                         ident.1 = "RELN+ Inhib_BD",
                         ident.2 = "RELN+ Inhib_Ctrl", verbose = T); gc()
reln$gene <- rownames(reln)
write.csv(reln,
          "reln_DEGs.csv")

FeaturePlot(batches, 
            features = c("TMEM176B"), 
            split.by = "Group", cols = c("grey",
                                        "red"), 
            reduction = "umap")
FeaturePlot(batches, 
            features = c("HSPA6"), 
            split.by = "Group", cols = c("grey",
                                         "red"), 
            reduction = "umap")
FeaturePlot(batches, 
            features = c("AC011586.2"), 
            split.by = "Group", cols = c("grey",
                                         "red"), 
            reduction = "umap")
FeaturePlot(batches, 
            features = c("BAG3"), 
            split.by = "Group", cols = c("grey",
                                         "red"), 
            reduction = "umap")
FeaturePlot(batches, 
            features = c("HSPA1B"), 
            split.by = "Group", cols = c("grey",
                                         "red"), 
            reduction = "umap")

FeaturePlot(batches, 
            features = c("NR4A3"), 
            split.by = "Group", cols = c("grey",
                                         "red"), 
            reduction = "umap")
FeaturePlot(batches, 
            features = c("CNR1"), 
            split.by = "Group", cols = c("grey",
                                         "red"), 
            reduction = "umap")
FeaturePlot(batches, 
            features = c("DAGLB"), 
            split.by = "Group", cols = c("grey",
                                         "red"), 
            reduction = "umap")
