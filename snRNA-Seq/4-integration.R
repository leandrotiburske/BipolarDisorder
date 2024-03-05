library(Seurat)

################################# Integration #################################

batch1 <- readRDS("batch1.rds")
batch1[["RNA"]] <- JoinLayers(batch1[["RNA"]]); gc()

batch2 <- readRDS("batch2.rds")
batch2[["RNA"]] <- JoinLayers(batch2[["RNA"]]); gc()

batch3 <- readRDS("batch3.rds")
batch3[["RNA"]] <- JoinLayers(batch3[["RNA"]]); gc()

batch4 <- readRDS("batch4.rds")
batch4[["RNA"]] <- JoinLayers(batch4[["RNA"]]); gc()

batch5 <- readRDS("batch5.rds")
batch5[["RNA"]] <- JoinLayers(batch5[["RNA"]]); gc()

batch6 <- readRDS("batch6.rds")
batch6[["RNA"]] <- JoinLayers(batch6[["RNA"]]); gc()

batch7 <- readRDS("batch7.rds")
batch7[["RNA"]] <- JoinLayers(batch7[["RNA"]]); gc()

batch8 <- readRDS("batch8.rds")
batch8[["RNA"]] <- JoinLayers(batch8[["RNA"]]); gc()

batches <- merge(batch1, y = c(batch2, batch3, batch4,
                               batch5, batch6, batch7, batch8),
                 project = "Batches")
gc()

rm(batch1, batch2, batch3, batch4,
   batch5, batch6, batch7, batch8); gc()


Idents(batches) <- "hash.ID"
table(Idents(batches))

batches <- subset(batches, 
                  idents = c("unmapped"),
                  invert = TRUE); gc()

table(Idents(batches))

# run standard anlaysis workflow
batches <- NormalizeData(batches); gc()
batches <- FindVariableFeatures(batches); gc()
batches <- ScaleData(batches); gc()
batches <- RunPCA(batches); gc()

batches <- FindNeighbors(batches,
                         dims = 1:30,
                         reduction = "pca"); gc()
batches <- FindClusters(batches,
                        resolution = 0.8,
                        cluster.name = "unintegrated_clusters"); gc()

batches <- RunUMAP(batches,
                   dims = 1:30,
                   reduction = "pca",
                   reduction.name = "umap.unintegrated"); gc()

DimPlot(batches,
        reduction = "umap.unintegrated", 
        group.by = "seurat_clusters")

batches <- IntegrateLayers(object = batches,
                           method = RPCAIntegration,
                           orig.reduction = "pca",
                           new.reduction = "integrated.rpca",
                           verbose = T); gc()

# re-join layers after integration
batches[["RNA"]] <- JoinLayers(batches[["RNA"]]); gc()

batches <- FindNeighbors(batches,
                         reduction = "integrated.rpca",
                         dims = 1:30)
batches <- FindClusters(batches,
                        resolution = 0.8); gc()

batches <- RunUMAP(batches,
                   dims = 1:30,
                   reduction = "integrated.rpca")

# Visualization
DimPlot(batches,
        reduction = "umap",
        group.by = "hash.ID", raster = F) + NoLegend()

DimPlot(batches,
        reduction = "umap.unintegrated",
        group.by = "hash.ID", raster = F) + NoLegend()

DimPlot(batches,
        reduction = "umap",
        group.by = "seurat_clusters") + NoLegend()

DimPlot(batches,
        reduction = "umap.unintegrated",
        group.by = "seurat_clusters") + NoLegend()

DimPlot(batches,
        reduction = "umap",
        group.by = "hash.ID") + NoLegend()

DimPlot(batches,
        reduction = "umap.unintegrated",
        group.by = "hash.ID") + NoLegend()

saveRDS(batches, "integrated.rds")

batches <- readRDS("integrated.rds")

batches@meta.data$Group[startsWith(batches$hash.ID, "BD")] <- "BD"
batches@meta.data$Group[startsWith(batches$hash.ID, "CON")] <- "Ctrl"
unique(batches$Group)

Idents(batches) <- "seurat_clusters"

markers <- FindAllMarkers(batches, verbose = T)

head(markers)

dir.create("markers")
setwd("markers/")

markers <- markers[markers$p_val_adj < 0.05,]
markers <- markers[markers$avg_log2FC >= 0.5,]

for(cluster in unique(markers$cluster)){
  subset <- markers[markers$cluster == cluster,]
  subset <- subset %>% arrange(desc(avg_log2FC))
  write.table(subset,
              file = paste0("c", cluster,
                            "_markers.tsv",
                            sep="\t"),row.names = F)
}

setwd("../")

############################## Annotation ##############################

# load libraries
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = batches@assays$RNA$scale.data, 
                      scaled = TRUE, 
                      gs = gs_list$gs_positive, 
                      gs2 = gs_list$gs_negative) 


# merge by cluster
cL_resutls = do.call("rbind", 
                     lapply(unique(batches@meta.data$seurat_clusters),
                            function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(batches@meta.data[batches@meta.data$seurat_clusters==cl, ])]), 
                   decreasing = !0)
  head(data.frame(cluster = cl, 
                  type = names(es.max.cl),
                  scores = es.max.cl, 
                  ncells = sum(batches@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

batches@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  batches@meta.data$customclassif[batches@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

write.csv(sctype_scores,
          "sctypesc.csv")

DimPlot(batches,
        reduction = "umap", 
        label = F,
        group.by = 'customclassif', raster = F,
        cols = c("red", "blue", "green", "purple", 
                 "orange", "cyan", "magenta", "darkgreen",
                 "yellow", "brown"))  + NoLegend()      

DimPlot(batches,
        reduction = "umap", 
        label = TRUE,
        repel = TRUE,
        group.by = 'seurat_clusters') 


# Visualize cell markers

# OPC and oligos
FeaturePlot(batches,
            features = c("PDGFRA",
                         "MBP"),
            #     split.by = "stim",
            cols = c("grey","red"),
            reduction = "umap")

# Astrocyes and neurons
FeaturePlot(batches,
features = c("GFAP", 
             "RBFOX3"),
#     split.by = "stim",
cols = c("grey","red"),
reduction = "umap")

# GABAergic and glutamatergic
FeaturePlot(batches,
            features = c("SLC6A1",
                         "SLC17A7"),
            #     split.by = "stim",
            cols = c("grey","red"),
            reduction = "umap")

# Interneurons

FeaturePlot(batches,
            features = c("NPY",
                         "SST"),
            #     split.by = "stim",
            cols = c("grey","red"),
            ncol = 2,
            reduction = "umap")
FeaturePlot(batches,
            features = c("VIP",
                         "PVALB"),
            #     split.by = "stim",
            cols = c("grey","red"),
            ncol = 2,
            reduction = "umap")

FeaturePlot(batches,
            features = c("RELN",
                         "LAMP5"),
            #     split.by = "stim",
            cols = c("grey","red"),
            ncol = 2,
            reduction = "umap")


# Microglia and endothelial
FeaturePlot(batches,
            features = c("PTPRC",
                         "VWF"),
            #     split.by = "stim",
            cols = c("grey","red"),
            reduction = "umap")


Idents(batches) <- batches$seurat_clusters
new_clusters_id <- c("Oligo",
                     "Exc",
                     "Exc",
                     "Oligo",
                     "Exc",
                     "Exc",
                     "Astrocytes",
                     "Exc",
                     "PVALB+ Inhib",
                     "VIP+ Inhib",
                     "Exc",
                     "Exc",
                     "OPC",
                     "Exc",
                     "Exc",
                     "Exc",
                     "NPY/SST+ Inhib",
                     "NPY/SST+ Inhib",
                     "Exc",
                     "RELN+ Inhib",
                     "Exc",
                     "LAMP5+ Inhib",
                     "Exc",
                     "Microglia",
                     "Exc",
                     "Exc",
                     "Microglia",
                     "LAMP5+ Inhib",
                     "Astrocytes",
                     "Astrocytes",
                     "PVALB+ Inhib",
                     "Endothelial",
                     "Astrocytes",
                     "Astrocytes",
                     "Exc",
                     "Exc",
                     "OPC",
                     "Microglia",
                     "RELN+ Inhib",
                     "Astrocytes",
                     "VIP+ Inhib",
                     "NPY/SST+ Inhib",
                     "PVALB+ Inhib",
                     "Unknown",
                     "Exc")

names(new_clusters_id) <- levels(batches)
batches <- RenameIdents(batches, new_clusters_id)


batches$celltype <- Idents(batches)

DimPlot(batches,
        reduction = "umap", 
        label = F,
        group.by = 'celltype', raster = F,
        cols = c("red", "blue", "green", "purple", 
                 "orange", "cyan", "magenta", "darkgreen",
                 "yellow", "brown", "pink", "gray")) + NoLegend()

CellChat::StackedVlnPlot(batches,
                         color.use = c("red", "blue", "green", "purple", 
                                              "orange", "cyan", "magenta", "darkgreen",
                                              "yellow", "brown", "pink", "gray"),
                        features = c("RBFOX3", "SLC17A7", "SLC6A1", "NPY", "SST",
                                     "VIP", "PVALB", "RELN", "LAMP5", "PDGFRA",
                                     "MBP", "GFAP", "PTPRC", "VWF"),
                        raster = F)

# Visualize QC metrics as a violin plot
VlnPlot(batches, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"), ncol = 3, raster = F,pt.size = 0,
        group.by = "orig.ident")

saveRDS(batches,
        "integrated.rds")
