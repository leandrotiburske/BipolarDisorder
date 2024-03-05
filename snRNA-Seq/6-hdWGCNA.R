setwd("~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/hdWGCNA/")

# single-cell analysis package
library(Seurat)
options(Seurat.object.assay.version = 'v4')

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)
library(fgsea)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 4)

# load the snRNA-seq dataset
seurat_obj <- readRDS('~/Documents/IC/BipolarDisorder/snRNAseq/PsychENCODE/hdWGCNA/seurat_astro.rds')

seurat_obj[["RNA3"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay")
DefaultAssay(seurat_obj) <- "RNA3"

DimPlot(seurat_obj, group.by='cluster', label=TRUE, reduction = "umap",) +
  umap_theme() + ggtitle('Astrocytes') + NoLegend()

DimPlot(seurat_obj, group.by='cluster', split.by = "Group", label=TRUE, reduction = "umap") +
  umap_theme() + ggtitle('Astrocytes') + NoLegend()


seurat_obj@assays$RNA <- NULL
seurat_obj@assays$HTO <- NULL
gc()


FeaturePlot(seurat_obj, features = c("HSPB1", "DNAJB1"), 
            reduction = "umap", split.by = "Group", ncol = 2, cols = c("grey",
                                                                       "red"))
FeaturePlot(seurat_obj, features = c("HSPH1", "HSP90AA1"), 
            reduction = "umap", split.by = "Group", ncol = 2, cols = c("grey",
                                                                       "red"))


FeaturePlot(seurat_obj, features = c("CHI3L1", "SERPINA3"), 
            reduction = "umap", split.by = "Group", ncol = 2, cols = c("grey",
                                                                       "red"))
FeaturePlot(seurat_obj, features = c("VIM", "GFAP"), 
            reduction = "umap", split.by = "Group", ncol = 2, cols = c("grey",
                                                                       "red"))

FeaturePlot(seurat_obj, features = c("ALDH1A1"), 
            reduction = "umap", split.by = "Group", ncol = 2, cols = c("grey",
                                                                       "red"))

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.01, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Astrocytes" # the name of the hdWGCNA experiment
); gc()

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cluster", "Group"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 20, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cluster' # set the Idents of the metacell seurat object
); gc()

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c(0:17), # the name of the group of interest in the group.by column
  group.by='cluster', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA3', # using RNA assay
  slot = 'data' # using normalized data
)

gc()

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power = 7,
  tom_name = 'Astrocytes', overwrite_tom = T # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='Astrocytes hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

saveRDS(seurat_obj, "hdWGCNA.rds")

#seurat_obj <- readRDS("hdWGCNA.rds")

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(seurat_obj, 
                               verbose = T,
                               group.by.vars = "Group"); gc()

MEs <- GetMEs(seurat_obj, harmonized=T)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj, group_name = 'ASTRO'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "ASTRO-M"
)

# plot genes ranked by kME for each module
PlotKMEs(seurat_obj)

# get the module assignment table:
MEs <- GetMEs(seurat_obj, harmonized=T)
modules <- GetModules(seurat_obj)

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)

saveRDS(seurat_obj, file='hdWGCNA.rds')

#seurat_obj <- readRDS("hdWGCNA.rds")

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=T)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
gc()
# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'cluster')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

# enrichr databases to test
dbs <- c('GO_Biological_Process_2021')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)

# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# compute cell-type marker genes with Seurat:
Idents(seurat_obj) <- seurat_obj$cluster
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1
)

write.csv(markers,
          "markers_astro.csv")

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork
wrap_plots(plot_list, ncol=4)

# plot odds ratio of the overlap as a dot plot
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = F, # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

plot_list$`ASTRO-M8`
plot_list$`ASTRO-M7`

saveRDS(seurat_obj, "hdWGCNA.rds")
#seurat_obj <- readRDS("hdWGCNA.rds")

set.seed(123)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs =5,
  n_neighbors=10,
  min_dist=0.5,
  spread=1,
  supervised=TRUE,
  target_weight=1
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.5,
  sample_edges=FALSE,vertex.label.cex = 0.25,
  edge_prop=0.075, # proportion of edges to sample (20% here)
  label_hubs=0,# how many hub genes to plot per module?
  label_genes = c("DOCK8", "ST18", "LHFPL3", "HSPH1", 
                  "OPCML", "AC012405.1", "FOS", "SLC1A2",
                  "CHI3L1", "AC008957.2", "NEAT1")
)

############### Plot animated network #########################

# different label weights to test
weights <- 0:10/10

# loop through different weights
df <- data.frame()
for(cur_weight in weights){
  set.seed(123)
  # make a module UMAP using different label weights
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 5,
    n_neighbors=10,
    exclude_grey = TRUE,
    min_dist=0.5,
    supervised=TRUE,
    target_weight = cur_weight
  )
  
  # add to ongoing dataframe
  cur_df <- GetModuleUMAP(seurat_obj)
  cur_df$weight <- cur_weight
  df <- rbind(df, cur_df)
}

# ggplot animation library
library(gganimate)

# plot with ggplot + gganimate
p <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(color=df$color, size=df$kME*2 ) +
  ggtitle("Supervised weight: {closest_state}") +
  transition_states(
    weight,
    transition_length = 2,
    state_length = 2,
    wrap = TRUE
  ) +
  view_follow() +
  enter_fade() +
  umap_theme()

animate(p, fps=60, duration=25)

library(gifski)
png_files <- list.files("./", pattern = ".*png$", full.names = TRUE)
gifski(png_files, gif_file = "animation.gif",loop = T,delay = 0.05)

# plot with ggplot
plot_df <- umap_df

library(igraph)

# individual module networks
ModuleNetworkPlot(
  seurat_obj,
  outdir = paste0('RGL_hubNetworks/')
)

# Show number of cells per cluster
table(seurat_obj@meta.data$cluster)

# Show number of cells whitin cluster 13 per group
c13 <- subset(seurat_obj, subset = cluster == 13)
table(c13@meta.data$Group)
rm(c13); gc()

# Show number of cells whitin cluster 16 per group
c16 <- subset(seurat_obj, subset = cluster == 16)
table(c16@meta.data$Group)
rm(c16); gc()

saveRDS(seurat_obj, "hdWGCNA.rds")
