library(Seurat)
# library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(harmony)
library(scCustomize)
library(MERINGUE)
library(qs)
# library(future)
# plan("multicore", workers = 10)

## 1. load the data =====================
dt1 <- Load10X_Spatial("data/gma_CS2_A/outs/",  slice="slice1")
dt1$orig.ident <- "CSA"
dt2 <- Load10X_Spatial("data/gma_CS2_B/outs/", slice="slice2")
dt2$orig.ident <- "CSB"

## 2. Quality control ================
# calculate the percentage of organelle genes
dt1 <- PercentageFeatureSet(dt1, pattern = "^GlmaCp", col.name = "percent_chlp")
dt2 <- PercentageFeatureSet(dt2, pattern = "^GlmaCp", col.name = "percent_chlp")
head(dtMerge, n=3)

###plot the molecular(RNA) counts 
plot1 <- VlnPlot_scCustom(dtMerge, features = c("nCount_Spatial", "nFeature_Spatial"), pt.size = 0.1, group.by = 'orig.ident') + NoLegend() & geom_hline(yintercept = 500)
# + scale_y_continuous(breaks = seq(0, 40000, by = 2000))
plot2 <- SpatialFeaturePlot(dtMerge, features = c("nCount_Spatial", "nFeature_Spatial")) + theme(legend.position = "right")
plots = wrap_plots(plot1, plot2)
ggsave("hp_results/QC_beforeFilter_hp.10.31.23.pdf", plots, device = 'pdf', width = 12, height = 8)

### Filter 
###Select all spots with 1000 detected genes. You must judge for yourself based on your knowledge of the tissue what are appropriate filtering criteria for your dataset.
dt1_filt = dt1[, dt1$nFeature_Spatial > 50]
SpatialFeaturePlot(dt1_filt, features = c("nCount_Spatial", "nFeature_Spatial")) + theme(legend.position = "right")

dt2_filt = dt2[, dt2$nFeature_Spatial > 50]
SpatialFeaturePlot(dt2_filt, features = c("nCount_Spatial", "nFeature_Spatial")) + theme(legend.position = "right")

p1 <- VlnPlot_scCustom(dt1_filt, features = "nCount_Spatial", pt.size = 0.1, group.by = "orig.ident") + NoLegend() & geom_hline(yintercept = 500)
p2 <- SpatialFeaturePlot(dt1_filt, features = c("nCount_Spatial", "nFeature_Spatial")) + theme(legend.position = "right") 
p3 <- VlnPlot_scCustom(dt2_filt, features = "nCount_Spatial", pt.size = 0.1, group.by = "orig.ident") + NoLegend() & geom_hline(yintercept = 500)
p4 <- SpatialFeaturePlot(dt2_filt, features = c("nCount_Spatial", "nFeature_Spatial")) + theme(legend.position = "right") 

plots <- wrap_plots(p1,p2,p3,p4,nrow = 2)
ggsave('QC/QC_filterFeature50_CS_11.15.23.pdf.pdf', plots, device = 'pdf', width = 12, height = 14)

## 3. Analysis: Dimensionality reduction, clustering, and visualization====================
## SCTransform will select variable genes and normalize in one step
dt1_filt <- SCTransform(dt1_filt, assay = "Spatial", verbose = FALSE)
dt2_filt <- SCTransform(dt2_filt, assay = "Spatial", verbose = FALSE)
dtMerge <- merge(dt1_filt, y = dt2_filt, add.cell.ids = c("CSA", "CSB"), project = 'gma_CS')

## integrate different datasets
features <- SelectIntegrationFeatures(object.list = c(dt1_filt, dt2_filt), nfeatures = 3000)
VariableFeatures(dtMerge) <- features

## Dimensionality reduction and clustering
dtMerge <- dtMerge %>% 
  RunPCA(assay = "SCT", verbose = FALSE) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(verbose = FALSE, resolution = 0.50) %>% 
  RunUMAP(reduction = "pca", dims = 1:30)

pdf('clustering/clustering_without_batch_removal_CS_11.15.23.pdf', width = 12, height = 10)
p1 <- DimPlot(dtMerge, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"),label = TRUE)
p2 <- SpatialDimPlot(dtMerge, label = TRUE, label.size = 3,alpha = c(0.6,1), repel = T)
wrap_plots(p1, p2, ncol = 1)
dev.off()

### Visualize using interactive plots 
LinkedDimPlot(dtMerge, alpha = c(0.1, 0.8), image = "slice1")

## 4. Harmony integration ================================
dtMerge <- dtMerge %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, reduction.save = "harmony")

### UMAP and Clustering using harmoney, which will generate a new $umap reduction based on harmony
dtMerge <- dtMerge %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.55) %>% 
    identity()
dtMerge <- FindClusters(dtMerge,resolution = 0.65)
LinkedDimPlot(dtMerge)

head(dtMerge)
# check the clusters again 
DimPlot(object = dtMerge, pt.size = .1, group.by = c("ident", "orig.ident"))
pdf('clustering_harmoney_batch_removal_CS7.17.23.pdf', height = 12, width = 12)
p1 <- DimPlot(dtMerge, reduction = "umap", group.by = c("ident", "orig.ident"),label = TRUE)
p2 <- SpatialDimPlot(dtMerge, label = TRUE, label.size = 3,alpha = c(0.6,1))
wrap_plots(p1, p2, ncol = 1)
dev.off()

LinkedDimPlot(dtMerge)

qsave(dtMerge, '03_spatial_visium/saved_obj/sp_sobj_cs_AB_res0.5_v3_endo_subcluster.qs')


## 5. Annotate the cell types ===========================
### subclustering (if necessary)###########
sobj_CS = qread('03_spatial_visium/saved_obj/sp_sobj_cs_AB_res0.5_v3_endo_subcluster.qs')
## subset root vasculature
sobj_CS_sub = FindSubCluster(sobj_CS, graph.name = 'SCT_snn', resolution = 0.1, cluster = 'chalazal endosperm', subcluster.name = 'endosperm_sub')

LinkedDimPlot(sobj_CS_sub, alpha = c(0.1, 0.8), image = "slice2", group.by = 'endosperm_sub')

sobj_CS_sub = SetIdent(sobj_CS_sub, value = sobj_CS_sub$endosperm_sub)
sobj_CS_sub@meta.data['cell_types'] <- Idents(sobj_CS_sub)

### 5.1. by tissue histology ################

## rename clusters, save it to a new object
sobj_CS_sub_ann <- RenameIdents(sobj_CS_sub, 
                                 "0"="embryo cotyledon", "1"="inner integument", "2"="embryo vascular bundles",
                                 "3"="chalazal endosperm", "4"="outer integument/tracheid bar", "5"="outer integument parenchyma",
                                 "6"="outer integument parenchyma", "7"="seed coat epidermis", "8"="embryo axis",
                                 "9"="hilum", "10"="micropylar endosperm", "11"="seed coat vascular bundles",
                                 "12"="embryo epidermis")
## add cell name to metadata
sobj_CS_sub_ann@meta.data['cell_types'] <- Idents(sobj_CS_sub_ann)

## save annotated seurat object
qsave(sobj_CS_sub_ann, file = '03_spatial_visium/saved_obj/sp_sobj_cs_AB_res0.5_v3_endo_subcluster_ann.qs')

## plot the cluster again with cell type annotation
pdf('03_spatial_visium/02_clustering/clustering_harmoney_batch_removal_rename_cluster_CS_res0.5_052024.pdf', height = 12, width = 16)
p1 <- DimPlot(sobj_CS_sub_ann, reduction = "umap", group.by = c("ident", "orig.ident"), label.size = 4)
p2 <- SpatialDimPlot(sobj_CS_sub_ann, label.size = 3,alpha = c(0.7,1))
wrap_plots(p1, p2, ncol = 1)
dev.off()

DimPlot(sobj_CS_sub_ann, reduction = "umap", group.by = c("ident", "orig.ident"),label = T)

## plot each clusters
pdf('annotated_clusters_separate_CS.pdf', height = 12, width = 16)
SpatialDimPlot(sobj_CS_sub_ann, cells.highlight = CellsByIdentities(dtMerge_rename), facet.highlight = TRUE, images = 'slice1', ncol = 5, alpha = c(0.2, 0.7))
dev.off()

### 5.2 find the cell type markers ================================
all_de_markers.cs  <- FindAllMarkers(sobj_CS_sub_ann, recorrect_umi = FALSE, test.use='wilcox', logfc.threshold = 1, only.pos = T, min.pct = 0.1)
write.csv(all_de_markers.cs , 'sp_all_markers.CS.csv')

### 5.3 correlation with other datasets############
avg_exp = AverageExpression(sobj_CS_sub_ann, slot="data", group.by = 'cell_types', assays = 'SCT')  
df.avg_exp <- avg_exp$SCT
# save the spatial cell type average expression for correlation analysis
write.csv(df.avg_exp, 'spatial_avg_exp_CS.csv')

# plot some markers.
pdf('markers_CS.pdf', height = 18, width = 12)
SpatialFeaturePlot(object = sobj_CS_sub_ann, features = rownames(all_de_markers.cs )[8:14], alpha = c(0.1, 1), ncol = 2)
dev.off()

# plot specific markers on spatial 
mkr_lst <- c('Glyma.01G044600', 'Glyma.02G019100', 'Glyma.08G298400', 'Glyma.18G123500')
pdf('sc_cluster0_marker.pdf')
SpatialFeaturePlot(object = sobj_HS, features = mkr_lst, 
                   alpha = c(0.1, 1), ncol = 4)
dev.off()

# selected = c('Glyma.12G056300', 'Glyma.19G186100', 'Glyma.03G185900', 'Glyma.08G341000')
pdf('markers/hilum_denovo_markers_CS.pdf', height = 18, width = 12)
selected = subset(all_de_markers.cs , all_de_markers.cs$cluster=='hilum')
selected_gene = selected$gene[1:10]
SpatialFeaturePlot(object = sobj_CS_sub_ann, features = selected_gene, alpha = c(0.1, 1), ncol = 4)
dev.off()

## 6 Find Spatially Variable genes ##########
## using markvariogram and moran'si 
sobj_CS_sub_ann <- FindSpatiallyVariableFeatures(sobj_CS_sub_ann, assay = 'SCT',
                                       slot = "scale.data",
                                       features =VariableFeatures(sobj_CS_sub_ann)[1:1000], 
                                       selection.method = 'moransi')

top.features <- head(SpatiallyVariableFeatures(sobj_CS_sub_ann, selection.method = "moransi"), 6)
# SpatiallyVariableFeatures(sobj_CS_sub_ann,selection.method = "moransi")
# plot the top features
SpatialFeaturePlot(sobj_CS_sub_ann, features = top.features, ncol = 4, alpha = c(0.1, 1))
