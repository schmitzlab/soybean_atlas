# Downstream analysis of the clusters / refinement 

##=============================== Load package ===============================
library(Seurat)
library(dplyr)
library(ggplot2)
library(monocle3)
library(scCustomize)
library(qs)   #quickly writing and reading any R object to and from disk.
library(psych)  ## for stats describeby()
library(DoubletFinder)
library(harmony)
library(ComplexHeatmap)


filter_outlier <- function(sobj, nfold, org) {
  lower_gene = median(sobj@meta.data$nGene) - nfold * mad(sobj@meta.data$nGene)
  lower_UMI = quantile(sobj@meta.data$nUMI, 0.15)
  filtObj = subset(x = sobj, subset= (nUMI >= lower_UMI) & (nGene >= lower_gene & nGene <= 17700)  & (percent.organelle < org) & (log10GenesPerUMI > 0.85)) 
  message('Filtering outlier cells for ', substitute(sobj))
  message('filter cell min nGene:', lower_gene)
  message('filter cell min UMI:', lower_UMI)
  message('Cells filter: ', nrow(sobj@meta.data), " ---> ", nrow(filtObj@meta.data), '\n')
  return(filtObj)
}

##=============================== 1. Load data ===============================
### processing cs seeds as example ################
qload('/saved_obj/gma_filtered_merged_seeds.res0.4.qs', nthread = 20)

##=============================== 2. Remove clusters/subset ===============================

cs_sobj_sub1 = subset(x = filter_cs_merge, subset= (seurat_clusters !=0)) # remove questionable cluster0

## re-clustering after cell removal
### the PCs/variable features may be different after removing the bad cluster
cs_sobj_sub1_1 <- cs_sobj_sub1 %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(assay = "SCT",reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(assay = "SCT",reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.2)

## adjust cluster resolution
cs_sobj_sub1_2 = FindClusters(cs_sobj_sub1_1, resolution = 0.3)

## plot after removing cells
pl = DimPlot_scCustom(cs_sobj_sub1_2, reduction = "umap", label = TRUE, group.by = c("seurat_clusters"), figure_plot = T, ggplot_default_colors = T)
ggsave('cluster/clustering_cs_rm_clst0_re_analysis.pdf', pl, device = 'pdf', height = 7)

## save the seurat object 
save(cs_sobj_sub1_2, 'saved_obj/gma_filtered_merged_cs_res0.4_ann.qs', nthreads = 20)

#check stats per cluster
p1 = QC_Plots_Genes(seurat_object = cs_sobj_sub1_2, pt.size = 0,
                    group.by = 'seurat_clusters', plot_median = TRUE, y_axis_log = T) 
p2 = QC_Plots_UMIs(seurat_object = cs_sobj_sub1_2, pt.size = 0,
                   group.by = 'seurat_clusters', plot_median = TRUE, y_axis_log = T)
plts = wrap_plots(p1,p2, nrow = 2)
ggsave('vlnplot_cell_per_cluster_cs_merge.pdf', plts, device = 'pdf', height = 7, width = 9) 


#### find subcluster (if necessary) #####
cs_sobj_sub1_3 <- FindSubCluster(cs_sobj_sub1_2, "0",graph.name = 'SCT_snn',  subcluster.name = "sub_clst0",  resolution = 0.3, algorithm = 1)
##Your subcluster is now saved in metadata. Use SetIdent to save new cluster assignment to the main Ident in your object.
cs_sobj_sub1_3 <- SetIdent(cs_sobj_sub1_3, value = cs_sobj_sub1_3@meta.data$sub_clst0)


#### find markers for adjusted cluster #########
cs_sobj_sub1_3 <- PrepSCTFindMarkers(cs_sobj_sub1_3,assay = "SCT")
cs.DEGs <- FindAllMarkers(cs_sobj_sub1_3, recorrect_umi = FALSE, test.use='wilcox', logfc.threshold = 1, only.pos = T, min.pct = 0.1)
write.csv(cs.DEGs, "all_markers_cs_merged_rmClst0_reFilter_res0.3.csv", quote = F)

##Visualize the main object without choosing meta data to see if it worked.
pdf('cluster/clustering_cs_merged_rm_clst0_re_filter_res0.3_subclust0.pdf', height = 10, width = 10)
DimPlot_scCustom(cs_sobj_sub1_3, label = TRUE, group.by = 'sub_clst0',label.size = 6, repel = T, ggplot_default_colors = T)
dev.off()

