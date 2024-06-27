library(Seurat)
library(qs)
library(scran)    # for scoreMarkers()
library(SPOTlight)
library(STdeconvolve)
library(Matrix)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

##======================= Load the data ===========================================
## spatial 
sobj_CS_spatial_ann <- qread('03_spatial_visium/saved_obj/sp_sobj_cs_AB_res0.5_v3_endo_subcluster_ann.qs')

## scRNA
sobj_CS_RNA_ann <- qread('01_scRNA/saved_obj/gma_filtered_merged_cs_res0.4_ann.qs')

##======================= Seurat integration (anchor based)================================
## Function for ploting prediction heatmap after Seurat integration

plot_prediction_heatmapt <- function(sobj){
  ## get the prediction matrix
  mtx1 <- t(as.matrix(sobj@assays$predictions$data))
  mtx1 <- mtx1[,-ncol(mtx1)]    # remove last row ('max')
  
  ## aggregate the cell type prediction for spatial spots 
  sp_cell_ann <- sobj@meta.data['cell_types']
  merged_mat <- merge(mtx1, sp_cell_ann, by="row.names")
  head(merged_mat)
  
  ## Calculate average values by CellType
  avg_values <- aggregate(. ~ cell_types, data = merged_mat[, -1], FUN = mean)
  rownames(avg_values) <- avg_values$cell_types
  avg_values <- avg_values[, -1]
  
  ## sort the rows and cols
  n_avg_values <- avg_values[order(row.names(avg_values),decreasing = T), ]
  n_avg_values <- n_avg_values[, order(colnames(n_avg_values))]
  
  col_fun = colorRamp2(c(-0, 0.8), c("lightyellow", "maroon"))
  
  htmp <- Heatmap(n_avg_values, heatmap_legend_param = list(legend_gp = gpar(fontsize = 15)),
                  show_row_dend = F, show_column_dend = F,cluster_rows = F, cluster_columns = F,
                  row_names_side = "left",row_title_side = "left", 
                  column_names_side ="bottom", column_title_side = "top", name = 'Cell fraction',
                  column_names_rot = 45, col = col_fun,
                  row_title_gp = gpar(fontsize = 20),
                  row_names_gp = gpar(fontsize = 18),
                  column_title_gp = gpar(fontsize = 20), 
                  column_names_gp = gpar(fontsize = 18),
                  row_title = "spRNA-seq",
                  column_title = "snRNA-seq")
  return(htmp)
}

### working on CS seeds #################################################
##check the single cell dataset
DimPlot(sobj_CS_RNA_ann, label = T, repel = T)

## use only 1 object
# sobj_CS_spatial_ann1 <- SplitObject(sobj_CS_spatial_ann, split.by = 'orig.ident')[['CSA']]
sobj_CS_RNA_ann1 <- SplitObject(sobj_CS_RNA_ann, split.by = 'sampleID')[['cs2']]
anchors <- FindTransferAnchors(reference = sobj_CS_RNA_ann1, 
                               query = sobj_CS_spatial_ann,
                               normalization.method = "SCT")

predictions.assay <- TransferData(anchorset = anchors, refdata = sobj_CS_RNA_ann1$cell_types, prediction.assay = TRUE,
                                  weight.reduction = sobj_CS_spatial_ann[["pca"]], dims = 1:30)
sobj_CS_spatial_ann[["predictions"]] <- predictions.assay

DefaultAssay(sobj_CS_spatial_ann) <- "predictions"
table(Idents(sobj_CS_spatial_ann))

SpatialFeaturePlot(sobj_CS_spatial_ann,
                   features = c("Em proper initials", 'Em epidermis', 'Em dividing cell',
                                'En micropylar', 'En peripheral','En chalazal','SC epidermis',
                                'SC parenchyma', 'SC vascular parenchyma', 'SC phloem', 'SC hillum vasculature',
                                'SC endothelium','SC hillum parenchyma', 'SC xylem'
                                ),
                   pt.size.factor = 1.6, ncol = 2, crop = TRUE, alpha = 0.8, images = 'slice2')
ggsave('03_spatial_visium/deconvolution/deconv_seurat_scRNA_celltype_spatial_plot_cs_LowResAnn_052924.pdf', height = 30, width = 8)
LinkedDimPlot(sobj_CS_spatial_ann)

## draw the cell fraction heatmap
htmp <- plot_prediction_heatmapt(sobj_CS_spatial_ann)

pdf('03_spatial_visium/deconvolution/deconv_seurat_cell_frac_htmap_gma_cs_LowResAnn_052924.pdf', width = 8, height = 7)
draw(htmp)
dev.off()

## alternatively, plot the fraction score following Seurat ATAC-RNA integration tutorial, which gives similar results
sobj_CS_spatial_ann$predicted.id <- GetTransferPredictions(sobj_CS_spatial_ann, score.filter = 0.5)   # need to predict the snRNA id on spatial, which is missing in spatial integration pipeline
predictions <- table(sobj_CS_spatial_ann$cell_types, sobj_CS_spatial_ann$predicted.id)
predictions <- predictions/rowSums(predictions) 
predictions <- as.data.frame(predictions)

p1 <- ggplot(predictions, aes(Var2, Var1, fill = Freq)) + geom_tile() + 
  scale_fill_gradient(name = "Fraction of cells",low = "#ffffc8", high = "#7d0025") + 
  xlab("Cell type annotation (snRNA-seq)") + ylab("Predicted cell type label (spRNA-seq)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


## Find spatially variable celltypes/clusters (from scRNA) based on the prediction socre
sobj_HP_spatial_ann <- FindSpatiallyVariableFeatures(sobj_HP_spatial_ann, 
                                                                 assay = "predictions", selection.method = "moransi",
                                                                 features = rownames(sobj_HP_spatial_ann), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(sobj_HP_spatial_ann, selection.method = "moransi"), 4)
SpatialPlot(object = sobj_HP_spatial_ann, features = top.clusters, ncol = 2)


