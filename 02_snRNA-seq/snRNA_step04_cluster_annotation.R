#######################################
# correlation with LCM/sp datasets
# plot known markers
# find new markers
# cluster enrichment
#######################################

library(dplyr)
library(tidyr)
library(Seurat)
library(scCustomize)
library(ComplexHeatmap)
library(ggplot2)
library(qs)


## Define function to calculate zscore correlation
z_score <- function(data_matrix) {
  mean_values <- rowMeans(data_matrix)
  sd_values <- apply(data_matrix, 1, sd)   #1: apply function sd() on rows
  z_matrix <- (data_matrix - mean_values) / sd_values
  return(z_matrix)
}

cor_zscore <- function(df1, df2, variable_gene_pct){
  ##filter low variation gene/row
  var1 <- apply(df1, 1, mad)
  cutoff1 <- quantile(var1,probs = variable_gene_pct:1) #remove lower 50% variable genes
  filter1 <- which(apply(df1, 1, mad) > cutoff1)
  df1.filter <- df1[filter1,]

  var2 <- apply(df2, 1, mad)
  cutoff2 <- quantile(var2,probs = variable_gene_pct:1)
  filter2 <- which(apply(df2, 1, mad) > cutoff2)
  df2.filter <- df2[filter2,]
  
  # get the common genes 
  common_genes <- intersect(rownames(df1.filter), rownames(df2.filter))
  df1.common <- df1[common_genes,]
  df2.common <- df2[common_genes,]
  
  ngene= length(common_genes)
  top_pct = 1-variable_gene_pct
  message('Used top ', top_pct, ' variable gene (', ngene, ') for correlation.')
  
  # calculate z scores (across samples/cluster of each gene) for both data
  z1 <- z_score(df1.common)
  z2 <- z_score(df2.common)
  
  # colnames(z2) = c('cluster0','cluster1', 'cluster2', 'cluster3', 'cluster4', 'cluster5','cluster6', 'cluster7', 'cluster8', 'cluster9', 'cluster10','cluster11', 'cluster12', 'cluster13')
  z1.m = as.matrix(z1)
  z2.m = as.matrix(z2)
  correlation <- cor(z1.m, z2.m, method="spearman")
  correlation <- t(correlation)
  return(correlation)
}


## ============== read the single cell & marker list/other dataset.========================================================
### single cell (df1) #########
sobj_rna_cs_0.4 = qread('01_scRNA/saved_obj/gma_filtered_merged_cs_res0.4_ann.qs', nthreads=20)
DimPlot(sobj_rna_es_0.2)

sc_mtx = AverageExpression(sobj_rna_cs_0.4, slot="data", assays = 'SCT' )  #give list of expression for different assays. (normally only 1 assay "RNA")
df1 <- sc_mtx$SCT

### spatial/LCM RNAseq (df2) ########
## spatial
load('03_spatial_visium/saved_obj/CS_AB_harmony_merged_obj.RData')
df2_raw <- AverageExpression(sobj_CS_spatial_ann_sub1, slot="data", assays = 'SCT' )  #give list of expression for different assays. (normally only 1 assay "RNA")
## LCM
# df2_raw <- read.csv('01_scRNA/soyseed_expression_profile_tpm_group_mean.csv', header = T, row.names = 1,check.names = F)
df2 <- df2_raw$SCT

### ATAC-seq (df3) ########
cs_atac_mtx = as.matrix(readRDS('Gm_atlas_Cotyledon_stage_seeds_metav3_gene_sparse.rds'))   # this is the raw reads count matrix, cannot be used for ploting
## read imputed ATAC matrix (processed by 04_plot_marker_accessibility.R) for ploting
cs_atac_mtx = as.matrix(cs_rna_imp$impute.activity)

##============== 1. Correlation with other (known) data.==========================================
correlation <- cor_zscore(df1,df2, 0.8)
correlation2 <- cor_zscore(df3,df2,0.8)
#edit the labels
rownames(correlation) <- gsub('\\.',' ', rownames(correlation))
colnames(correlation) <- gsub('_',' ', colnames(correlation))
rownames(correlation2) <- gsub('\\.',' ', rownames(correlation2))
colnames(correlation2) <- gsub('_',' ', colnames(correlation2))

pdf('01_scRNA/correlation_with_other_data/CS_final_zscore_correlation_vs_SpHS_res0.7_top20var.5.20.24.pdf', width = 9, height = 8)
# pheatmap(correlation, show_rownames = T, show_colnames = T)
htmap <- Heatmap(correlation, name = 'Z-score Correlation',
                 heatmap_legend_param = list(legend_direction = 'horizontal'),
                 row_names_max_width = unit(20, "cm"),
                 column_names_rot = -30,
                 # row_title = "scRNA-seq",
                 # column_title = "spatial",
                 show_row_dend = F,
                 show_column_dend = F
                 )
draw(htmap, heatmap_legend_side = "bottom")
dev.off()

####draw 2 correlation heatmap
htmp1 <- Heatmap(correlation, heatmap_legend_param = list(legend_direction = 'horizontal'),
                  show_row_dend = F, show_column_dend = F,cluster_rows = T, 
                  row_names_side = "left",row_title_side = "left",
                  column_names_side ="bottom", column_title_side = "top", name = 'r (scRNA-seq)',
                  # column_names_rot = -30, 
                  row_title = "spatial RNA-seq",
                  column_title = "single cell RNA-seq")
# draw(htmp1)
htmp2 <- Heatmap(correlation2, heatmap_legend_param = list(legend_direction = 'horizontal'),
                  show_row_dend = F, show_column_dend = F,cluster_rows = T, 
                  row_names_side = "left",row_title_side = "left",
                  column_names_side ="bottom", column_title_side = "top",name = 'r (scATAC-seq)',
                  # column_names_rot = -30,
                  row_title = "spatial",
                  column_title = "single cell ATAC-seq")
# draw(htmp2)
pdf('figure2_spatial_heatmap_cs.pdf', height = 8, width = 10)
draw(htmp1+ htmp2,heatmap_legend_side = "bottom")
dev.off()

##============== 2. Plot marker gene==========================================
marker.lst <- read.csv('./known_markers_list.cvs')
markersID <- unique(marker.lst$id[is.na(marker.lst$id)==FALSE])

#### 2.1 loop markers by dotplot/vlnplot ####
#loop to plot all the markers sets(20)
loop = ceiling(length(markersID)/10)
startn = 1
for (i in 1:loop){
  # fig_name = paste("marker_validation/HS.sc_knwon_marker.green_markers.dotplot.cs", i, ".pdf", sep = "")
  fig_name = paste("marker_validation/HS.sc_knwon_marker.vlnplot.cs", i, ".pdf", sep = "")
  endn = i * 10
  print(paste(startn, endn))
  # DotPlot(filter_cs_merge_sub, features = markersID[startn:endn]) + coord_flip()
  VlnPlot_scCustom(filter_hs_merge, features = markersID[startn:endn])
  # FeaturePlot_scCustom(filter_cs_merge, features = markersID[startn:endn])
  ggsave(fig_name, height = 20, width = 20, limitsize = FALSE)
  startn = endn
}

#### 2.2 marker heatmap of all by z-score  ####
# calculate z-score for df2 (use all genes)
z1 <- z_score(df1)    # re-calculate z-score for all the genes (without filtering)
# extract the marker genes
z1_marker = as.data.frame(subset(z1, row.names(z1) %in% marker.lst$id))
# transfer the good.marker cell label to the matrix
z1_marker$id = row.names(z1_marker)
z1_marker_celltype = merge(z1_marker, marker.lst[,c('id','name')],on = id)
row.names(z1_marker_celltype) = z1_marker_celltype$name
z1_marker_celltype <- na.omit(subset(z1_marker_celltype, select=-c(id, name)))

## plot the marker z-score in all clusters
pdf('marker_validation/marker_heatmap_zscore_ES_rmClst0_reFilter_res0.3.pdf', height = 23, width = 7)
htmap<- Heatmap(z1_marker_celltype, name = "Marker Expression",  
                cluster_columns = TRUE,
                show_column_dend = TRUE,
                cluster_column_slices = TRUE,
                cluster_rows = TRUE,
                show_row_dend = TRUE,
                # col = col_fun,
                row_names_gp = gpar(fontsize = 6),
                column_title_rot = 90,
                heatmap_legend_param = list(direction = "horizontal"),
                use_raster = TRUE,
                row_names_rot = -15,
                raster_quality = 4)
draw(htmap, heatmap_legend_side = "bottom")
dev.off()


### plot specific markers on clusters
pdf('marker_validation/CS.sp_denovo_marker.em_cotyledon.pdf',height = 10, width = 10)
FeaturePlot(sobj_rna_cs_0.4, features = c('Glyma.13G356000', 'Glyma.03G185900', 'Glyma.08G341000'))
dev.off()


##============== 3. Annotate the clusters by GO enrichment ==========================================
require(ggplot2)
require(clusterProfiler)

# read GO info, 
go2gene <- read.csv("~/OneDrive - University of Georgia/1.Research/my_database/soybean_GO_KO_annotation/go2gene_gma.csv", header = FALSE)
go2name <- read.csv("~/OneDrive - University of Georgia/1.Research/my_database/soybean_GO_KO_annotation/go2term.csv", header = TRUE)
names(go2gene) <- c("ID","Gene")
names(go2name) <- c("ID","Description","Ontology")
go_anno = merge(go2gene, go2name, by = "ID")

# subset to only analyze the bilogical process
go_anno_BP <- subset(go_anno, Ontology%in%'biological_process')

#read genes, technically not all genes because findmarkers filtered genes showing less 0.1 FC with other clusters

#### 3.1 GO enrichment analysis ####
# enrichr for GO analysis
#### define a function to do enrichment for each genes list (cluster gene)
## provide the marker table without filtering
enrichment_analysis <- function(marker_table) {
  # Get unique cluster numbers
  cluster_names <- unique(marker_table$cluster)
  
  for (cluster_name in cluster_names) {
    # Read the gene id (rownames)
    message('Working on cluster: ', cluster_name)
    
    sample_name <- my_string <- deparse(substitute(marker_table))
    markers_cluster <- subset(marker_table, cluster == cluster_name&p_val_adj <0.05)
    markers_cluster <- gsub('ann1.', '', rownames(markers_cluster))
    # print(head(markers_cluster))
    message('...DEGs in cluster ',cluster_name, ': ', length(markers_cluster))
    
    # Enrichr for GO analysis
    go_enrich <- enricher(gene = markers_cluster, TERM2GENE = go_anno[c('ID', 'Gene')], TERM2NAME = go_anno[c('ID', 'Description')], pvalueCutoff = 0.05, pAdjustMethod = 'BH')
    message('finished enricher')
    
    # Output result
    print(head(go_enrich@result))
    go_enrich@result <- merge(go_enrich@result, go2name[c('ID', 'Ontology')], by = "ID", sort = FALSE)
    file_name <- paste("GO_enrichment", sample_name,'.', cluster_name, ".csv", sep = "")
    print(file_name)
    write.csv(go_enrich@result, file_name, row.names = FALSE)
    
    # Plot the result
    pic_name <- paste("GO_enrichment_dotplot", sample_name,'.', cluster_name,".pdf", sep = "_")
    dp.enrich <- dotplot(go_enrich, showCategory = 8, title = paste("Cluster", cluster_name, "genes"))
    ggsave(dp.enrich, file = pic_name)
    
    #### GSEA ####
    # prepare the data
    sub_marker <- subset(marker_table, cluster == cluster_name)  
    glist <- sub_marker[, "avg_log2FC"]
    names(glist) <- as.character(gsub('ann1.', '', rownames(sub_marker)))
    glist <- sort(glist, decreasing = TRUE)
    
    message('...number of Genes used for gsea: ',length(glist))
    
    # run analysis
    gsea =GSEA(glist,TERM2GENE = go_anno[c('ID', 'Gene')], TERM2NAME = go_anno[c('ID', 'Description')], pvalueCutoff = 0.5 )
    file_name2 = paste('GSEA_result.', sample_name,'.', cluster_name, '.csv',sep = '')
    write.csv(gsea@result, file_name2)
  }
  
  return(list(go_enrichment = go_enrich, gsea_list = gsea))
}

enrich_res <- enrichment_analysis(all_spatial_markers)

##============== 4. Rename cell types ==========================================
### GS ############
table(Idents(sobj_rna_cs_0.4))
sobj_rna_cs_0.4_ann <- RenameIdents(sobj_rna_cs_0.4,
                  "0"="SC parenchyma", "1_0"="SC epidermis",'1_1'='SC epidermis', '1_2'='SC hillum parenchyma',
                  '1_3'='SC xylem',"2"="SC phloem",'3'='SC parenchyma','4'='SC hillum vasculature',
                  '5'='SC endothelium',"6"="SC inner integument", "7"="Dividing cell", "8"="Endosperm",
                   "9"="SC vascular parenchyma", "10"="Endosperm", "11"="Endosperm",
                  "12"="Emb initials","13"="Endosperm", '14'='SC unknown')
## add renamed celltype to metadata
sobj_rna_cs_0.4_ann$cell_types <- Idents(sobj_rna_cs_0.4_ann)
## save celltype annotated seurat object
qsave(sobj_GS_RNA_ann,'01_scRNA/saved_obj/gma_cs_merged_celltype_anno.qs')

