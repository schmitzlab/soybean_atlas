###################################################################################################
##                  transfer scATAC with annotated snRNA-seq  with Seurat                        ##
###################################################################################################
# load libraries ----------------------------------------------------------------------------------
library(Seurat)
library(Matrix)
library(qs)
library(parallel)
library(ggplot2)
library(patchwork)
library(cowplot)
library(RColorBrewer)

# load arguments ----------------------------------------------------------------------------------
args <- commandArgs(T)
if(length(args)!=4){stop("Rscript Label_transfer_snRNA_scATAC_Seurate.R <scRNA_seurate.obj.qs> <scATAC.soc.processed.rds> <scatac.gene.sparse.rds>  <prefix>")}

rna <- as.character(args[1])
atac <- as.character(args[2])
atac.gene.rds <- as.character(args[3])
prefix <- as.character(args[4])


#---load rna---
rna <- readRDS(rna)

#updata rna
rna$tech <- "rna"
rna$celltype <- Idents(rna)

#---load atac---
atac <- readRDS(atac)
atac.gene.rds <- readRDS(atac.gene.rds)
atac.gene.rds <- as.matrix(atac.gene.rds)
row.names(atac.gene.rds) <- paste0("ann1.",row.names(atac.gene.rds))

#update socrate object...
atac$counts <- atac$counts[,colnames(atac$counts) %in% rownames(atac$h.Clusters)]
atac$counts <- atac$counts[,colnames(atac$counts) %in% colnames(atac.gene.rds)]
atac$counts <- atac$counts[Matrix::rowMeans(atac$counts)>0,]
atac$h.Clusters <- atac$h.Clusters[colnames(atac$counts),]
atac$h.UMAP <- atac$h.UMAP[colnames(atac$counts),]
atac$l2.PCA <- atac$l2.PCA[colnames(atac$counts),]
atac.gene.rds <- atac.gene.rds[,colnames(atac.gene.rds) %in% colnames(atac$counts)]

#create atac seurate obj...
atac.seurobj <- CreateSeuratObject(counts = atac$counts, assay = "ATAC", project = prefix)
atac.seurobj <- AddMetaData(atac.seurobj, metadata = atac$h.Clusters)
colnames(atac$l2.PCA) <- paste0("PC_", 1:(ncol(atac$l2.PCA)))
atac.seurobj[['pca']] <- CreateDimReducObject(embeddings=as.matrix(atac$l2.PCA), assay = "ATAC", key = "PC_")
colnames(atac$h.UMAP) <- paste0("UMAP_", 1:(ncol(atac$h.UMAP)))
atac.seurobj[['UMAP']] <- CreateDimReducObject(embeddings = as.matrix(atac$h.UMAP), key = "UMAP_", global = T, assay = "ATAC")
atac.seurobj[["ACTIVITY"]] <- CreateAssayObject(counts = atac.gene.rds)

#---Transfer labels---
#plot raw umap.
p1 <- DimPlot(atac.seurobj, reduction = "UMAP") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
ggsave(paste0(prefix,"_RNA_ATAC_raw_UMAP.pdf"), plot = p1+p2, width = 12, height = 6)

#find anchors.
transfer.anchors <- FindTransferAnchors(reference = rna, query = atac.seurobj, features = VariableFeatures(object = rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

#might need to test if ignore the first PCA or not.
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna$celltype,
  weight.reduction = atac.seurobj[['pca']],
  dims = 1:(ncol(atac$l2.PCA))
)

atac.seurobj <- AddMetaData(object = atac.seurobj, metadata = predicted.labels)

#-plot the prediction score
#summary(atac.seurobj$prediction.score.max)
high <- sum(atac.seurobj$prediction.score.max >= 0.5)
low <- sum(atac.seurobj$prediction.score.max < 0.5)
med <- round(median(atac.seurobj$prediction.score.max),digits = 3)
pdf(file=paste0(prefix,"_AllPCA.prediction.score.max.pdf"),width=4, height=4)
hist(atac.seurobj$prediction.score.max)
abline(v = 0.5, col = "red")
mtext(paste0("#High=",high,"; #Low=",low,"; Median=",med), side=3, cex = 1)
dev.off()

#-plot all cells.
atac.seurobj$predicted.id <- factor(atac.seurobj$predicted.id, levels = levels(rna))
p1 <- DimPlot(atac.seurobj, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells")
ggsave(paste0(prefix,"_RNA_ATAC_AllMatched_UMAP.pdf"), plot = p1+p2, width = 16, height = 6)

#-plot the cells with high and low prediction score.
atac.high <- subset(atac.seurobj, subset = prediction.score.max >= 0.5)
atac.high$predicted.id <- factor(atac.high$predicted.id, levels = levels(rna))

p1 <- DimPlot(atac.high, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("HighQ scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
ggsave(paste0(prefix,"_RNA_ATAC_High_UMAP.pdf"), plot = p1, width = 12, height = 6)

#-plot high quality cell-
predictions <- table(atac.high@meta.data$LouvainClusters, atac.high$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() +
  scale_fill_gradient(name = "Fraction of cells",low = "#ffffc8", high = "#7d0025") + 
  xlab("ATAC-seq LouvainClusters") + 
  ylab("Predicted RNA cell type label") + 
  ggtitle("High quality predicted cells") +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(hjust = 0.5))
ggsave(paste0(prefix,"_RNA_ATAC_HighMatched_proportion_heatmap.pdf"), plot = p1, width = 10, height = 4)

#---Output files---
write.table(atac.meta, file=paste0(prefix, ".atac.labelTran.metadata.txt"), quote=F, row.names=T, col.names=T, sep="\t")
saveRDS(atac.seurobj, file=paste0(prefix,".atac.seurat.labelTran.rds"))
