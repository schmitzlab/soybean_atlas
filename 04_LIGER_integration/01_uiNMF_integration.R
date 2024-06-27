## uinmf ##

# load libraries
library(rliger)
library(Seurat)
library(Socrates)
library(RColorBrewer)

# data
atac.g <- readRDS("Gm_atlas_seed_4_stages.embryo.gene_sparse.rds")
atac.p <- readRDS("Gm_atlas_seed_4_stages.embryo.scATAC_ACR_sparse.rds")
rna <- readRDS("Gm_atlas_seed_4_stages.embryo.snRNA_gene_sparse.rds")
atac.meta <- read.table("Gm_atlas_all_Embryo_atac.celltype.metadata.txt", comment.char="")
rna.meta <- read.table("Gm_atlas_all_Embryo_rna.celltype.metadata.txt", comment.char="")

# output
input <- "soybean_embryo"

# filter cells
atac.shared <- intersect(rownames(atac.meta), colnames(atac.g))
atac.shared <- intersect(atac.shared, colnames(atac.p))
atac.meta <- atac.meta[atac.shared,]
atac.p <- atac.p[,atac.shared]
atac.g <- atac.g[,atac.shared]
atac.p <- atac.p[Matrix::rowSums(atac.p) > 0,]
atac.g <- atac.g[Matrix::rowSums(atac.g) > 0,]
rna.match <- intersect(rownames(rna), rownames(atac.g))
atac.g <- atac.g[rna.match,]
rna <- rna[rna.match,]
atac.g <- atac.g[,Matrix::colSums(atac.g) > 0]
atac.ids <- intersect(colnames(atac.g), rownames(atac.meta))
atac.ids <- intersect(atac.ids, colnames(atac.p))
atac.meta <- atac.meta[atac.ids,]
atac.p <- atac.p[,atac.ids]

# select unshared features
soc.obj <- list(counts=atac.p, meta=atac.meta)
norm <- tfidf(soc.obj, doL2=T)$residuals
liger <- createLiger(list(peaks = norm))
liger <- rliger::normalize(liger)
liger@norm.data$peaks <- liger@raw.data$peaks %*% Diagonal(x=1/Matrix::colSums(liger@raw.data$peaks))
rm(soc.obj)

# top 5000 variable features
se <- CreateSeuratObject(norm)
vars_5000 <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
top5000 <- head(VariableFeatures(vars_5000),2000)
top5000 <- gsub("-","_",top5000)
liger <- selectGenes(liger)
liger@var.genes <- top5000
liger <- scaleNotCenter(liger)
unshared_feats = liger@scale.data$peaks

# clean-up
rm(se)
rm(norm)
rm(vars_5000)
rm(top5000)

# preprocessing
liger <- createLiger(list(rna = rna, atac = atac.g))
liger <- rliger::normalize(liger)
liger <- selectGenes(liger, 
                     var.thresh = 0.1, 
                     datasets.use =1 , 
                     unshared = TRUE,  
                     unshared.datasets = list(2), 
                     unshared.thresh= 0.2)
liger <- scaleNotCenter(liger)
peak_names <- colnames(unshared_feats)
liger@var.unshared.features[[2]] = peak_names
liger@scale.unshared.data[[2]] = t(unshared_feats)

# integrate
liger <- optimizeALS(liger, k=30, use.unshared = TRUE, max_iters =30,thresh=1e-10)
liger <- quantile_norm(liger, ref_dataset="rna")

# cluster
liger <- louvainCluster(liger)

# plot
liger <- runUMAP(liger)
umap_plots <-plotByDatasetAndCluster(liger, axis.labels = c("UMAP1","UMAP2"), return.plots = TRUE)

# plot
pdf(paste0(input,".LIGER_integration.pdf"), width=16, height=15)
umap_plots[[1]]
dev.off()

pdf(paste0(input,".LIGER_clusters.pdf"), width=16, height=15)
umap_plots[[2]]
dev.off()

# add library info
atac.libdata <- c("Duke15_Gm_Globular_stage_seeds", "Duke7_Gm_Middle_maturation_stage_seeds",
                  "Heart_stage_seeds", "NCS2_Gm_Cotyledon_stage_seeds", "NCS3_Gm_Early_maturation_stage_seeds",
                  "NCS3_Gm_Globular_stage_seeds", "NCS3_Gm_Pod")
atac.tissue <- c("Globular", "Middle_Maturation", "Heart", "Cotyledon",
                 "Early_Maturation", "Globular", "Pod")
names(atac.tissue) <- atac.libdata
rna.libdata <- c("RNA_cs1", "RNA_cs2", "RNA_es1", "RNA_es2", "RNA_gs1", "RNA_gs2", "RNA_hs1", "RNA_hs2")
rna.tissue <- c("Cotyledon", "Cotyledon", "Early_Maturation", "Early_Maturation", 
                "Globular", "Globular", "Heart", "Heart")
names(rna.tissue) <- rna.libdata
atac.meta$TissueType <- atac.tissue[atac.meta$tissue]
rna.meta$TissueType <- rna.tissue[rna.meta$library]
TissueType <- c(atac.meta$TissueType, rna.meta$TissueType)
names(TissueType) <- c(rownames(atac.meta), rownames(rna.meta))
liger@cell.data$tissue <- TissueType[rownames(liger@cell.data)]

# save output
inmf <- liger@H.norm
colnames(inmf) <- paste0("NMF_",seq(1:ncol(inmf)))
meta <- cbind(liger@cell.data, liger@clusters)
colnames(meta)[5] <- "clusters"
meta <- cbind(meta, liger@tsne.coords)
colnames(meta)[6:7] <- c("umap1", "umap2")
meta$tissue <- factor(meta$tissue, levels=c("Globular", "Heart", "Cotyledon", "Early_Maturation"))
liger@cell.data <- meta
saveRDS(liger, file="uiNMF_integration.soybean_embryo.rds")

# plot results
cols <- brewer.pal(4, "Spectral")
cols2 <- c("grey75", "dodgerblue3")
pdf("seed_stage.pdf", width=6, height=6)
plot(meta$umap1, meta$umap2, pch=16, cex=0.3, col=cols[as.numeric(meta$tissue)],
     xlab="UMAP1", ylab="UMAP2")
legend("bottomleft", legend=levels(meta$tissue), fill=cols)
dev.off()

# differential analysis
dg <- runWilcoxon(liger, data.use = 'all', compare.method = 'clusters')
saveRDS(dg, file="uiNMF_integration.cluster_wilcoxon.rds")