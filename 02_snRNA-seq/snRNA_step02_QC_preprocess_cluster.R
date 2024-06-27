# ml R/4.1.2-foss-2021b
##########################################
#### QC the data & initial clustering ####
##########################################

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

#===============================Load the dataset===============================
rt3.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/R1Solo.out/GeneFull/filtered/") # new batch
rt4.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/R2Solo.out/GeneFull/filtered/") # new batch
hp1.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/hp1Solo.out/GeneFull/filtered/")
hp2.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/hp2Solo.out/GeneFull/filtered/")
gs1.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/gs1Solo.out/GeneFull/filtered/")
gs2.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/gs2Solo.out/GeneFull/filtered/")
hs1.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/hs1Solo.out/GeneFull/filtered/")
hs2.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/hs2Solo.out/GeneFull/filtered/")
cs1.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/cs1Solo.out/GeneFull/filtered/")
cs2.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/cs2Solo.out/GeneFull/filtered/")
es1.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/es1Solo.out/GeneFull/filtered/")
es2.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/es2Solo.out/GeneFull/filtered/")
ms1.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/MS1Solo.out/GeneFull/filtered/") # new batch
ms2.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/MS2Solo.out/GeneFull/filtered/") # new batch
en1.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/EN1Solo.out/GeneFull/filtered")
en2.data <- Read10X(data.dir = "/scratch/zl57208/soybean/scRNAseq/starsolo/EN2Solo.out/GeneFull/filtered")

# Initialize the Seurat object with the raw (non-normalized data).
rt3 <- CreateSeuratObject(counts = rt3.data, project = "rt3", min.cells = 50, min.features = 100)
rt4 <- CreateSeuratObject(counts = rt4.data, project = "rt4", min.cells = 50, min.features = 100)
hp1 <- CreateSeuratObject(counts = hp1.data, project = "hp1", min.cells = 50, min.features = 100)
hp2 <- CreateSeuratObject(counts = hp2.data, project = "hp2", min.cells = 50, min.features = 100)
gs1 <- CreateSeuratObject(counts = gs1.data, project = "gs1", min.cells = 50, min.features = 100)
gs2 <- CreateSeuratObject(counts = gs2.data, project = "gs2", min.cells = 50, min.features = 100)
hs1 <- CreateSeuratObject(counts = hs1.data, project = "hs1", min.cells = 50, min.features = 100)
hs2 <- CreateSeuratObject(counts = hs2.data, project = "hs2", min.cells = 50, min.features = 100)
cs1 <- CreateSeuratObject(counts = cs1.data, project = "cs1", min.cells = 50, min.features = 100)
cs2 <- CreateSeuratObject(counts = cs2.data, project = "cs2", min.cells = 50, min.features = 100)
es1 <- CreateSeuratObject(counts = es1.data, project = "es1", min.cells = 50, min.features = 100)
es2 <- CreateSeuratObject(counts = es2.data, project = "es2", min.cells = 50, min.features = 100)
ms1 <- CreateSeuratObject(counts = ms1.data, project = "ms1", min.cells = 50, min.features = 100)
ms2 <- CreateSeuratObject(counts = ms2.data, project = "ms2", min.cells = 50, min.features = 100)
en1 <- CreateSeuratObject(counts = en1.data, project = "en1", min.cells = 50, min.features = 100)
en2 <- CreateSeuratObject(counts = en2.data, project = "en2", min.cells = 50, min.features = 100)

# combine the data
# combine the Seurat objects based on the raw count matrices, erasing any previously normalized and scaled data matrices
# merge.data = TRUE, will enable the normalized data
gma_merge <- merge(lf3, y=c(lf4,lf5,rt3,rt4,
                            hp1,hp2,gs1,gs2,hs1,
                            hs2,cs1,cs2,cs3,es1,es2,
                            ms1,ms2,sam1,sam2,
                            pd1,pd2,en1,en2),
                   add.cell.ids = c("rt3","rt4",
                                    "hp1","hp2","gs1","gs2",
                                    "hs1","hs2","cs1","cs2",
                                    "es1","es2","ms1","ms2",
                                    "en1", "en2"), project= "soybean")

##=============================== QC===============================
# QC Note:
# important features:
# unique genes / Cell (too little -> broken; too many -> doublets)
# molecules / Cell
# Organelle reads%
# rule of thumb:
# A general rule of thumb when performing QC is to set thresholds for individual metrics to be as permissive as possible, and always consider the joint effects of these metrics.

## ===============================1 edit metadata for QC===============================
gma_merge[["percent.organelle"]] <- PercentageFeatureSet(gma_merge, pattern = "^Glma") / 100
gma_merge[["percent.plstd"]] <- PercentageFeatureSet(gma_merge, pattern = "^GlmaCp") / 100
gma_merge[["percent.mt"]] <- PercentageFeatureSet(gma_merge, pattern = "^GlmaxMt") / 100

# Add number of genes per UMI for each cell to metadata
# this socre shows the complexity of the RNA, Generally, we expect the novelty score to be above 0.80 for good quality cells.
gma_merge$log10GenesPerUMI <- log10(gma_merge$nFeature_RNA) / log10(gma_merge$nCount_RNA)

# Create metadata dataframe
mtdt <- gma_merge@meta.data
mtdt$cells <- rownames(mtdt)
head(mtdt)
# Rename columns
mtdt <- mtdt %>%
  dplyr::rename(sampleID = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# # Create sample column
# mtdt$sample <- NA
# mtdt$sample[which(str_detect(mtdt$cells, "^ctrl_"))] <- "ctrl"
# mtdt$sample[which(str_detect(mtdt$cells, "^stim_"))] <- "stim"

# Add mtdt back to Seurat object
gma_merge@meta.data <- mtdt
# Create .RData object to load at any time
# save(gma_merge, file="gma_merged_24sample_rawdata_sObj.RData")
qsave(gma_merge, 'gma_merged_24sample_rawdata_sObj.qs', nthreads=20)
gma_merge = qread('gma_merged_24sample_rawdata_sObj.qs', nthreads = 30)

##=============================== 2 plot the knee plot===============================
umi <- gma_merge@meta.data[order(gma_merge@meta.data$nUMI, decreasing=T),]

pdf('kneeplots_all.pdf', height = 12, width = 20)
ggplot(umi, aes(x = log10(1:nrow(umi)), y = log10(nUMI))) +
    geom_line() +
    labs(x = "Sorted Rows", y = "nUMI") +
    ggtitle("Knee Plots") +
    # theme_minimal() +
    xlim(1,NA) + ylim(2,NA) +
    theme(text = element_text(size = 20)) +
    facet_wrap(~ sampleID, nrow = 3)   # Split the plots by "sampleID" column, 2 plots per row
    # scale_x_continuous(breaks = log10(c(1, 100, 10000)), labels = c("1", "100", "1000")
    # )       # Custom x-axis scale
dev.off()

##=============================== 3 plot the features by violin plot===============================
pdf("QC_violin_preFilter2.pdf", height = 16, width = 26)
p1 <- ggplot(gma_merge@meta.data, aes(y= nUMI, x= sampleID)) + geom_violin(aes(fill=sampleID), show.legend = F, width = 1) + geom_boxplot(width = 0.08, outlier.size = 0.5) +
      ylim(NA, 5000) + theme(axis.text=element_text(size=14), axis.title = element_text(size = 16) )
p2 <- ggplot(gma_merge@meta.data, aes(y= nGene, x= sampleID)) + geom_violin(aes(fill=sampleID), show.legend = F, width = 1) + geom_boxplot(width = 0.08,outlier.size = 0.5) +
      ylim(NA, 4000) + theme(axis.text=element_text(size=14), axis.title = element_text(size = 16) )
p3 <- ggplot(gma_merge@meta.data, aes(y= percent.organelle, x= sampleID)) + geom_violin(aes(fill=sampleID), show.legend = F, width = 1.2) + geom_boxplot(width = 0.1,outlier.size = 0.5) +
      ylim(NA, 0.5) + theme(axis.text=element_text(size=14), axis.title = element_text(size = 16) )
wrap_plots(p1, p2, p3,  nrow = 3)
dev.off()


##=============================== 4 plot the feature correlation===============================

# # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
png('gene_umi_cor_preFilter.png', height = 1500, width = 2000)
gma_merge@meta.data %>%
    ggplot(aes(x=nUMI, y=nGene, color=percent.organelle)) +
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "green") +
    stat_smooth(method=lm) +
    scale_x_log10() + scale_y_log10() +
    theme_classic() +
    xlim(1,NA) +
    geom_vline(xintercept = 500) + geom_hline(yintercept = 250) +
    facet_wrap(~sampleID)
dev.off()

# ggplot(obj_ls$ms1, aes(x=log10GenesPerUMI)) + geom_histogram()

##=============================== 5 check statistics===============================
# use psych package to summarize the statistics
sum <- describeBy(gma_merge@meta.data, group = gma_merge@meta.data$sampleID, mat=TRUE)
write.csv(sum,'QC_stats_beforeFilter.csv')

##=============================== 6 Filter out low quality cells===============================
## using selected thresholds - these will change with experiment
## filter each datasets individually&differently.
obj_ls <- SplitObject(object = gma_merge, split.by = 'sampleID')

## filter cells with over 30% x Total gene # (59k) = 17,700, filter gene< median - MAD, nUMI 2x(median - MAD), organelle <0.20
## log10GenesPerUMI > 0.80 filter low complexity gene (too many UMI with too little genes, probably due to contamination)
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

# Filter
filter.gs1 = filter_outlier(obj_ls$gs1, 1, 1, 0.15)
filter.gs2 = filter_outlier(obj_ls$gs2, 1, 1, 0.15)
filter.hs1 = filter_outlier(obj_ls$hs1, 1, 1, 0.15)
filter.hs2 = filter_outlier(obj_ls$hs2, 1, 1, 0.15)
filter.cs1 = filter_outlier(obj_ls$cs1, 1, 1, 0.15)
filter.cs2 = filter_outlier(obj_ls$cs2, 1, 1, 0.15)
filter.es1 = filter_outlier(obj_ls$es1, 1, 1, 0.15)
filter.es2 = filter_outlier(obj_ls$es2, 1, 1, 0.15)
filter.ms1 = filter_outlier(obj_ls$ms1, 1, 0, 0.15)
filter.ms2 = filter_outlier(obj_ls$ms2, 1, 1, 0.15)
filter.hp1 = filter_outlier(obj_ls$hp1, 1, 1, 0.15)
filter.hp2 = filter_outlier(obj_ls$hp2, 1, 1, 0.15)
filter.rt3 = filter_outlier(obj_ls$rt3, 1, 1, 0.15)
filter.rt4 = filter_outlier(obj_ls$rt4, 1, 1, 0.15)
filter.en1 = filter_outlier(obj_ls$en1, 1, 1, 0.15)
filter.en2 = filter_outlier(obj_ls$en2, 1, 1, 0.15)

# Remove the organelle genes
filter.rt3 = filter.rt3[!grepl("^Glma", rownames(filter.rt3)),]
filter.rt4 = filter.rt4[!grepl("^Glma", rownames(filter.rt4)),]
filter.hp1 = filter.hp1[!grepl("^Glma", rownames(filter.hp1)),]
filter.hp2 = filter.hp2[!grepl("^Glma", rownames(filter.hp2)),]
filter.gs1 = filter.gs1[!grepl("^Glma", rownames(filter.gs1)),]
filter.gs2 = filter.gs2[!grepl("^Glma", rownames(filter.gs2)),]
filter.hs1 = filter.hs1[!grepl("^Glma", rownames(filter.hs1)),]
filter.hs2 = filter.hs2[!grepl("^Glma", rownames(filter.hs2)),]
filter.cs1 = filter.cs1[!grepl("^Glma", rownames(filter.cs1)),]
filter.cs2 = filter.cs2[!grepl("^Glma", rownames(filter.cs2)),]
filter.es1 = filter.es1[!grepl("^Glma", rownames(filter.es1)),]
filter.es2 = filter.es2[!grepl("^Glma", rownames(filter.es2)),]
filter.ms1 = filter.ms1[!grepl("^Glma", rownames(filter.ms1)),]
filter.ms2 = filter.ms2[!grepl("^Glma", rownames(filter.ms2)),]
filter.en1 = filter.en1[!grepl("^Glma", rownames(filter.en1)),]
filter.en2 = filter.en2[!grepl("^Glma", rownames(filter.en2)),]

filter_gma_merge <- merge(filter.lf3, y=c(filter.lf4,filter.lf5, filter.rt3, filter.rt4,
                                           filter.sd1, filter.hp1, filter.hp2, filter.gs1, filter.gs2,
                                           filter.hs1, filter.hs2, filter.cs1, filter.cs2, filter.es1,
                                           filter.es2, filter.ms1, filter.ms2,filter.sam1,filter.sam2,
                                           filter.pd1, filter.pd2,filter.en1, filter.en2), project= "soybean")

sum2 <- describeBy(filter_gma_merge@meta.data, group = filter_gma_merge@meta.data$sampleID, mat=TRUE)
write.csv(sum2,'QC_stats_filterGeneUMI.csv')

qsavem(filter.lf3, filter.lf4,filter.lf5, filter.rt3, filter.rt4,
    filter.sd1, filter.hp1, filter.hp2, filter.gs1, filter.gs2,
    filter.hs1, filter.hs2, filter.cs1, filter.cs2, filter.es1,
    filter.es2, filter.ms1, filter.ms2,filter.sam1,filter.sam2,
    filter.pd1, filter.pd2,filter.en1, filter.en2, file = "gma_all_filterGeneUMI.qs", nthreads = 20)
qload("gma_all_filterGeneUMI.qs", nthreads = 30)

##=============================== 7 Find doublets===============================
#normalize, scale, and run UMAP of each data
#this is needed to identify the doublets in each dataset
ls(pattern='filter\\.') ## check the seurat objects loaded

## process all the data
for (i in ls(pattern='filter\\.')){
    message('... Working on ', i)
    sobj = get(i)
    assign(i, sobj %>% 
             SCTransform(vars.to.regress = 'percent.organelle') %>% 
             RunPCA(verbose = F, npcs = 20) %>% 
             RunUMAP(reduction = "pca", dims = 1:20)
  )
}

#find doublets, 10X expected 0.8% doublets in 1k targeted cells (load~16000 cell). 5k cell has ~4%, 10k cell has ~8% doublets, 15k ~12%
for (i in ls(pattern='filter\\.')){
    message('... Finding Doubletes for ', i)
    sobj = get(i)                   # get the current looped object
    nCell = round(length(sobj$sampleID))  #calculate cell number in 1k unit
    nExpRate = round(nCell/1000) * 0.008
    nExp <- round(nCell * nExpRate)         # expected number of doublets  (0.8% double rate, 8 doublets in 1k cells called)
    message('...Cell number: ', nCell, ". Expected doublets: ", nExp, "(", nExpRate, ")")
    sobj <- doubletFinder_v3(sobj, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct=T)   # find doublet
    colnames(sobj@meta.data)[ncol(sobj@meta.data)] <- "DoubletFinder"   # change the colname of the last column (DoubletFinder result)
    assign(i,sobj)    # assign the processed object back to its name
}

message('... finished finding doublets')

#filter doublets
# table(filter.lf1$DoubletFinder)
filter.rt3 = filter.rt3[, filter.rt3@meta.data[, "DoubletFinder"] == "Singlet"]
filter.rt4 = filter.rt4[, filter.rt4@meta.data[, "DoubletFinder"] == "Singlet"]
filter.hp1 = filter.hp1[, filter.hp1@meta.data[, "DoubletFinder"] == "Singlet"]
filter.hp2 = filter.hp2[, filter.hp2@meta.data[, "DoubletFinder"] == "Singlet"]
filter.gs1 = filter.gs1[, filter.gs1@meta.data[, "DoubletFinder"] == "Singlet"]
filter.gs2 = filter.gs2[, filter.gs2@meta.data[, "DoubletFinder"] == "Singlet"]
filter.hs1 = filter.hs1[, filter.hs1@meta.data[, "DoubletFinder"] == "Singlet"]
filter.hs2 = filter.hs2[, filter.hs2@meta.data[, "DoubletFinder"] == "Singlet"]
filter.cs1 = filter.cs1[, filter.cs1@meta.data[, "DoubletFinder"] == "Singlet"]
filter.cs2 = filter.cs2[, filter.cs2@meta.data[, "DoubletFinder"] == "Singlet"]
filter.es1 = filter.es1[, filter.es1@meta.data[, "DoubletFinder"] == "Singlet"]
filter.es2 = filter.es2[, filter.es2@meta.data[, "DoubletFinder"] == "Singlet"]
filter.ms1 = filter.ms1[, filter.ms1@meta.data[, "DoubletFinder"] == "Singlet"]
filter.ms2 = filter.ms2[, filter.ms2@meta.data[, "DoubletFinder"] == "Singlet"]
filter.en1 = filter.en1[, filter.en1@meta.data[, "DoubletFinder"] == "Singlet"]
filter.en2 = filter.en2[, filter.en2@meta.data[, "DoubletFinder"] == "Singlet"]

# merge_back for qc again
filter_gma_merge2 <- merge(filter.lf3, y=c(filter.lf4,filter.lf5, filter.rt3, filter.rt4,
                                           filter.sd1, filter.hp1, filter.hp2, filter.gs1, filter.gs2,
                                           filter.hs1, filter.hs2, filter.cs1, filter.cs2, filter.es1, 
                                           filter.es2, filter.ms1, filter.ms2,filter.sam1,filter.sam2,
                                           filter.pd1, filter.pd2,filter.en1, filter.en2), project= "soybean")
# check how many genes left
table(filter_gma_merge2$sampleID)

sum3 <- describeBy(filter_gma_merge2@meta.data, group = filter_gma_merge2@meta.data$sampleID, mat=TRUE)
write.csv(sum3,'QC_stats_filterGene_doubletRm.csv')

pdf("QC_violin_filterGene_doubletRm2.pdf", height = 16, width = 26)
p1 <- ggplot(filter_gma_merge2@meta.data, aes(y= nUMI, x= sampleID)) + geom_violin(aes(fill=sampleID), show.legend = F, width = 1) + geom_boxplot(width = 0.08, outlier.size = 0.5) +
  ylim(NA, 6000) + theme(axis.text=element_text(size=14), axis.title = element_text(size = 16) )
p2 <- ggplot(filter_gma_merge2@meta.data, aes(y= nGene, x= sampleID)) + geom_violin(aes(fill=sampleID), show.legend = F, width = 1) + geom_boxplot(width = 0.08,outlier.size = 0.5) +
  ylim(NA, 5000) + theme(axis.text=element_text(size=14), axis.title = element_text(size = 16) )
p3 <- ggplot(filter_gma_merge2@meta.data, aes(y= percent.organelle, x= sampleID)) + geom_violin(aes(fill=sampleID), show.legend = F, width = 1.2) + geom_boxplot(width = 0.1,outlier.size = 0.5) +
  ylim(NA, 0.2) + theme(axis.text=element_text(size=14), axis.title = element_text(size = 16) )
wrap_plots(p1, p2, p3,  nrow = 3)
dev.off()

#=============================== Clustering & Annotation===============================
# ## load the saved filtered data (merged) if the previously processed objects are not saved
qload('gma_all_filterGeneUMI_doubletRm.qs', nthread=20)
# 
##=============================== 1 Harmony integration===============================
## create the list for finding common variable features

rt_list <- list(filter.rt3, filter.rt4)
hp_list <- list(filter.hp1, filter.hp2)
gs_list <- list(filter.gs1, filter.gs2)
hs_list <- list(filter.hs1, filter.hs2)
cs_list <- list(filter.cs1, filter.cs2)
es_list <- list(filter.es1, filter.es2)
ms_list <- list(filter.ms1 )
en_list <- list(filter.en1, filter.en2)
## merge the same tissue datasets for harmony integration
filter_rt_merge = merge(filter.rt3, y=filter.rt4, project= "soybean")
filter_hp_merge = merge(filter.hp1, y=filter.hp2, project= "soybean")
filter_gs_merge = merge(filter.gs1, y=filter.gs2, project= "soybean")
filter_hs_merge = merge(filter.hs1, y=filter.hs2, project= "soybean")
filter_cs_merge = merge(filter.cs1, y=filter.cs2, project= "soybean")
filter_es_merge = merge(filter.es1, y=filter.es2, project= "soybean")
filter_ms_merge = merge(filter.ms1, project= "soybean")
filter_sam_merge = merge(filter.sam1, y=filter.sam2, project= "soybean")
filter_pd_merge = merge(filter.pd1, y=filter.pd2, project= "soybean")
filter_en_merge = merge(filter.en1, y=filter.en2, project= "soybean")

# # method1: run SCT on individual, find common feature, merge, add common feature to merge data, run normal PCA, run harmony
## find the highly variable gene present in both datasets
features <- SelectIntegrationFeatures(object.list = rt_list, nfeatures = 3000)
VariableFeatures(filter_rt_merge) <- features
features <- SelectIntegrationFeatures(object.list = hp_list, nfeatures = 3000)
VariableFeatures(filter_hp_merge) <- features
features <- SelectIntegrationFeatures(object.list = gs_list, nfeatures = 3000)
VariableFeatures(filter_gs_merge) <- features
features <- SelectIntegrationFeatures(object.list = hs_list, nfeatures = 3000)
VariableFeatures(filter_hs_merge) <- features
features <- SelectIntegrationFeatures(object.list = cs_list, nfeatures = 3000)
VariableFeatures(filter_cs_merge) <- features
features <- SelectIntegrationFeatures(object.list = es_list, nfeatures = 3000)
VariableFeatures(filter_es_merge) <- features
features <- SelectIntegrationFeatures(object.list = ms_list, nfeatures = 3000)
VariableFeatures(filter_ms_merge) <- features
features <- SelectIntegrationFeatures(object.list = sam_list, nfeatures = 3000)
features <- SelectIntegrationFeatures(object.list = en_list, nfeatures = 3000)
VariableFeatures(filter_en_merge) <- features
message('... Finished finding anchor features for merging datasets')

## run the integration
for (i in ls(pattern='filter.*merge')){
  message('... Merging dataset ', i)
  sobj = get(i)                   # get the current looped object
  sobj = sobj %>% 
    RunPCA(npcs = 20, verbose = FALSE) %>%
    RunHarmony(group.by.vars = "sampleID", assay.use = "SCT")
  message('... finished harmony - ', i)
  sobj <- sobj %>%
    RunUMAP(assay = "SCT",reduction = "harmony", dims = 1:20) %>%
    FindNeighbors(assay = "SCT",reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 0.4)
  message('... finished clustering - ', i)
  assign(i,sobj)    # assign the processed object back to its name
}

##=============================== 2 plot the clusters===============================
pdf('FeaturePlotOnUMAP_scCus.pdf', width = 14)
FeaturePlot_scCustom(filter_rt_merge, features=c('nUMI', 'percent.organelle'))
FeaturePlot_scCustom(filter_hp_merge, features=c('nUMI', 'percent.organelle'))
FeaturePlot_scCustom(filter_gs_merge, features=c('nUMI', 'percent.organelle'))
FeaturePlot_scCustom(filter_hs_merge, features=c('nUMI', 'percent.organelle'))
FeaturePlot_scCustom(filter_cs_merge, features=c('nUMI', 'percent.organelle'))
FeaturePlot_scCustom(filter_es_merge, features=c('nUMI', 'percent.organelle'))
FeaturePlot_scCustom(filter_ms_merge, features=c('nUMI', 'percent.organelle'))
FeaturePlot_scCustom(filter_en_merge, features=c('nUMI', 'percent.organelle'))
dev.off()

#merge all figure in 1
pdf('clustering_all_merged_res0.4_10.19.23.pdf', width = 20, height = 22)
p1 = DimPlot_scCustom(filter_rt_merge, reduction = "umap", label = TRUE, group.by = c("seurat_clusters"), figure_plot = T, ggplot_default_colors = T) + ggtitle('rt')
p2 = DimPlot_scCustom(filter_hp_merge, reduction = "umap", label = TRUE, group.by = c("seurat_clusters"), figure_plot = T, ggplot_default_colors = T) + ggtitle('hp')
p3 = DimPlot_scCustom(filter_gs_merge, reduction = "umap", label = TRUE, group.by = c("seurat_clusters"), figure_plot = T, ggplot_default_colors = T) + ggtitle('gs')
p4 = DimPlot_scCustom(filter_hs_merge, reduction = "umap", label = TRUE, group.by = c("seurat_clusters"), figure_plot = T, ggplot_default_colors = T) + ggtitle('hs')
p5 = DimPlot_scCustom(filter_cs_merge, reduction = "umap", label = TRUE, group.by = c("seurat_clusters"), figure_plot = T, ggplot_default_colors = T) + ggtitle('cs')
p6 = DimPlot_scCustom(filter_es_merge, reduction = "umap", label = TRUE, group.by = c("seurat_clusters"), figure_plot = T, ggplot_default_colors = T) + ggtitle('es')
p7 = DimPlot_scCustom(filter_ms_merge, reduction = "umap", label = TRUE, group.by = c("seurat_clusters"), figure_plot = T, ggplot_default_colors = T) + ggtitle('ms')
p8 = DimPlot_scCustom(filter_en_merge, reduction = "umap", label = TRUE, group.by = c("seurat_clusters"), figure_plot = T, ggplot_default_colors = T) + ggtitle('en')
wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, ncol=3)
dev.off()

## plot the the cell number by cluster
for (i in ls(pattern='filter.*merge')){
    message('... Plot cell # for dataset ', i)
    temp_obj = get(i)
    ### plot the cell#/cluster
    cluster_freq = as.data.frame(table(temp_obj@meta.data$seurat_clusters,temp_obj@meta.data$sampleID))
    bp_name = paste0('cluster/barplot_cell_per_cluster.', i, '.pdf')
    p1 = ggplot(cluster_freq, aes(x = Var1, y = Freq,fill = Var2)) + 
      geom_bar(stat = "identity", position = "dodge") +
      labs(
          title = paste0("Cells Per Cluster: ", i),
          x = "Cluster",
          y = "Number of Cells") +
      theme_minimal() +
      theme(text = element_text())
    ggsave(bp_name,p1, device = pdf)
}

# plot for single smaples, split by replicates
cluster_freq = as.data.frame(table(filter_es_merge@meta.data$seurat_clusters,filter_es_merge@meta.data$sampleID))
p1 = ggplot(cluster_freq, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = paste0("Cells Per Cluster"),
    x = "Cluster",
    y = "Number of Cells") +
  theme_minimal() +
  theme(text = element_text())
ggsave('barplot_cell_per_cluster.pdf',p1, device = pdf)


###2.3 save the clustering stats
cluster_stats = Cluster_Stats_All_Samples(filter_ms_merge, group_by_var = 'sampleID')
write.csv(cluster_stats, "clustering_stats.filter_ms_merge.csv")

# save the obj 


