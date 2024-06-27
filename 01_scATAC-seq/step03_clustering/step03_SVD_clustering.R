#run Socrates on seperate socrates object #
#module load R/4.3.1-foss-2022a

# load arguments
args <- commandArgs(T)
if(length(args) != 2){stop("Rscript SVD_clustering_xz.R <tissue>")}
tissue <- as.character(args[1])



# libraries
library(Seurat)
library(Socrates)
library(harmony)
library(symphony)
library(igraph)
library(Matrix)
library(RcppML)
library(tidyverse)

#functions
############mergeSocObjects--------------------------------------------------------
mergeSocObjects <- function(obj.list){
  
  # functions
  .merge.sparse <- function(cnt.list) {
    
    cnnew <- character()
    rnnew <- character()
    x <- vector()
    i <- numeric()
    j <- numeric()
    
    for (M in cnt.list) {
      
      cnold <- colnames(M)
      rnold <- rownames(M)
      
      cnnew <- union(cnnew,cnold)
      rnnew <- union(rnnew,rnold)
      
      cindnew <- match(cnold,cnnew)
      rindnew <- match(rnold,rnnew)
      ind <- summary(M)
      i <- c(i,rindnew[ind[,1]])
      j <- c(j,cindnew[ind[,2]])
      x <- c(x,ind[,3])
    }
    
    sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
  }
  
  # separate counts and meta
  counts <- lapply(obj.list, function(x){
    x$counts
  })
  counts <- .merge.sparse(counts)
  
  # meta
  metas <- lapply(obj.list, function(x){
    x$meta
  })
  metas <- do.call(rbind, metas)
  rownames(metas) <- metas$cellID
# metas$library <- data.frame(do.call(rbind, str_split(rownames(metas), "-", n=2)))[,2]
  metas <- metas[colnames(counts),]
  
  # new object
  new.obj <- list(counts=counts, meta=metas)
  return(new.obj)
}


# load rds files and pre-processing --------------------------------------

meta <- c("Gm_Root-rep3.updated_metadata.txt","Gm_Root-rep4.updated_metadata.txt","Gm_Root-rep5.updated_metadata.txt")
rds <- c("Gm_Root-rep3.raw.soc.rds","Gm_Root-rep4.raw.soc.rds","Gm_Root-rep5.raw.soc.rds")

out <- tissue

#merge meta
meta <- lapply(meta, function(x){read.table(x)})
meta <- do.call(rbind, meta)
message("read data****")
meta$library <- data.frame(do.call(rbind, str_split(rownames(meta),"-",n=2)))[,2]
meta$tissue <- data.frame(do.call(rbind, str_split(rownames(meta),"-",n=3)))[,2]
meta$rep <- data.frame(do.call(rbind, str_split(rownames(meta),"-",n=3)))[,3]
meta <- meta[meta$call==1,] #filter bad cell or not.

message("finished read data****")

#filter meta and merge rds files.
objs <- lapply(rds, function(x){
  df <- readRDS(x)
  df$counts <- df$counts[,colnames(df$counts) %in% rownames(meta)]
  df$counts <- df$counts[Matrix::rowMeans(df$counts)>0,]
  df$meta <- df$meta[colnames(df$counts),]
  return(df)
})
soc.obj <- mergeSocObjects(objs)

meta <- meta[rownames(soc.obj$meta),]
soc.obj$meta$library <- meta$library
soc.obj$meta$tissue <- meta$tissue
soc.obj$meta$rep <- meta$rep
soc.obj$meta$pOrg <- meta$pOrg
soc.obj$meta$pTSS <- meta$pTSS
soc.obj$meta$log10nSites <- meta$log10nSites
soc.obj$meta$FRiP <- meta$FRiP

soc.obj$counts <- soc.obj$counts[Matrix::rowSums(soc.obj$counts)>0,]


# get per cell feature counts --------------------------------------------
cell.counts <- log10(Matrix::colSums(soc.obj$counts))  # count number of features with Tn5 insertions per cell
cell.counts.z <- as.numeric(scale(cell.counts)) # convert features counts into Z-scores
cell.counts.threshold <- max(c((10^cell.counts[cell.counts.z < -3]), 100)) # minimum feature counts (greater of 1 std or 1000)


# clean sparse counts matrix ---------------------------------------------
soc.obj <- cleanData(soc.obj, 
                     min.c=cell.counts.threshold,  # minimum number of accessible features per cell
                     min.t=0.01, 		   # minimum feature frequency across cells
                     max.t=0,		     	   # maximum feature frequency across cells
                     verbose=T)


# normalize with TFIDF ---------------------------------------------------
soc.obj <- tfidf(soc.obj, doL2=T)



# project with SVD -------------------------------------------------------
soc.obj <- reduceDims(soc.obj,
                      method="SVD", 
                      n.pcs=25,
                      cor.max=0.4,
                      num.var=100000,
                      verbose=T,
                      scaleVar=T,
                      doSTD=F,
                      doL1=F,
                      doL2=T,
                      refit_residuals=F)

soc.obj <- projectUMAP(soc.obj, metric="cosine", k.near=15, svd_slotName="PCA", umap_slotName="UMAP")

#detect and remove potential doublets-----------------------------------
soc.obj <- detectDoublets(soc.obj, threads=8, nTrials=5, nSample=1000, rdMethod = "SVD", svd_slotName="PCA")
soc.obj <- filterDoublets(soc.obj, filterRatio=1.0, embedding="UMAP", libraryVar=c("library"), umap_slotname="UMAP", svd_slotname="PCA", removeDoublets=T)

# remove batch effects with harmony --------------------------------------
ids <- rownames(soc.obj$PCA)
ref.obj <- HarmonyMatrix(soc.obj$PCA, meta_data=soc.obj$meta, 
                         vars_use=c("library"), 
                         do_pca=F,
                         theta=c(2),
                         tau=c(5),
                         sigma=0.1,
                         lambda=c(0.1),
                         nclust=50,
                         max.iter.cluster=100,
                         max.iter.harmony=50,
                         return_object=T)

# process
soc.obj$h.PCA <- t(ref.obj$Z_corr) 
colnames(soc.obj$h.PCA) <- paste0("PC_", 2:(ncol(soc.obj$PCA)+1))
rownames(soc.obj$h.PCA) <- ids

# l2 normalize
soc.obj$l2.PCA <- t(apply(soc.obj$h.PCA, 1, function(x){(x/sqrt(sum(x^2)))}))

# reduce to 2-dimensions with UMAP ---------------------------------------
soc.obj <- projectUMAP(soc.obj, metric="correlation", k.near=15, svd_slotName="l2.PCA", umap_slotName="h.UMAP")

# identify clusters using neighborhood graph -----------------------------
soc.obj <- callClusters(soc.obj, 
                        res=0.5,
                        k.near=15,
                        verbose=T,
                        cleanCluster=T,
                        cl.method=3,
                        e.thresh=3,
                        m.clst=25,
						min.reads=5e5,
                        svd_slotName="l2.PCA",
                        umap_slotName="h.UMAP",
                        cluster_slotName="h.Clusters")

# plot cluster membership on UMAP embedding ------------------------------

pdf(paste0(out,".UMAP.harmony.clusters.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="h.Clusters", cex=0.5)
dev.off()

pdf(paste0(out,".UMAP.harmony.library.pdf"), width=16, height=16)
plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="library", cex=0.5, main = "Library")
dev.off()

pdf(paste0(out,".UMAP.harmony.QC.pdf"), width=12, height=3)
layout(matrix(c(1:4), nrow=1, ncol = 4))

plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="log10nSites", cex=0.2, main = "Total nSites(log10)")

plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="pOrg", cex=0.2, main = "Fraction reads in organelle")

plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="pTSS", cex=0.2, main = "Fraction reads in TSS")

plotUMAP(soc.obj, cluster_slotName="h.Clusters", column="d.type", cex=0.2, main = "Doublet distribution")

dev.off()


# save data --------------------------------------------------------------
saveRDS(soc.obj, file=paste0(out,".soc.processed.rds"))

