## infer developmental time

# load libraries
library(glmnet)
library(Matrix)
library(caret)
library(RColorBrewer)
library(mgcv)
library(scales)
library(gplots)
library(vioplot)
library(parallel)

# functions
smooth.data<- function(x, k=15, step=3, npcs=30, df=NULL, rds=NULL, verbose=F){
  
  # verbose
  if(verbose){message(" - imputing gene activity ...")}
  
  # input
  data.use <- x
  
  # verbose
  if(!is.null(rds)){
    
    if(!is.null(df)){
      if(verbose){message("   * using UMAP manifold for smoothing ...")}
      pcs <- df[,c("umap1","umap2")]
    }else{
      if(verbose){message("   * using prior PC space as manifold ...")}
      pcs <- rds[colnames(x),c(1:npcs)]
    }
  }else{
    
    # LSI
    if(verbose){message("   * PC manifold set to NULL, running LSI (TFIDF)...")}
    x[x>0] <- 1
    tf.idf <- tfidf(x)
    
    # get PCS
    if(verbose){message("   * PC manifold set to NULL, running LSI ...")}
    pc <- irlba(t(tf.idf), npcs)
    pcs <- pc$u 
    rownames(pcs) <- colnames(x)
    colnames(pcs) <- paste0("PC_", seq(1:ncol(pcs)))
    
    # do l2-norm
    pcs <- apply(pcs, 2, function(x){x/sqrt(sum(x^2))})
  }
  
  # get KNN
  if(verbose){message("   * finding knn graph ...")}
  knn.graph <- nn2(pcs, k=k, eps=0)$nn.idx
  j <- as.numeric(x = t(x = knn.graph))
  i <- ((1:length(x = j)) - 1) %/% k + 1
  edgeList = data.frame(i, j, 1)
  A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])
  
  # Smooth graph
  if(verbose){message("   * smoothing graph ...")}
  A = A + t(A)
  A = A / Matrix::rowSums(A)
  step.size = step
  if(step.size > 1){
    for(i in 1:step.size){
      if(verbose){message("     ~ step ",i)}
      A = A %*% A
    }
  }
  
  # smooth data
  if(verbose){message("   * smoothing activity ...")}
  impute.activity <- t(A %*% t(data.use))
  colnames(impute.activity) <- colnames(x)
  rownames(impute.activity) <- rownames(x)
  
  # clean empty rows (if present) and round to two decimal places
  if(verbose){message("   * remove empty rows/columns and scale to per million ...")}
  impute.activity <- impute.activity[Matrix::rowSums(impute.activity)>0,]
  impute.activity <- impute.activity[,Matrix::colSums(impute.activity)>0]
  impute.activity <- impute.activity %*% Diagonal(x=1e6/Matrix::colSums(impute.activity))
  impute.activity@x <- round(impute.activity@x, digits=2)
  
  # return sparse Matrix
  return(impute.activity)
}
plotTrajGN <- function(obj, pt, prefix="temp", threads=1){
  
  # subset traj
  message(" - Sort cells ...")
  pt <- pt[order(pt$inferred_age, decreasing=F),]
  binary <- obj[,rownames(pt)]
  binary <- binary[,Matrix::colSums(binary)>0]
  binary <- binary[Matrix::rowSums(binary)>0,]
  
  # generalized additive model for logistic regression
  message(" - running generalized additive model ...")
  newdat <- data.frame(p.time=seq(from=min(pt$inferred_age, na.rm=T), to=max(pt$inferred_age, na.rm=T), length.out=500))
  fit <- mclapply(seq(1:nrow(binary)),function(x){
    df <- data.frame(acc=as.numeric(binary[x,]), p.time=pt$inferred_age)
    mod <- gam(acc~s(p.time, bs="cr"), data=df)
    pred <- predict(mod, newdat, type="response")
    zscore <- (pred-mean(pred, na.rm=T))/sd(pred, na.rm=T)
  }, mc.cores=threads)
  names(fit) <- rownames(binary)
  fit <- do.call(rbind, fit)
  write.table(fit, file=paste0(prefix,".Gene_pt.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  
  # filter by variances
  fit <- t(apply(fit, 1, function(x){
    rescale(x, c(0,1))
  }))
  
  # reformat output
  row.o <- apply(fit, 1, which.max)
  fit <- fit[order(row.o, decreasing=F),]
  
  # plot
  message(" - plotting cell trajectory ...")
  cols <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
  pdf(paste0(prefix,".trajectoryGene.pdf"), width=10, height=10)
  heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=NA, dendrogram="none",
            scale="none", labRow = NA, labCol=NA, useRaster=T,
            ylab=paste("Genes", paste0("(n=",nrow(fit),")"), sep=" "))
  dev.off()
  
  # return
  return(fit)
  
}

# load data
atacm <- read.table("ATAC.metadata.annotated.txt", comment.char="")
rnam <- read.table("soybean_embryo_snRNA.metaData.txt", comment.char="")
atacs <- readRDS("Gm_atlas_seed_4_stages.embryo.scATAC_ACR_sparse.rds")
rnas <- readRDS("Gm_atlas_seed_4_stages.embryo.snRNA_gene_sparse.rds")

# normalize sparse matrices
c.ids <- colnames(rnas)
rnas@x <- log1p(rnas@x)
rnas <- rnas %*% Diagonal(x=1e5/Matrix::colSums(rnas))
colnames(rnas) <- c.ids

# remove seed coat
atacm <- atacm[atacm$emb_celltype != 'Seed_coat',]
rnam <- rnam[rnam$emb_celltype != 'Seed_coat',]

# subsample cells by stage
min.atac <- min(table(atacm$TissueType))
min.rna <- min(table(rnam$TissueType))

# select equal numbers of cells
atac.sub <- lapply(unique(atacm$TissueType), function(z){
  message(z)
  df <- subset(atacm, atacm$TissueType==z)
  df[sample(nrow(df), min.atac),]
})
atac.sub <- do.call(rbind, atac.sub)
rna.sub <- lapply(unique(rnam$TissueType), function(z){
  message(z)
  df <- subset(rnam, rnam$TissueType==z)
  df[sample(nrow(df), min.rna),]
})
rna.sub <- do.call(rbind, rna.sub)
atacss <- atacs[,rownames(atac.sub)]
rnass <- rnas[,rownames(rna.sub)]
atacss <- atacss[Matrix::rowSums(atacss)>0,]
rnass <- rnass[Matrix::rowSums(rnass)>0,]
colnames(rnass) <- rownames(rna.sub)

# set dev time
atac.sub$devtime <- as.numeric(factor(atac.sub$TissueType, levels=c("Globular", "Heart", "Cotyledon", "Early_Maturation")))
rna.sub$devtime <- as.numeric(factor(rna.sub$TissueType, levels=c("Globular", "Heart", "Cotyledon", "Early_Maturation")))

# split
atac.trainidx <- createDataPartition(atac.sub$TissueType, p=10/11, list=F, times=1)
rna.trainidx <- createDataPartition(rna.sub$TissueType, p=10/11, list=F, times=1)
atac.train <- atac.sub[atac.trainidx,]
atac.test <- atac.sub[-atac.trainidx,]
rna.train <- rna.sub[rna.trainidx,]
rna.test <- rna.sub[-rna.trainidx,]

# train
atac.model <- cv.glmnet(t(atacss[,rownames(atac.train)]), atac.train$devtime)
rna.model <- cv.glmnet(t(rnass[,rownames(rna.train)]), rna.train$devtime)

# test
atac.pred <- as.data.frame(predict(atac.model, newx=t(atacss[,rownames(atac.test)]), s='lambda.min', type="response"))
rna.pred <- as.data.frame(predict(rna.model, newx=t(rnass[,rownames(rna.test)]), s='lambda.min', type="response"))
atac.pred$obs <- atac.test$devtime
rna.pred$obs <- rna.test$devtime

# plot test boxplot
pdf("vioplot_obs_inferred.RNA.test.pdf", width=5, height=5)
vioplot(rna.pred$lambda.min~rna.pred$obs, xlab="Predicted developmental time", ylab="Observed developmental time", main="RNA")
dev.off()

# plot test
pdf("scatter_obs_inferred.RNA.test.pdf", width=5, height=5)
plot(rna.pred$lambda.min, rna.pred$obs, pch=16, xlab="Predicted developmental time", ylab="Observed developmental time", main="RNA")
dev.off()

# collect dev time for all cells
atac.age <- as.data.frame(predict(atac.model, newx=t(atacs[rownames(atacss),]), s='lambda.min', type="response"))
rna.age <- as.data.frame(predict(rna.model, newx=t(rnas[rownames(rnass),]), s='lambda.min', type="response"))
atacm$inferred_age <- atac.age[rownames(atacm),1]
rnam$inferred_age <- rna.age[rownames(rnam),1]
atacm$devtime <- as.numeric(factor(atacm$TissueType, levels=c("Globular", "Heart", "Cotyledon", "Early_Maturation")))
rnam$devtime <- as.numeric(factor(rnam$TissueType, levels=c("Globular", "Heart", "Cotyledon", "Early_Maturation")))

# plot all predictions
pdf("scatter_obs_inferred.RNA.ALL.pdf", width=5, height=5)
plot(rnam$devtime, rnam$inferred_age, pch=16, xlab="Obs dev age", ylab="Inferred dev age", main="RNA")
dev.off()

pdf("inferred_vs_observed_developmental_time.pdf", width=5, height=5)
boxplot(rnam$inferred_age~factor(rnam$devtime), outline=F)
dev.off()

# coefficients ----------------------------------
rna.coef <- coefficients(rna.model, s='lambda.min')
rna.coe <- rna.coef[-1,1]
rna.coe <- rna.coe[order(rna.coe, decreasing=F)]
rna.coe <- rna.coe[rna.coe != 0]

pdf("age_coefficients.pdf", width=5, height=5)
cols <- colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)
plot(seq(1:length(rna.coe)), rna.coe, pch=16, cex=0.5, col=cols[cut(rna.coe, breaks=seq(from=-0.015, to=0.015, length.out=101))])
grid(lty=1)
dev.off()

# plot dev age gene programs
plotTrajGN(rnas, rnam, prefix="Inferred_AGE", threads=10)


# plot UMAP -------------------------------------

# parameters
pdf("Inferred_developmental_age.pdf", width=16, height=8)
layout(matrix(c(1:2), nrow=1))

# ATAC
cols <- colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)
cols <- cols[cut(atacm$inferred_age, breaks=101)]
plot(atacm$umap1, atacm$umap2, pch=16, cex=0.5, col=cols, main="ATAC")

# RNA
cols <- colorRampPalette(c("grey70","grey75", brewer.pal(9, "RdPu")[3:9]), bias=0.7)(100)
cols <- cols[cut(rnam$inferred_age, breaks=101)]
plot(rnam$umap1, rnam$umap2, pch=16, cex=0.3, col=cols, main="RNA")

# device off
dev.off()


