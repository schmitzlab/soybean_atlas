
# load libraries
library(Matrix)
library(nabor)
library(dplyr)
library(viridis)
library(data.table)
library(scales)
library(mgcv)
library(gplots)
library(RColorBrewer)
library(parallel)
library(speedglm)
library(gtools)
library(uwot)
library(splines)
library(RANN)
library(dtwclust)

# functions
loadData   <- function(sm, mt, gn, ma, rd, doL2=0){
  
  # load sparse/dense data
  message(" - loading data ...")
  sm <- readRDS(sm)
  sm@x <- rep(1, length(sm@x))
  mt <- readRDS(mt)
  gn <- readRDS(gn)
  
  # meta/svd
  ma <- read.table(ma, comment.char="")
  rd <- read.table(rd)
  ids <- intersect(rownames(ma), colnames(sm))
  ids <- intersect(ids, colnames(gn))
  ids <- intersect(ids, colnames(mt))
  
  # filter on IDs
  sm <- sm[,ids]
  gn <- gn[,ids]
  mt <- mt[,ids]
  ma <- ma[ids,]
  
  # filter on counts
  sm <- sm[Matrix::rowSums(sm)>0,]
  sm <- sm[,Matrix::colSums(sm)>0]
  gn <- gn[Matrix::rowSums(gn)>0,]
  gn <- gn[,Matrix::colSums(gn)>0]
  ids <- intersect(rownames(ma), colnames(sm))
  ids <- intersect(ids, colnames(gn))
  ids <- intersect(ids, colnames(mt))
  sm <- sm[,ids]
  gn <- gn[,ids]
  mt <- mt[,ids]
  ma <- ma[ids,]
  rd <- rd[ids,]
  
  # L2
  if(doL2==1){
    rd <- t(apply(rd, 1, function(x) x/sqrt(sum(x^2))))
  }else if(doL2==2){
    rd <- apply(rd, 2, function(x) x/sqrt(sum(x^2)))
  }
  
  # add meta data
  ma$log10ATAC <- log10(Matrix::colSums(sm))
  ma$log10RNA <- log10(Matrix::colSums(gn))
  
  # return list
  return(list(ATAC=sm, Motifs=mt, RNA=gn, meta=ma, svd=rd))
}
loadTests <- function(){
  
  # load data
  maindir <- "./pseudotime/"
  axis <- read.table(paste0(maindir,"Axis_Parenchyma/Axis_Parenchyma.Genes.diffTestsQVALs.txt"))
  coty <- read.table(paste0(maindir,"Cotyledon_Parenchyma/Cotyledon_Parenchyma.Genes.diffTestsQVALs.txt"))
  epid <- read.table(paste0(maindir,"epidermis_pseudotime/epidermis.Genes.diffTestsQVALs.txt"))
  prov <- read.table(paste0(maindir,"provascular/provascular.Genes.diffTestsQVALs.txt"))
  samr <- read.table(paste0(maindir,"SAMRAM/SAMRAM.Genes.diffTestsQVALs.txt"))
  
  # filter by fdr
  #axis <- subset(axis, axis$qval < 0.01)
  #coty <- subset(coty, coty$qval < 1e-4)
  #epid <- subset(epid, epid$qval < 0.05)
  #prov <- subset(prov, prov$qval < 0.05)
  #samr <- subset(samr, samr$qval < 0.05)
  axis <- head(axis[order(axis$pval, decreasing=F),], n=5000)
  coty <- head(coty[order(coty$pval, decreasing=F),], n=5000)
  epid <- head(epid[order(epid$pval, decreasing=F),], n=5000)
  prov <- head(prov[order(prov$pval, decreasing=F),], n=5000)
  samr <- head(samr[order(samr$pval, decreasing=F),], n=5000)
  
  # add geneID column
  axis$geneID <- rownames(axis)
  coty$geneID <- rownames(coty)
  epid$geneID <- rownames(epid)
  prov$geneID <- rownames(prov)
  samr$geneID <- rownames(samr)
  
  # change rownames
  rownames(axis) <- seq(1:nrow(axis))
  rownames(coty) <- seq(1:nrow(coty))
  rownames(epid) <- seq(1:nrow(epid))
  rownames(prov) <- seq(1:nrow(prov))
  rownames(samr) <- seq(1:nrow(samr))
  
  # append cell type info
  axis$branch <- "Axis_Parenchyma"
  coty$branch <- "Cotyledon_Parenchyma"
  epid$branch <- "Epidermis"
  prov$branch <- "Provascular"
  samr$branch <- "Apical_Meristem"
  
  # join
  all <- rbind(axis, coty, epid, prov, samr)
  
  # filter for unique
  all <- all[order(all$pval, decreasing=F),]
  all <- all[!duplicated(all$geneID),]
  
  # return
  return(all)
  
}
pseudoScores <- function(obj, pt, data.type="ATAC", threads=1, featureMin=0){
  
  # select data type:
  # ATAC
  # RNA
  # Motifs
  binary <- obj[[data.type]]
  
  # filter features
  binary <- binary[Matrix::rowSums(binary > 0)>featureMin,]
  binary <- binary[,Matrix::colSums(binary > 0)>0]
  message(" - nuclei x features: ",ncol(binary)," | ", nrow(binary))
  
  # iterate over each trajectory
  ptt <- lapply(names(pt), function(zz){
    
    # verbose
    message(" - running generalized additive model for ", zz, "...")
    z <- pt[[zz]]
    
    # pt ids
    n.ids <- rownames(z)
    bb <- binary[,n.ids]
    
    # append meta data to pt meta
    log10ATAC <- obj$meta$log10ATAC
    log10RNA <- obj$meta$log10RNA
    names(log10ATAC) <- rownames(obj$meta)
    names(log10RNA) <- rownames(obj$meta)
    z$log10ATAC <- log10ATAC[rownames(z)]
    z$log10RNA <- log10RNA[rownames(z)]
    
    # generalized additive model for logistic regression
    newdat <- data.frame(p.time=seq(from=min(z$trajectory, na.rm=T), to=max(z$trajectory, na.rm=T), length.out=500))
    fit <- mclapply(seq(1:nrow(bb)),function(x){
      
      # verbose
      if((x %% 1000)==0){message(" - iterated over ", x, " records...")}
      
      # create input data frame
      if(data.type=='ATAC'){
        df <- data.frame(acc=as.numeric(bb[x,]), p.time=z$trajectory, lib=z$log10ATAC)
        fam <- binomial()
      }else if(data.type=='RNA'){
        df <- data.frame(acc=as.numeric(bb[x,]), p.time=z$trajectory, lib=z$log10RNA)
        fam <- poisson()
      }else{
        df <- data.frame(acc=as.numeric(bb[x,]), p.time=z$trajectory, lib=z$log10ATAC)
        fam <- gaussian()
      }
      
      # residuals
      mod.scores <- glm(acc~lib, family=fam, data=df)
      df$res.acc <- residuals(mod.scores)
      
      # smooth fit
      mod <- gam(res.acc~s(p.time, bs="cr"), data=df)
      
      # get predictions 
      pred <- predict(mod, newdat, type="response")
      #zscore <- (pred - mean(pred, na.rm=T))/sd(pred, na.rm=T)
      
      # return
      return(pred)
    }, mc.cores=threads)
    
    # merge scores
    fit <- do.call(rbind, fit)
    rownames(fit) <- rownames(bb)
  
    # sort rows
    #row.o <- apply(fit, 1, which.max)
    #fit <- fit[order(row.o, decreasing=F),]
      
    # return
    return(fit)
  
  })
  names(ptt) <- names(pt)
  
  
  # return
  return(ptt)
  
}
pseudoScores2 <- function(obj, pt, data.type="ATAC", threads=1, featureMin=0){
  
  # select data type:
  # ATAC
  # RNA
  # Motifs
  binary <- obj[[data.type]]
  
  # filter features
  binary <- binary[Matrix::rowSums(binary > 0)>featureMin,]
  binary <- binary[,Matrix::colSums(binary > 0)>0]
  message(" - nuclei x features: ",ncol(binary)," | ", nrow(binary))
  
  # normalize
  if(data.type=="RNA" | data.type=="ATAC"){
    c.ids <- colnames(binary)
    binary <- binary %*% Diagonal(x=1e5/Matrix::colSums(binary))
    colnames(binary) <- c.ids
    
  }
  
  # iterate over each trajectory
  ptt <- lapply(names(pt), function(zz){
    
    # verbose
    message(" - running generalized additive model for ", zz, "...")
    z <- pt[[zz]]
    
    # pt ids
    n.ids <- rownames(z)
    bb <- binary[,n.ids]
    
    # append meta data to pt meta
    log10ATAC <- obj$meta$log10ATAC
    log10RNA <- obj$meta$log10RNA
    names(log10ATAC) <- rownames(obj$meta)
    names(log10RNA) <- rownames(obj$meta)
    z$log10ATAC <- log10ATAC[rownames(z)]
    z$log10RNA <- log10RNA[rownames(z)]
    
    # generalized additive model for logistic regression
    newdat <- data.frame(p.time=seq(from=min(z$trajectory, na.rm=T), to=max(z$trajectory, na.rm=T), length.out=500))
    fit <- mclapply(seq(1:nrow(bb)),function(x){
      
      # verbose
      if((x %% 1000)==0){message(" - iterated over ", x, " records...")}
      
      # create input data frame
      if(data.type=='ATAC'){
        df <- data.frame(acc=as.numeric(bb[x,]), p.time=z$trajectory, lib=z$log10ATAC)
        fam <- binomial()
      }else if(data.type=='RNA'){
        df <- data.frame(acc=as.numeric(bb[x,]), p.time=z$trajectory, lib=z$log10RNA)
        fam <- poisson()
      }else{
        df <- data.frame(acc=as.numeric(bb[x,]), p.time=z$trajectory, lib=z$log10ATAC)
        fam <- gaussian()
      }
      
      # residuals
      #mod.scores <- glm(acc~lib, family=fam, data=df)
      #df$res.acc <- residuals(mod.scores)
      
      # smooth fit
      mod <- gam(acc~s(p.time, bs="cr"), data=df)
      
      # get predictions 
      pred <- predict(mod, newdat, type="response")
      #zscore <- (pred - mean(pred, na.rm=T))/sd(pred, na.rm=T)
      
      # return
      return(pred)
    }, mc.cores=threads)
    
    # merge scores
    fit <- do.call(rbind, fit)
    rownames(fit) <- rownames(bb)
    
    # sort rows
    #row.o <- apply(fit, 1, which.max)
    #fit <- fit[order(row.o, decreasing=F),]
    
    # return
    return(fit)
    
  })
  names(ptt) <- names(pt)
  
  
  # return
  return(ptt)
  
}
clustPseudo <- function(ptsc){
  
  clusts <- lapply(names(ptsc), function(x){
    mat <- t(apply(ptsc[[x]],1,function(z){(z-mean(z, na.rm=T))/sd(z, na.rm=T)}))
    mat <- t(apply(mat, 1, rescale, c(-1,1)))
    km <- kmeans(mat, centers=35)
    mat <- mat[order(km$cluster, decreasing=F),]
    
    
  })
  names(clusts) <- names(ptsc)
  return(clusts)
  
}
plotSC <- function(fit, sites=NULL, prefix="ATAC", col="ATAC", nonhclust=F,
                   optimalOrder=F, k.means=F, k=5, filt.thresh=F){  
  
  # if optimal order
  if(optimalOrder){
    nonhclust <- T
    outs <- lapply(names(fit), function(z){
      zz <- fit[[z]]
      
      # z-score
      zz <- t(apply(zz, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
      
      # plot top 50% by variance
      zz <- t(apply(zz, 1, function(x){rescale(x, c(-1,1))}))
      
      row.o <- t(apply(zz, 1, rank))
      row.o 
    })
    outs <- Reduce(`+`, outs) / length(outs)
    ave <- apply(outs, 1, which.max)
    
  }
  
  # join
  fit <- do.call(cbind, fit)
  
  # optimal order
  if(optimalOrder){
    fit <- fit[order(ave, decreasing=F),]
  }
  
  # sites
  if(!is.null(sites)){
    fit <- fit[sites,]
  }
  
  # z-score
  fit <- t(apply(fit, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
  
  if(filt.thresh){
    fit <- fit[rowSums(fit > filt.thresh)>0,]
    message(" - ",nrow(fit), " features remaining after z-score threshold ...")
  }
  
  # plot top 50% by variance
  fit <- t(apply(fit, 1, function(x){rescale(x, c(-1,1))}))
  
  # kmeans
  if(k.means){
    nonhclust <- T
    cluster <- kmeans(fit, centers=k)
    fit <- fit[order(cluster$cluster, decreasing=F),]
    
  }
  
  # plot
  message(" - plotting cell trajectory ...")
  if(col=="ATAC"){
    cols <- colorRampPalette(c("paleturquoise4", "white","palevioletred3"))(100)
  }else if(col=="RNA"){
    cols <- colorRampPalette(c("grey80", "grey75",brewer.pal(7, "YlGnBu")[2:7]))(100)  
  }else{
    cols <- colorRampPalette(c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3"))(100)
  }
  
  pdf(paste0(prefix,".trajectoryHM.pdf"), width=10, height=10)
  heatmap.2(fit, trace="none", col=cols, Colv=NA, Rowv=ifelse(nonhclust==T, NA, T), 
            dendrogram="none",
            scale="none", labCol=NA, useRaster=T, 
            ylab=paste(col, paste0("(n=",nrow(fit),")"), sep=" "))
  dev.off()
  
}
plotIndSC <- function(fit, prefix=NULL, col=NULL){
  
    # order rows
    outs <- lapply(names(fit), function(z){
      zz <- fit[[z]]
      
      # z-score
      zz <- t(apply(zz, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
      
      # plot top 50% by variance
      zz <- t(apply(zz, 1, function(x){rescale(x, c(-1,1))}))
      
      row.o <- apply(zz, 1, which.max)
      zz[order(row.o, decreasing=F),]
    })
    names(outs) <- gsub("_pseudotime\\.trajectory\\.txt","",names(fit))
    
    # plot cols
    message(" - plotting cell trajectory ...")
    if(col=="ATAC"){
      cols <- colorRampPalette(c("paleturquoise4", "white","palevioletred3"))(100)
    }else if(col=="RNA"){
      cols <- colorRampPalette(c("grey80", "grey75",brewer.pal(7, "YlGnBu")[2:7]))(100)  
    }else{
      cols <- colorRampPalette(c("dodgerblue4", "deepskyblue","grey85", "darkorange", "firebrick3"))(100)
    }
  
    # plot
    for(i in names(outs)){
      ftt <- outs[[i]]
      pdf(paste0(prefix,".",i,".trajectoryHM.pdf"), width=10, height=10)
      heatmap.2(ftt, trace="none", col=cols, Colv=NA, Rowv=F, 
                dendrogram="none",
                scale="none", labCol=NA, useRaster=T, 
                ylab=paste(col, paste0("(n=",nrow(ftt),")"), sep=" "))
      dev.off()
      
    }
}

# load data
pt.data <- list.files(pattern="*.trajectory.txt")
pt.data <- lapply(pt.data, function(x){
  df <- read.table(x, comment.char="")
  df <- df[!is.na(df$trajectory),]
  df <- df[order(df$trajectory, decreasing=F),]
  return(df)
})
names(pt.data) <- list.files(pattern="*.trajectory.txt")

# sparse data
sparse <- "data/Gm_atlas_ACRs_imputed_on_RNA.rds"
motif <- "data/Gm_atlas_Motif_deviation_imputed_on_RNA.rds"
genes <- "RNA_sparse_gene"
meta <- "RNA_embryo_cells.metadata.inferred_age.PT.txt"
svd <- "soybean_embryo_snRNA.harmonyPCs.txt"

# load sparse data
obj <- loadData(sparse, motif, genes, meta, svd)

# tests
tests <- loadTests()

# get pseudotime scores
ptsc.motif <- pseudoScores2(obj, pt.data, data.type="Motifs", featureMin=0, threads=30)
ptsc.rna <- pseudoScores2(obj, pt.data, data.type="RNA", featureMin=10, threads=30)

# plot
plotSC(ptsc.motif, prefix="Motifs", col="Motifs")
plotSC(ptsc.rna, sites=tests, prefix="RNA", col="RNA", filt.thresh=F)

# plot individual trajectories
plotIndSC(ptsc.rna, prefix="Motif", col="Motif")
plotIndSC(ptsc.rna, prefix="RNA", col="RNA")

# plot predefined order
tests$branch <- factor(tests$branch, levels=c("Axis_Parenchyma","Cotyledon_Parenchyma",
                                              "Epidermis","Provascular","Apical_Meristem"))
tests <- tests[order(tests$branch,decreasing=F),]
its <- 0
row.ids <- lapply(names(ptsc.rna), function(z){
  its <<- its + 1
  zz <- ptsc.rna[[z]]
  shared <- intersect(tests$geneID, rownames(zz))
  zz <- zz[shared,]
  tests <- tests[tests$geneID %in% shared,]
  
  # cluster sub
  branch <- subset(tests, as.numeric(tests$branch)==its)
  s.geneIDs <- branch$geneID
  sub.zz <- zz[s.geneIDs,]
  row.o <- apply(sub.zz, 1, which.max)
  rownames(sub.zz)[order(row.o, decreasing=F)]
  
})
row.ids <- unlist(row.ids)
plotSC(ptsc.rna, sites=row.ids, prefix="RNA", col="RNA", nonhclust=T)




















# split scaled RNA
axis <- mptsc.rna[,1:500]
coty <- mptsc.rna[,501:1000]
epid <- mptsc.rna[,1001:1500]
prov <- mptsc.rna[,1501:2000]
meri <- mptsc.rna[,2001:2500]

# print top genes
for(i in unique(fgenes$celltype)){
  print(head(fgenes[fgenes$celltype==i,], n=25))
}


# markers
its <- 0
topgenes <- lapply(names(ptsc.rna), function(z){
  its <<- its + 1
  zz <- ptsc.rna[[z]]
  shared <- intersect(tests$geneID, rownames(zz))
  zz <- zz[shared,]
  tests2 <- tests[tests$geneID %in% shared,]
  
  # cluster sub
  branch <- subset(tests2, as.numeric(tests2$branch)==its)
  s.geneIDs <- branch$geneID
  sub.zz <- zz[s.geneIDs,]
  vals <- apply(sub.zz, 1, which.max)
  score <- apply(sub.zz, 1, max)
  val <- data.frame(geneID=names(vals), branchHit=vals, score=score, celltype=z)
  
  return(val)
})
topgenes <- do.call(rbind, topgenes)

# filter
fgenes <- subset(topgenes, topgenes$branchHit >=250)
zscore <- apply(mptsc.rna, 1, max)
fgenes$zscore <- zscore[rownames(fgenes)]
fgenes <- fgenes[order(fgenes$celltype, -fgenes$zscore),]
fgenes$AtName <- Atname[fgenes$geneID]
write.table(fgenes, file="topHits_embryogenesis_branches.ranked.txt", quote=F, row.names=T, col.names=T, sep="\t")

# plot ATML1 Glyma.10G251300
geneID <- "Glyma.17G256500"
geneName <- "MP"
atml1.axis <- as.numeric(axis[geneID,])
atml1.coty <- as.numeric(coty[geneID,])
atml1.epid <- as.numeric(epid[geneID,])
atml1.prov <- as.numeric(prov[geneID,])
atml1.meri <- as.numeric(meri[geneID,])
yran <- range(atml1.axis, atml1.coty, atml1.epid, atml1.prov, atml1.meri)
cols <- c("#529BA3", "#417691", "#FBB363", "#CDCFD6", "#BC81B8")

pdf(paste0(geneName,"-", geneID, ".pseudotimeProfile.pdf"), width=5, height=5)
plot(atml1.axis, type="l", col=cols[1], ylim=yran, ylab="Scaled gene expression", xlab="Pseudotime", main=geneName, lwd=2)
lines(atml1.coty, col=cols[2], lwd=2)
lines(atml1.epid, col=cols[3], lwd=2)
lines(atml1.prov, col=cols[4], lwd=2)
lines(atml1.meri, col=cols[5], lwd=2)
grid(lty=1)
dev.off()



# plot motif
mmotif <- do.call(cbind, motif)
mmotif <- t(apply(mmotif, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))



