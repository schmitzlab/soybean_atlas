###################################################################################################
## project motif scores on umap embeddings
###################################################################################################

# setup arguments
args <- commandArgs(T)
if(length(args) != 4){
	stop("Rscript project.motif.UMAP.R <reducedDims> <motif.scores> <pcs> <prefix>")
}

# variables
input <- as.character(args[1])
motifs <- as.character(args[2])
pc.dat <- as.character(args[3])
prefix <- as.character(args[4])

# load libs
library(RColorBrewer)
library(viridis)
library(scales)
library(pheatmap)
library(gplots)
library(RANN)
library(Matrix)
library(parallel)
library(png)

smoothByLouvain <- function(x, y, pc, threads=10){
    
    # transform to + values
    x <- t(apply(x, 1, function(z){
        z[is.na(z)] <- 0
	z <- z - min(z)
        return(z)
    }))
    
    # iterate over cluster
    clusts <- unique(y$LouvainClusters)
    motifs.ids <- rownames(x)
    out <- smooth.data(x, rds=pc, k=15, step=3, npcs=ncol(pc), threads=threads)
    out <- as.matrix(t(scale(t(out))))
    return(out)
}
smooth.data        <- function(x, 
                               k=25, 
                               step=3, 
                               npcs=30, 
                               df=NULL, 
                               rds=NULL,
                               n.perms=10,
                               threads=1){
    
    # verbose
    message(" - imputing gene activity ...")
    
    # input
    data.use <- x
    
    # hidden functions
    .markov_affinity <- function(com, dat.use, step=3, k=15){
        
        # get KNN
        knn.graph <- nn2(com, k=k, eps=0)$nn.idx
        j <- as.numeric(x = t(x = knn.graph))
        i <- ((1:length(x = j)) - 1) %/% k + 1
        edgeList = data.frame(i, j, 1)
        A = sparseMatrix(i = edgeList[,1], j = edgeList[,2], x = edgeList[,3])
        
        # Smooth graph
        A = A + t(A)
        A = A / Matrix::rowSums(A)
        step.size = step
        if(step.size > 1){
            for(i in 1:step.size){
                A = A %*% A
            }
        }
        
        # smooth data
        im.activity <- t(A %*% t(dat.use))
        colnames(im.activity) <- colnames(dat.use)
        rownames(im.activity) <- rownames(dat.use)
        
        # return sparse Matrix
        return(im.activity)
        
    }
    
    # verbose
    if(!is.null(rds)){
        
        if(!is.null(df)){
            message("   * using UMAP manifold for smoothing ...")
            pcs <- df[,c("umap1","umap2")]
        }else{
            message("   * using prior PC space as manifold ...")
            if(npcs > ncol(rds)){
                npcs <- ncol(rds)
            }
            pcs <- rds[colnames(x),c(1:npcs)]
        }
    }else{
        
        # LSI
        message("   * PC manifold set to NULL, running LSI (TFIDF)...")
        x[x>0] <- 1
        tf.idf <- tfidf(x)
        
        # get PCS
        message("   * PC manifold set to NULL, running LSI ...")
        pc <- irlba(t(tf.idf), npcs)
        pcs <- pc$u
        rownames(pcs) <- colnames(x)
        colnames(pcs) <- paste0("PC_", seq(1:ncol(pcs)))
        
        # do l2-norm
        pcs <- apply(pcs, 2, function(x){x/sqrt(sum(x^2))})
    }
    
    # check number of cells
    num.cells <- nrow(pcs)
    if(num.cells > 50000){
        message(" - number of cells > 50K, initiating multiple runs ...")
        max.val <- 10000
        x.vec <- seq(1:nrow(pcs))
        perms <- mclapply(seq(1:n.perms), function(z){
            message(" - initialized ", z, " runs ...")
            x.vec.r <- x.vec[sample(length(x.vec))]
            pcs.s <- split(pcs, ceiling(x.vec.r/max.val))
            outs <- lapply(pcs.s, function(y){
                ids <- rownames(y)
                dat <- data.use[,ids]
                .markov_affinity(y, dat, step=step, k=k)
            })
            outs <- do.call(cbind, outs)
            outs <- outs[,colnames(data.use)]
            return(outs)
        }, mc.cores=threads)
        message(" - averaging runs ...")
        imputed.activity <- Reduce("+", perms) / length(perms)
        return(imputed.activity)
        
    }else{
        imputed.activity <- .markov_affinity(pcs, data.use, step=step, k=k)
        return(imputed.activity)
    }
}

# load data
message("loading data ...")
a <- read.table(input, comment.char="")
rownames(a) <- a$cellID
bb <- t(as.matrix(read.table(motifs)))
pcs <- read.table(pc.dat)
ids <- intersect(rownames(a), colnames(bb))
non.ids <- rownames(a)[!rownames(a) %in% ids]
mat0 <- matrix(NA, nrow=nrow(bb), ncol=length(non.ids))
rownames(mat0) <- rownames(bb)
colnames(mat0) <- non.ids
bb <- cbind(bb, mat0)
ids <- intersect(rownames(a), colnames(bb))
ids <- intersect(ids, rownames(pcs))
a <- a[ids,]
bb <- bb[,ids]
pcs <- pcs[ids,]

# smooth scores
bb[is.na(bb)] <- 0
b <- smoothByLouvain(bb, a, pcs)
write.table(t(b), file=paste0(prefix,".imputed_motifs.txt"), quote=F, row.names=T, col.names=T, sep="\t")

# project motifs
num_motifs <- 200
message("collecting top ",num_motifs, " motifs ...")
motif.scores <- apply(b, 1, function(x){var(x, na.rm=T)})
motif.scores <- motif.scores[order(motif.scores, decreasing=T)]
df <- as.data.frame(motif.scores)
write.table(df, file=paste0(prefix,".imputed.variance.txt"), quote=F, row.names=T, col.names=T, sep="\t")
topM <- head(motif.scores, n=num_motifs)
motif.sub.scores <- b[rownames(b) %in% names(topM),]
motif.sub.scores <- motif.sub.scores[,rownames(a)]

# plot
message("plotting motif scores on umap embeddings ...")
cols <- 5
nums <- ceiling(num_motifs/5)
tnums <- nums*cols
png(paste0(prefix,".motif.projections.imputed.UMAP.png"), width=(5*cols), height=(5*nums), units="in", res=150, type="cairo")
layout(matrix(c(1:tnums), nrow=nums, ncol=5, byrow=T))
par(mar=c(1,1,0.25,0.25))
for(i in 1:num_motifs){
    cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.75)(100)
    scores <- as.numeric(motif.sub.scores[i,])
    quants <- quantile(scores, c(0.02, 0.98))
    scores[scores > quants[2]] <- quants[2]
    scores[scores < quants[1]] <- quants[1]
    min.q <- quants[1] - abs(quants[1]*0.1)
    max.q <- quants[2] + abs(quants[2]*0.1)
    col.vec <- cols[cut(scores, breaks=seq(from=min.q, to=max.q, length.out=101))]
    plot(a$umap1, a$umap2, pch=16, col=col.vec, cex=1.0, 
         xlab="", ylab="", xaxt="n", yaxt="n", main=rownames(motif.sub.scores)[i],
         bty="none")
}
dev.off()

