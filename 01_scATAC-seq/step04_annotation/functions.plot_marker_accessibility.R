###################################################################################################
###################################################################################################
##                                                                                               ##
##                  functions for plotting accessibility of markers from cicero                  ##
##                                                                                               ##
###################################################################################################
###################################################################################################

# load libraries
library(viridis)
library(mclust)
library(irlba)
library(Matrix)
library(RANN)
library(reshape2)
library(gtools)
library(RColorBrewer)
library(gplots)
library(scales)
library(varistran)
library(edgeR)
library(parallel)
library(png)
      

###################################################################################################
###################################################################################################
###################################################################################################

# load data
loadData           <- function(meta, 
                               pcs, 
                               geneact, 
                               markers){

    # verbose
    message(" - loading data ...")

    # load
    b <- read.table(meta)
    b <- validateCluster(b)
    marker.info <- read.table(markers, header=T)
    marker.info <- marker.info[!duplicated(marker.info$geneID),]
    rownames(marker.info) <- marker.info$geneID
    h.pcs <- read.table(pcs)
    if(grepl(".gz$", geneact)){
	message(" - gene accessibility is gzipped, opening with gzfile")
	activity <- read.table(gzfile(geneact))
    }else if(grepl(".rds", geneact)){
	message(" - gene accessibility is RDS compressed")
	activity <- readRDS(geneact)
    }else{
	message(" - gene accessibility is non-compressed, reading as is")
        activity <- read.table(geneact)
    }
    b <- b[rownames(b) %in% rownames(h.pcs),]
    h.pcs <- h.pcs[rownames(h.pcs) %in% rownames(b),]

    # reformat sparse
    if(!grepl(".rds", geneact)){
	    activity <- sparseMatrix(i=as.numeric(activity$V1),
        	                     j=as.numeric(activity$V2),
                	             x=as.numeric(activity$V3),
                        	     dimnames=list(levels(activity$V1), levels(activity$V2)))
    }
    print(head(activity[,1:5]))
    activity <- activity[,colnames(activity) %in% rownames(b)]
    activity <- activity[,Matrix::colSums(activity>0) > 0]
    activity <- activity[Matrix::rowSums(activity>0) > 0,]
    activity <- activity[,Matrix::colSums(activity>0) > 0]
    b <- b[colnames(activity),]
    h.pcs <- h.pcs[colnames(activity),]

    # sanity check
    sanity <- 0
    if(sanity > 0){
	print("SANITY CHECK")
	print(head(b))
        print(head(marker.info))
        print(head(h.pcs[,1:5]))
        print(head(activity[,1:5]))
    }

    # output
    return(list(h.pcs=h.pcs, b=b, activity=activity, marker.info=marker.info))

}

# validate cluster ID
validateCluster    <- function(df){
    if(! "Cluster" %in% colnames(df)){
        if("leiden_clusters" %in% colnames(df)){
            df$Cluster <- paste0("cluster", as.numeric(as.factor(df$LouvainClusters)))
            return(df)
        }else if("iNMF_clusters" %in% colnames(df)){
            df$Cluster <- paste0("cluster",df$iNMF_clusters)
            return(df)
        }else if("LouvainClusters" %in% colnames(df)){
	    df$Cluster <- paste0("cluster",df$LouvainClusters)
	    return(df)
	}
    }else{
        if(length(grepl("cluster",df$Cluster)==T)==0){
            df$Cluster<- paste0("cluster",df$Cluster)
            return(df)
        }else{
            return(df)
        }
    }
}



###################################################################################################
###################################################################################################
###################################################################################################

# tfidf run
safe_tfidf         <- function(tf, 
                               idf,  
                               block_size=2000e6){
    result = tryCatch({
        result = tf * idf
        result
    }, error = function(e) {
        options(DelayedArray.block.size=block_size)
        DelayedArray:::set_verbose_block_processing(TRUE)

        tf = DelayedArray(tf)
        idf = as.matrix(idf)

        result = tf * idf
        result
    })
    return(result)
}

# do tfidf
tfidf              <- function(bmat, 
                               frequencies=F, 
                               log_scale_tf=T, 
                               scale_factor=100000){

    # Use either raw counts or divide by total counts in each cell
    if (frequencies) {
        # "term frequency" method
        tf = t(t(bmat) / Matrix::colSums(bmat))
    } else {
        # "raw count" method
        tf = bmat
    }

    # Either TF method can optionally be log scaled
    if (log_scale_tf) {
        if (frequencies) {
            tf@x = log1p(tf@x * scale_factor)
        } else {
            tf@x = log1p(tf@x * 1)
        }
    }

    # IDF w/ "inverse document frequency smooth" method
    idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))

    # TF-IDF
    tf_idf_counts = safe_tfidf(tf, idf)
    rownames(tf_idf_counts) = rownames(bmat)
    colnames(tf_idf_counts) = colnames(bmat)
    return(Matrix(tf_idf_counts, sparse=T))
}



###################################################################################################
###################################################################################################
###################################################################################################

# processing
sparse_apply       <- function(Sp_X, 
                               MARGIN, 
                               FUN, 
                               convert_to_dense, ...){
    if (convert_to_dense){
        if (MARGIN == 1){
            Sp_X <- Matrix::t(Sp_X)
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(as.matrix(Sp_X[,i]), ...)
            }, FUN, ...)
        }else{
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(as.matrix(Sp_X[,i]), ...)
            }, FUN, ...)
        }
    }else{
        if (MARGIN == 1){
            Sp_X <- Matrix::t(Sp_X)
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(Sp_X[,i], ...)
            }, FUN, ...)
        }else{
            res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
                FUN(Sp_X[,i], ...)
            }, FUN, ...)
        }
    }

    if(MARGIN == 1){
        new <- Matrix(t(as.matrix(do.call(cbind, res))), sparse=T)
        return(new)
    }else{
        new <- Matrix(as.matrix(do.call(cbind, res)), sparse=T)
        return(new)
    }

}

# normalize by size factors
normalizeSF        <- function(adj.act, 
                               verbose=F){

    # verbose
    if(verbose==T){
        message(" - unnormalized activities:")
        print(head(adj.act[,1:5]))
    }

    # get per sample norm factors
    norm.factor <- estimate_sf_sparse(adj.act)
    if(verbose==T){
        message(" - normalization factors:")
        print(head(norm.factor))
    }

    # normalized counts
    norm.act <- adj.act %*% Diagonal(x=1/norm.factor)
    colnames(norm.act) <- colnames(adj.act)
    rownames(norm.act) <- rownames(adj.act)
    if(verbose==T){
        message(" - normalized activities:")
        print(head(norm.act[,1:5]))
    }

    # output
    return(list(norm.act=norm.act, norm.factor=norm.factor))
}

# estimate sf
estimate_sf_sparse <- function(counts, 
                               round_exprs=T, 
                               method="mean-geometric-mean-total"){
    if (round_exprs)
        counts <- round(counts)

    if(method == 'mean-geometric-mean-total') {
        cell_total <- Matrix::colSums(counts)
        sfs <- cell_total / exp(mean(log(cell_total)))
    }else if(method == 'mean-geometric-mean-log-total') {
        cell_total <- Matrix::colSums(counts)
        sfs <- log(cell_total) / exp(mean(log(log(cell_total))))
    }

    sfs[is.na(sfs)] <- 1
    sfs
}

# estimate shannon entropy
shannon.entropy    <- function(p){
    if (min(p) < 0 || sum(p) <= 0) return(NA)
    p.norm <- p[p>0]/sum(p)
    -sum(log2(p.norm)*p.norm)
}

# estimate gene scores
gene.score         <- function(data=data, 
                               markers=markers, 
                               rd=rd){

    # reformat data into sparseMatrix
    if(is.data.frame(data)){
        message("reformating data into sparseMatrix ...")
        data <- sparseMatrix(i=as.numeric(data$V1),
                             j=as.numeric(data$V2),
                             x=as.numeric(data$V3),
                             dimnames=list(levels(data$V1), levels(data$V2)))
    }

    # reorder matrices
    data <- data[,rownames(rd)]
    cellave <- Matrix::colSums(data)

    # collect unique genes
    genes <- unique(markers$V4)
    mat <- matrix(ncol=length(genes), nrow=nrow(rd), dimnames=list(rownames(rd), genes))
    cnts <- 0

    # iterate over unique genes
    for(i in genes){
        cnts <- cnts+1
        message("estimating accessibility for gene: ",i, " ...")
        peak.sub <- subset(markers, markers$V4==i)
        gene.peak <- as.character(droplevels(peak.sub$V8))
        if(length(gene.peak) > 1){
            gene.scores <- data[rownames(data) %in% gene.peak,]
            summed.score <- log1p((Matrix::colSums(gene.scores)/cellave)*10000)
            summed.score <- summed.score - mean(summed.score, na.rm=T)
            summed.score <- rescale(summed.score, c(0,1))
        }else{
            summed.score <- log1p((data[rownames(data) %in% gene.peak,]/cellave)*10000)
            summed.score <- summed.score - mean(summed.score, na.rm=T)
            summed.score <- rescale(summed.score, c(0,1))
        }
        ordered.score <- summed.score
        mat[,cnts] <- ordered.score
    }

    # return new reduced dimension matrix
    df <- cbind(rd, mat)
    return(df)
}

# normalize activity
normalize.activity <- function(df, 
                               acts, 
                               n.random=NULL, 
                               output="output",
                               norm.clust.means=T,
                               logTransform=F, 
                               scaleP=F,
                               plotHeatmapRaw=F){

    # verbose
    message(" - normalizing libraries ...")

    # quick clean
    acts <- acts[Matrix::rowSums(acts>0)>0,]
    acts <- acts[,Matrix::colSums(acts>0)>0]

    # check data frames/matrices
    df <- df[rownames(df) %in% colnames(acts),]
    acts <- acts[,colnames(acts) %in% rownames(df)]
    df <- df[colnames(acts),]
    acts.o <- acts
	message(" - check acts.o 1 ...")
	print(head(colnames(acts.o)))

    # if select random
    if(!is.null(n.random)){
        acts <- acts[sample(nrow(acts), n.random),]
    }

    # find cluster means
    its <- 0
    message(" - extracting average gene activity across cells per cluster ...")
    clust.means <- lapply(unique(df$Cluster), function(x){
	    clust.ids <- rownames(subset(df, df$Cluster==x))
	    Matrix::rowMeans(acts[,clust.ids])
    })
    clust.means <- do.call(cbind, clust.means)
    colnames(clust.means) <- unique(df$Cluster)

    # plot raw cluster means
    clust.means <- as.matrix(clust.means)
    c.means <- clust.means[,mixedorder(colnames(clust.means))]
    row.o <- apply(c.means, 1, which.max)
    c.means <- c.means[order(row.o, decreasing=F),]
    row.o <- rownames(c.means)
    
    # if do plot heatmap
    if(plotHeatmapRaw==TRUE){
        pdf(paste0(output,".raw_heatmap.pdf"), width=5, height=6)
        heatmap.2(c.means, trace='none', scale='row', Colv=F, Rowv=F, dendrogram='none',
                  useRaster=T, col=colorRampPalette(c("dodgerblue4","white","firebrick4"))(100), labRow=F)
        dev.off()
    }
    #write.table(c.means, file=paste0(output,".RawClusterMeans.txt"),
    #            quote=F, row.names=T, col.names=T, sep="\t")


    # fit 2-guassian mixture model, return higher threshold
    message(" - fitting mixture model ...")
    thresholds <- apply(clust.means, 2, function(x){
        mod <- Mclust(x, G=2, verbose=F)
        top <- names(mod$parameters$mean)[which.max(mod$parameters$mean)]
        upper.cells <- x[which(mod$classification==top)]
        if(norm.clust.means){
            val <- min(upper.cells[upper.cells>0])
        }else{
            val <- quantile(upper.cells[upper.cells>0], c(0.05))
        }
        if(val == 0){
            val <- min(upper.cells[upper.cells>0])
        }
        return(val)
        #max(mod$parameters$mean)
    })

    # scale cell activity by cluster-specific threshold for expression
    message(" - scaling by cluster averages ...")
    adj.thresh <- thresholds[df$Cluster]
    print(head(adj.thresh, n=10))
    adj.act <- acts.o %*% Diagonal(x=1/adj.thresh)
	message(" - check adj.act 0 ...")
	colnames(adj.act) <- colnames(acts.o) #added by xuan to add column names, don't know if it is correct though.
	print(head(adj.act[,1:9]))
    adj.act@x[is.infinite(adj.act@x)] <- 0
    adj.act@x[is.na(adj.act@x)] <- 0
    adj.act@x <- round(adj.act@x, digits=0)
    adj.act <- adj.act[Matrix::rowSums(adj.act>0)>0,]
    adj.act <- adj.act[,Matrix::colSums(adj.act>0)>0]
	message(" - check adj.act 1 ...")
	print(head(adj.act[,1:11]))

    # un-normalized counts
    ua.out <- as.data.frame(summary(adj.act))
    ua.out$i <- rownames(adj.act)[ua.out$i]
    ua.out$j <- colnames(adj.act)[ua.out$j]
    ua.out <- ua.out[ua.out$x>0,]
    #write.table(ua.out, file=paste0(output,".countsActivity.sparse"),
    #            quote=F, row.names=F, col.names=F, sep="\t")

    # normalize by size factors
    if(logTransform==F){
        message(" - lib size before size factors = ")
        print(head(Matrix::colSums(adj.act)))
		message(" - check adj.act 2 ...")
		print(head(adj.act[,1:12]))
    }

    # do normalization
    message(" - estimating normalization factors ...")
    results <- normalizeSF(adj.act, verbose=T)
    norm.act <- results$norm.act
    norm.factor <- results$norm.factor

    # log transform?
    if(logTransform==T){
        message(" - square-root transforming counts activity ...")
        norm.act <- sqrt(norm.act)
        norm.act@x[is.na(norm.act@x)] <- 0
        norm.act@x[is.infinite(norm.act@x)] <- 0
        message(" - lib size afer square-root transformation = ")
        print(head(Matrix::colSums(norm.act)))
    }

    # write size factors to meta data file
    df$size_factors <- norm.factor[rownames(df)]
    #write.table(df, file=paste0(output,".size_factors.txt"),
    #            quote=F, row.names=T, col.names=T, sep="\t")

    # verbose
    message(" - lib size after size factors = ")
    print(head(Matrix::colSums(norm.act)))

    # remove empty cells/genes
    norm.act <- norm.act[Matrix::rowSums(norm.act>0)>0,]
    norm.act <- norm.act[,Matrix::colSums(norm.act>0)>0]
    message(" - cells = ",ncol(norm.act), " | genes = ", nrow(norm.act))
    print(head(norm.act[,1:10]))

    # print output to disk
    ia.out <- as.data.frame(summary(norm.act))
    ia.out$i <- rownames(adj.act)[ia.out$i]
    ia.out$j <- colnames(adj.act)[ia.out$j]
    ia.out <- ia.out[ia.out$x>0,]
    #write.table(ia.out, file=paste0(output,".normalizedActivity.sparse"), quote=F, row.names=F,
    #            col.names=F, sep="\t")

    # return
    message(" - returning normalized matrix ...")
    return(list(norm.act=norm.act, norm.factor=norm.factor, adj.act=adj.act, row.o=row.o))
}

# smooth data
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

# find cluster averages
clusterAves        <- function(df, 
                               acts, 
                               output, 
                               markers){

    # verbose start-up
    if(is.null(df) | is.null(acts)){
        stop("ERROR: must provide metadata and activity matrix")
    }
    if(! "Cluster" %in% colnames(df)){
        df$Cluster <- df$umapclust
    }
    rownames(markers) <- markers$geneID

    # estimate significant differences
    df <- df[colnames(acts),]
    clust.o <- mixedsort(unique(df$Cluster))
    df$Cluster <- factor(df$Cluster, levels=clust.o)

    # iterate over genes
    message(" - transversing gene activity ...")
    it <- 0
    df.acts <- as.data.frame(t(as.matrix(acts)))
    df.splits <- split(df.acts, df$Cluster)
    difs.mean <- lapply(df.splits, function(z){colMeans(z)})
    difs.mean <- as.matrix(do.call(cbind, difs.mean))

    # plot
    c.means <- difs.mean[,mixedorder(colnames(difs.mean))]
    c.means <- c.means[rownames(c.means) %in% rownames(markers),]
    markers <- markers[rownames(c.means),]
    row.o <- apply(c.means, 1, which.max)
    markers <-  markers[order(row.o, decreasing=F),]
    c.means <- c.means[order(row.o, decreasing=F),]

    # rowSide cols
    r.cols <- colorRampPalette(brewer.pal(8,"Set2"))(length(unique(markers$type)))

    # plot
    pdf(paste0(output,"_markers_heatmap.pdf"), width=10, height=12)
    heatmap.2(c.means, trace='none', scale='row', Colv=F, Rowv=F, dendrogram='none',
              useRaster=F, col=colorRampPalette(c("dodgerblue4","white","firebrick4"))(100),
              margins=c(10,10), cexRow = 0.7, RowSideColors = r.cols[factor(markers$type)],
              labRow=paste(markers$name,markers$type,sep="-"))
    dev.off()

    # estimate shannon
    sh <- apply(c.means, 1, shannon.entropy)
    total <- sum(sh)
    print(sh)
}



###################################################################################################
###################################################################################################
###################################################################################################

# plot activity scores
plot.act.scores    <- function(df, 
                               acts=acts, 
                               info=NULL, 
                               top=NULL,
                               logT=F,
                               marker.dist=NULL,
                               outname="markerActivityScores.pdf", 
                               lim=0.95){

    # prep data
    df <- df[rownames(df) %in% colnames(acts),]
    acts <- acts[,which(rownames(df) %in% colnames(acts))]

    # reorder rows
    rownames(info) <- info$geneID
    info <- info[order(info$type),]
    info.genes <- rownames(info)
    act.genes <- rownames(acts)
    rd.cells <- rownames(df)

    # common genes
    common <- intersect(info.genes, act.genes)
    info <- info[which(rownames(info) %in% common),]
    info.ordered <- rownames(info)
    sub.scores <- acts[info.ordered,]
    gids <- info.ordered

    # setup plot size
    nrows <- ceiling(length(gids)/6)
    totals <- nrows*6
    ratio <- nrows/6

    # params
    png(file=outname, width=12, height=ratio*12, units="in", res=500, type="cairo")
    layout(matrix(c(1:totals), ncol=6, byrow=T))
    par(mar=c(2,2,1,1))

    # adjust cluster IDs
    message("begin plotting pre-defined markers...")
    for (i in 1:length(gids)){

        # copy meta data
        gene.index <- which(rownames(sub.scores) == gids[i])
        acv <- sub.scores[gene.index,]

        # set up plot cols/sizes
        orderRow <- order(acv, decreasing=F)
        #cols <- colorRampPalette(c("grey75","grey75","goldenrod2","firebrick3"), bias=1)(100)
        #cols <- colorRampPalette(c("deepskyblue","goldenrod2","firebrick3"))(100)
        #cols <- inferno(100)
        #cols <- plasma(100)
        cols <- colorRampPalette(c("grey80","grey76","grey72",brewer.pal(9, "RdPu")[3:9]), bias=0.75)(100)
        acv <- as.numeric(acv[orderRow])
        if(logT==T){
            acv <- log2(acv+1)
        }
        df2 <- df[orderRow,]
        acv[is.na(acv)] <- 0
        acv[is.infinite(acv)] <- 0
        upper.lim <- quantile(acv, lim)
        acv[acv > upper.lim] <- upper.lim
        if(!is.null(marker.dist)){
            message(" - # cells = ", length(acv), "| min: ", marker.dist[[gids[i]]][1], " | max: ",marker.dist[[gids[i]]][2])
            colvec <- cols[cut(acv, breaks=seq(from=marker.dist[[gids[i]]][1], to=marker.dist[[gids[i]]][2], length.out=101))]
        }else{
            min.acv <- min(acv) - (1e-6*min(acv))
            max.acv <- max(acv) + (1e-6*max(acv))
            message(" - # cells = ", length(acv), "| min: ", min.acv, " | max: ",max.acv)
	    if(min.acv == max.acv){
		next
	    }
            colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv, length.out=101))]
        }
        colvec[is.na(colvec) & acv > mean(acv)] <- cols[length(cols)]
	colvec[is.na(colvec) & acv == 0] <- cols[1]
        sizes <- 0.2 #rescale(acv, c(0.25, 0.3))

        # # plot(original code)
        # plot(df2$umap1, df2$umap2, col=colvec,
             # main=info$name[i],
             # xlab="", ylab="", bty="n",
             # xaxt="n", yaxt="n", pch=16, cex=0.1)
		
		# plot(xuan fixed the missing gene id.)
        plot(df2$umap1, df2$umap2, col=colvec,
             main=gids[i],
             xlab="", ylab="", bty="n",
             xaxt="n", yaxt="n", pch=16, cex=0.3)

    }

    # turn device off
    dev.off()

}

# plot activity scores of de novo genes
plot.new.markers   <- function(df=NULL, 
                               acts=NULL, 
                               outname="ActivityScores.pdf",
                               row.o=NULL,
                               top=5, 
                               normT='prop.dif',
                               logT=F, 
                               threshold=0.01,
                               lim=0.95){

    # verbose start-up
    if(is.null(df) | is.null(acts)){
        stop("ERROR: must provide metadata and activity matrix")
    }
    if(! "Cluster" %in% colnames(df)){
        df$Cluster <- df$umapclust
    }

    # estimate significant differences
    df <- df[colnames(acts),]
    df$Cluster <- factor(df$Cluster)

    # iterate over genes
    message(" - transversing gene activity ...")
    it <- 0
    total.mean <- Matrix::rowMeans(acts)
    props <- acts
    props@x <- rep(1, length(props@x))
    total.prop <- Matrix::rowMeans(props)
    df.splits <- split(colnames(acts), df$Cluster)

    # estimate mean, proportions
    message(" - estimating difference in proportions and means ...")
    difs.prop <- lapply(df.splits, function(z){
        sub.act <- acts[,z]
        Matrix::rowMeans(sub.act >0)
    })
    difs.mean <- lapply(df.splits, function(z){
        sub.act <- acts[,z]
        rowMeans(sub.act)
    })
    difs.prop <- do.call(cbind, difs.prop)
    difs.mean <- do.call(cbind, difs.mean)
    colnames(difs.prop) <- names(df.splits)
    colnames(difs.mean) <- names(df.splits)

    # estimate pairwise log2 fold change
    message(" - estimating mean log2 fold changes")
    ave.difs.mean <- t(apply(difs.mean, 1, function(x){
        sapply(x, function(z){mean(log2((z+1)/(x+1)), na.rm=T)})
    }))

    # save output
    write.table(ave.difs.mean, file=paste0(outname,".NormalizedAveLog2FC.txt"),
                quote=F, row.names=T, col.names=T, sep="\t")
    write.table(difs.mean, file=paste0(outname,".NormalizedClusterMeans.txt"),
                quote=F, row.names=T, col.names=T, sep="\t")
    write.table(difs.prop, file=paste0(outname,".NormalizedClusterProportions.txt"),
                quote=F, row.names=T, col.names=T, sep="\t")

    # output options
    if(normT=="mean.dif"){
        difs <- ave.difs.mean
    }else if(normT=="prop.dif"){
        difs <- difs.prop
    }else if(normT=="adj.dif"){
        difs <- ave.difs.mean * difs.prop
    }
    difs <- difs[,mixedorder(colnames(difs))]

    # # plot normalized heatmap
    # c.means <- difs.mean[,mixedorder(colnames(difs.mean))]
    # pdf(paste0(output,".normalized_heatmap.pdf"), width=5, height=6)
    # heatmap.2(c.means, trace='none', scale='row', Colv=F, Rowv=F, dendrogram='none',
    #           useRaster=T, col=colorRampPalette(c("dodgerblue4","white","firebrick4"))(100), labRow=F)
    # dev.off()

    # filter out genes by proportion in at least one cluster
    message(" - selecting genes passing minimum proportion threshold ...")
    pass.genes <- apply(difs.prop, 1, function(x){length(x[x>threshold])})
    difs <- difs[pass.genes > 0,]

    # melt
    reduced <- melt(difs)
    reduced <- reduced[order(reduced$value, decreasing=T),]
    top.genes <- Reduce(rbind, by(reduced, reduced$Var2, head, n=top))
    top.genes.out <- subset(reduced, reduced$value >= 1)
    write.table(top.genes.out, file=paste0(outname,".log2FC1genes.txt"), quote=F, row.names=T, col.names=T, sep="\t")
    topN <- top.genes$Var1
    message(" - # genes = ", nrow(top.genes), " | # per cluster: ", top)
    print(head(top.genes))

    # plot top 50 -----------------------------------------------------------------------------
    nrows <- ceiling(length(topN)/top)
    totals <- nrows*top
    ratio <- nrows/top

    # params
    pdf(file=paste0(outname,".denovo.pdf"), width=20, height=ratio*20)
    layout(matrix(c(1:totals), ncol=top, byrow=T))
    par(mar=c(2,2,1,1))

    # adjust cluster IDs
    message("begin plotting ...")
    for (i in 1:nrow(top.genes)){

        # copy meta data
        dff <- df
        geneID <- top.genes$Var1[i]
        acv <- as.numeric(acts[rownames(acts) %in% geneID,])
        if(logT==T){
            acv <- log2(acv+1)
        }
        message(" - gene: ", geneID)

        # set up plot cols/sizes
        orderRow <- order(acv, decreasing=F)
        #cols <- colorRampPalette(c("grey85",brewer.pal(9,"YlGnBu")), bias=1)(100)
        #cols <- colorRampPalette(c("deepskyblue","goldenrod2","firebrick3"))(100)
        #cols <- inferno(100)
        #cols <- plasma(100)
        cols <- colorRampPalette(viridis(100), bias=0.5)(100)
        #cols <- colorRampPalette(c("grey75","darkorange","firebrick3"), bias=0.75)(100)
        acv <- acv[orderRow]
        df2 <- df[orderRow,]

        # set up upper limits
        upper.lim <- quantile(acv, lim)
        acv[acv > upper.lim] <- upper.lim

        min.acv <- -0.1
        max.acv <- max(acv) + (0.05*max(acv))
        colvec <- cols[cut(acv, breaks=seq(min.acv, max.acv,length.out=101))]
        colvec[is.na(colvec)] <- cols[1]
        sizes <- rescale(acv, c(0.35, 0.4))

        # plot
        plot(df2$umap1, df2$umap2,
             col=colvec,
             main=paste(top.genes$Var1[i],top.genes$Var2[i],sep="-"),
             xlab="", ylab="", bty="n", xaxt="n", yaxt="n", pch=16, cex=0.5)

    }

    # turn device off
    dev.off()

}



###################################################################################################
###################################################################################################
###################################################################################################

# parallel implementation over major clusters
runMajorPriori     <- function(all.b, 
                               all.activity,
                               all.hpcs,
                               marker.info, 
                               threads=1,
                               output="all",
                               smooth.markers=F,
                               plotHeatmaps=F){

#out <- runMajorPriori(b.meta, activity.all, h.pcs1, marker.info.dat, threads=threads, output=prefix, smooth.markers=T)


    # align cell barcodes
    all.activity <- all.activity[,rownames(all.b)] #all.activity <- activity.all[,rownames(b.meta)]
    all.activity <- all.activity[Matrix::rowSums(all.activity)>0,]
    all.activity <- all.activity[,Matrix::colSums(all.activity)>0]
    ids <- intersect(rownames(all.b), colnames(all.activity)) #ids <- intersect(rownames(b.meta), colnames(all.activity))
    ids <- intersect(ids, rownames(all.hpcs))#ids <- intersect(ids, rownames(h.pcs1))
    all.b <- all.b[ids,] #b.meta <- b.meta[ids,]
    all.activity <- all.activity[,ids]
    all.hpcs <- all.hpcs[ids,]#h.pcs1 <- h.pcs1[ids,]

    # normalize per cell activity by cluster average and size factors
    results <- normalize.activity(all.b, 
                                  all.activity, 
                                  output=output, 
                                  logTransform=F, 
                                  scaleP=F, 
                                  plotHeatmapRaw=F)
								  
	# results <- normalize.activity(b.meta, 
                                  # all.activity, 
                                  # output=output, 
                                  # logTransform=F, 
                                  # scaleP=F, 
                                  # plotHeatmapRaw=F)
								  
    activity <- results$norm.act
    row.o <- results$row.o
    rm(results)
    gc()
    
    # line-up ids ###################################No problem here
    barcode.ids <- intersect(colnames(activity), rownames(all.hpcs))#barcode.ids <- intersect(colnames(activity), rownames(all.hpcs))
    barcode.ids <- intersect(barcode.ids, rownames(all.b))
    activity <- activity[,barcode.ids]
    all.hpcs <- all.hpcs[barcode.ids,]
    all.b <- all.b[barcode.ids,]
    
    # if smooth genes only
    if(smooth.markers){
        activity <- activity[rownames(activity) %in% as.character(marker.info$geneID),]
    }

    # impute gene accessibility scores
    impute.activity <- smooth.data(activity, 
                                   k=25,
                                   step=3, 
                                   npcs=ncol(all.hpcs), 
                                   df=NULL,
                                   rds=as.data.frame(all.hpcs), 
                                   threads=threads)
    
    # reformat sparse matrix
    #ia <- as.data.frame(summary(impute.activity))
    #ia$i <- rownames(impute.activity)[as.numeric(ia$i)]
    #ia$j <- colnames(impute.activity)[as.numeric(ia$j)]
    #write.table(ia, file=paste0(output,".smoothed.sparse"), quote=F, row.names=F, col.names=F, sep="\t")

    # ranges
    marker.impact <- impute.activity[rownames(impute.activity) %in% as.character(marker.info$geneID),]
    marker.ranges <- lapply(rownames(marker.impact), function(x){
        x[is.na(x)] <- 0
        min.x <- min(marker.impact[x,], na.rm=T)
        max.x <- max(marker.impact[x,], na.rm=T)
        rng <- c(min.x,max.x)
        rng[is.na(rng)] <- 0
        rng[1] <- rng[1] - 0.1
        rng[2] <- rng[2] + (0.05*rng[2])
        return(rng)
    })
    names(marker.ranges) <- rownames(marker.impact)

    # plot all
    #plot.act.scores(all.b, acts=activity,
    #                info=marker.info,
    #                logT=T,
    #                lim=0.99,
    #                marker.dist=NULL,
    #                outname=paste0("combined.",output,".normalized.known.Markers.png"))

    plot.act.scores(all.b, acts=impute.activity,
                    info=marker.info,
                    logT=F,
                    lim=0.999,
                    marker.dist=NULL,
                    outname=paste0("combined.",output,".impute.known.Markers.png"))

    # return
    return(list(b=all.b, activity=activity, impute.activity=impute.activity))

}

# run major de novo gene markers
runMajorDeNovo     <- function(b, 
                               activity,
                               impute.activity, 
                               output="all",
                               marker.info){

    # verbose
    message(" - finding de novo enriched markers ...")

    # find de novo markers
    plot.new.markers(df=b,
                     acts=activity,
                     top=5,
                     normT='mean.dif',
                     logT=T,
                     lim=0.99,
                     outname=paste0(output,".normalized"))

    plot.new.markers(df=b,
                     acts=impute.activity,
                     top=5,
                     normT='mean.dif',
                     logT=F,
                     lim=0.99,
                     outname=paste0(output,".imputed"))

    # cluster aves
    clusterAves(b, activity, paste0(output,".normalized"), marker.info)
    clusterAves(b, impute.activity, paste0(output,".imputed"), marker.info)

}


