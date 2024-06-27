## analyze branch-defining TFs ##

# load libraries
library(RColorBrewer)
library(gplots)
library(scales)

# load functions
loadTests <- function(){
  
  # load data
  maindir <- "/scratch/apm25309/single_cell/ATACseq/embryo_development_soybean/step04_pseudotime/"
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
kneePlot <- function(mydata, crange=2:15){
  
  wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in crange) wss[i] <- sum(kmeans(mydata,
                                       centers=i)$withinss)
  plot(1:max(crange), wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  
}

# load data
rna <- readRDS("PTSC.RNA.03.12.2024.rds")
motif <- readRDS("PTSC.Motif.03.12.2024.rds")
gene.ann <- read.delim("/scratch/apm25309/single_cell/ATACseq/embryo_development_soybean/data/annotations/soybean_arabidopsis.protein_alignment.top_hit.geneName.txt", header=F)
motif.ann <- read.table("/scratch/apm25309/single_cell/ATACseq/embryo_development_soybean/data/annotations/JASPAR_2024_motif_data.v1.txt", header=T)
tests <- loadTests()

# get shared TFs
shared <- intersect(motif.ann$geneID, gene.ann$V14)
motif.ann <- motif.ann[motif.ann$geneID %in% shared,]
gene.ann.TF <- gene.ann[gene.ann$V14 %in% shared,]
gene.ann.TF$V15 <- ifelse(is.na(gene.ann.TF$V15), gene.ann.TF$V14, gene.ann.TF$V15)

# rename lists
names(rna) <- gsub("_pseudotime\\.trajectory\\.txt","",names(rna))
names(motif) <- gsub("_pseudotime\\.trajectory\\.txt","",names(motif))
names(rna) <- gsub("RNA_","",names(rna))
names(motif) <- gsub("RNA_","",names(motif))

# get gene orders
tests$branch <- factor(tests$branch, levels=c("Axis_Parenchyma","Cotyledon_Parenchyma",
                                              "Epidermis","Provascular","Apical_Meristem"))
tests <- tests[order(tests$branch,decreasing=F),]
its <- 0
row.ids <- lapply(names(rna), function(z){
  its <<- its + 1
  zz <- rna[[z]]
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

# scale RNA
mrna <- do.call(cbind, rna)
mmotif <- do.call(cbind, motif)
mrna <- t(apply(mrna, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
mmotif <- t(apply(mmotif, 1, function(x){(x-mean(x, na.rm=T))/sd(x, na.rm=T)}))
mrna <- mrna[row.ids,]
raxis <- mrna[,1:500]
rcoty <- mrna[,501:1000]
repid <- mrna[,1001:1500]
rprov <- mrna[,1501:2000]
rmeri <- mrna[,2001:2500]
s.rna <- list(raxis,rcoty,repid,rprov,rmeri)
names(s.rna) <- names(rna)

# compare gene expression between branches
branches <- names(s.rna)
gene.cors <- lapply(seq(1:(length(branches)-1)), function(x){
  comp <- lapply(seq(from=(x+1), to=length(branches)), function(y){
    bID1 <- branches[x]
    bID2 <- branches[y]
    message(" - getting correlations between ", bID1, " and ", bID2, " | idx ", x, " and ",y)
    cors <- diag(cor(t(s.rna[[branches[x]]]), t(s.rna[[branches[y]]])))
    names(cors) <- rownames(s.rna[[1]])
    return(cors)
  })
  names(comp) <- paste(rep(branches[x], length(comp)), branches[(x+1):length(branches)], sep=":")
  comp <- do.call(cbind, comp)
  return(comp)
})
all <- do.call(cbind, gene.cors)

# plot axis vs. cotyledon
a.c <- all[1:(2456+4074),1]
a.c[is.na(a.c)] <- 0
cols <- colorRampPalette(rev(brewer.pal(9, "Spectral")))(100)

pdf("gene_correlations.axis.cotyledon.par.pdf", width=10, height=5)
sizes <- rep(0.5, length(a.c))
sizes[4844] <- 10
plot(a.c, pch=16, col=cols[cut(a.c, breaks=101)],
     cex=sizes)
fit <- smooth.spline(seq(1:length(a.c)), a.c)
lines(fit, col="grey75", lwd=2)
abline(v=2456, col="black")
dev.off()

# find branch TFs when expression diverges
axis.ranks <- which(a.c[1350:2456] < 0)
coty.ranks <- which(a.c[4600:length(a.c)] < 0)
axis.tfs <- as.data.frame(axis.ranks[names(axis.ranks) %in% gene.ann.TF$V13])
coty.tfs <- as.data.frame(coty.ranks[names(coty.ranks) %in% gene.ann.TF$V13])
colnames(axis.tfs) <- "rank"
colnames(coty.tfs) <- "rank"
ids <- gene.ann.TF$V15
names(ids) <- gene.ann.TF$V13
axis.tfs$AtID <- ids[rownames(axis.tfs)]
coty.tfs$AtID <- ids[rownames(coty.tfs)]
axis.tfs$branch <- "axis"
coty.tfs$branch <- "cotyledon"
df <- rbind(axis.tfs, coty.tfs)
write.table(df, file="divergent_TFs.axis.cotyledon.txt", quote=F, row.names=T, col.names=T, sep="\t")


# correlate transcriptome and motif deviation patterns
TF.trx.cor <- cor(t(mrna), t(mmotif))
TF.pt.trx <- cor(t(mrna[row.ids,]), t(mmotif))

# cluster TF motifs
kneePlot(t(TF.pt.trx), crange=2:20)
centers <- 8
tf.modules <- kmeans(t(TF.pt.trx), centers=centers)
gene.modules <- kmeans(TF.pt.trx, centers=centers)
c.o <- unlist(lapply(seq(1:centers), function(z){
  
  ids <- names(tf.modules$cluster)[tf.modules$cluster==z]
  sub.df <- TF.pt.trx[,ids]
  col.o <- hclust(as.dist(1-cor(sub.df)))$order
  ids[col.o]
  
}))
TF.pt.trx <- TF.pt.trx[,c.o]
col.cols <- brewer.pal(centers, "YlGnBu")
col.col <- col.cols[tf.modules$cluster[c.o]]

pdf("gene_motif_modules.PT.kmeans_TFs.03.19.2024.pdf", width=10, height=10)
heatmap.2(TF.pt.trx, trace='none', Rowv=T, Colv=F, dendrogram="row",
          ColSideColors = col.col,
          col=colorRampPalette(c(rev(brewer.pal(9, "Blues")),"white",hcl.colors(9, "YlOrRd", rev = T)))(100),
          labRow=NA, labCol=NA,
          useRaster=T)
dev.off()

# AtHB13 compare between TF motif deviation and gene expression
motif.df <- as.data.frame(do.call(rbind, strsplit(motif.ann$AC,"\\.")))
motif.ann$motifID <- as.character(motif.df$V1)
df <- as.data.frame(do.call(rbind, strsplit(rownames(mmotif),"_")))
ddf <- as.data.frame(do.call(rbind, strsplit(as.character(df$V1),"\\.")))
jaspar <- cbind(ddf$V1,df,rownames(mmotif))
colnames(jaspar) <- c("motifID","AC", "ID", "jaspID")
sharedIDs <- intersect(jaspar$motifID, motif.ann$motifID)
rownames(jaspar) <- jaspar$motifID
rownames(motif.ann) <- motif.ann$motifID
motif.ann <-  motif.ann[sharedIDs,]
jaspar <- jaspar[sharedIDs,]
jaspar$AtID <- motif.ann$geneID

mrna <- t(apply(mrna, 1, rescale, c(-1,1)))
mmotif <- t(apply(mmotif, 1, rescale, c(-1,1)))

tfid <- "Glyma.08G298400"
geneID <- "AT1G69780"
jaspID <- "MA1212.1_ATHB.13"
tf.exp <- mrna[tfid,501:1000]
pdf("ATHB13_dynamics.pdf", width=5, height=5)
plot(seq(1:500),mrna[tfid,501:1000], type="l", col="dodgerblue4",
     xlab="Pseudotime", ylab="Scaled gene expression/motif deviation")
lines(seq(1:500),mmotif[jaspID,501:1000], col="firebrick4")
grid(lty=1)
dev.off()

