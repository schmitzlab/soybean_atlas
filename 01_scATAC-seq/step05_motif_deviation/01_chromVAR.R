###################################################################################################
###                             chromVAR analysis of scATAC data                                ###
###################################################################################################

# arguments
args <- commandArgs(TRUE)
if(length(args) != 6){stop("Rscript chromVAR.R <threads> <sparseMatrix.rds> <meta> <peaks.bed> <fai> <output_prefix>")}
threads <- as.numeric(args[1])
input.sp <- as.character(args[2])
metadata <- as.character(args[3])
peakFile <- as.character(args[4])
FAI <- as.character(args[5])
prefix <- as.character(args[6])

#load libraries
library(chromVAR)
library(motifmatchr)
library(BiocParallel)
library(BSgenome.Gmax.a4.v1) 
library(Matrix)
library(SummarizedExperiment)
library(GenomicAlignments)
library(dplyr)
library(TFBSTools)
library(JASPAR2022)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(stringr)

# functions
loadPeaks <- function(x, y, peaks, extra_cols=4){

	# create ref
	fai <- lapply(as.character(y$V1), function(z){
		return(as.numeric(as.character(y$V2)[y$V1==z]))
	})
	names(fai) <- as.character(y$V1)

	# load bed
	bed <- as.data.frame(do.call(rbind, strsplit(rownames(x),"_")))	
	newbed <- read.table(peaks)
	
	rownames(newbed) <- paste(newbed$V1,newbed$V2,newbed$V3,sep="_")
	#newbed <- subset(newbed, newbed$V4 != "exons")
	rownames(bed) <- rownames(x)
	bed <- bed[rownames(newbed),]
	rownames(bed) <- seq(1:nrow(bed))

	# convert 2 GR
	colnames(bed) <- c("chr", "start", "end")
	bed$chr <- as.character(bed$chr)
	bed$start <- as.numeric(as.character(bed$start))
	bed$end <- as.numeric(as.character(bed$end))
	bed$keep <- ifelse(bed$start > fai[bed$chr] | bed$end > fai[bed$chr], 0, 1)
	x <- x[bed$keep > 0, ]
	bed <- bed[bed$keep > 0,]
	bed$keep <- NULL
	bed[, "start"] <- bed[, "start"]
	bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = F)

	# sort
	sorted_bed <- sortSeqlevels(bed)
	sorted_bed <- sort(sorted_bed, ignore.strand = TRUE)
	sbeddf <- as.data.frame(sorted_bed)
	s.ids <- paste(sbeddf$seqnames,sbeddf$start,sbeddf$end,sep="_")
	shared <- intersect(s.ids, rownames(x))
	x <- x[shared,]
	sorted_bed <- subset(sorted_bed, c(s.ids %in% shared))
	
	return(list(bed=sorted_bed, cnts=x))
}
getJasparMotifs2 <- function(species = "Homo sapiens", collection = "CORE", ...){
    opts <- list()
    opts["species"] <- species
    opts["collection"] <- collection
    opts <- c(opts, list(...))
    out <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)
    if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
        names(out) <- paste(names(out), TFBSTools::name(out),
            sep = "_")
    return(out)
}

# set number of cores
register(MulticoreParam(threads))

# verbose
message("########################################")
message("########################################")
message("")
message("============================")
message("     running chromVAR       ")
message("============================")
message("")


###################################################################################################
### load and process data									   
###################################################################################################

# build counts matrix
message("Loading count matrix ...")
a <- readRDS(input.sp)

# input files
message("Loading peak information ...")
ref <- read.table(FAI)
obj <- loadPeaks(a, ref, peaks=peakFile)
peaks <- obj$bed
a <- obj$cnts

# load meta.data
message("Loading meta data ...")
meta <- read.table(metadata,comment.char = "")
a <- a[,colnames(a) %in% rownames(meta)]
meta <- meta[colnames(a),]
meta$depth <- Matrix::colSums(a)
message("cells = ",ncol(a), " | peaks = ", nrow(a))

rownames(a) <- NULL

# create frag counts object
message("Creating experiment object ...")
fragment_counts <- SummarizedExperiment(assays = list(counts = a),
                                        rowRanges = peaks,
                                        colData = meta)

# add GC data
message("Estimating GC bias ...")
fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.Gmax.a4.v1)

# filter cells
message("Filtering samples ...")
filtered_counts <- filterSamples(fragment_counts, min_depth=100, min_in_peaks=0.1, shiny=F)

# filter peaks
message("Filtering peaks ...")
filtered_counts <- filterPeaks(filtered_counts, non_overlapping=T, min_fragments_per_peak=10) #Here I give it the final filtered peaks so setting a small number for min_fragments_per_peak 

###############################################################################
## motif deviation
###############################################################################

# estimate deviations
message("Running motif analysis ...")
jaspmotifs.at       <- getJasparMotifs2(species = "Arabidopsis thaliana")
jaspmotifs.zm       <- getJasparMotifs2(species = "Zea mays")
jaspmotifs.os       <- getJasparMotifs2(species = "Oryza sativa")
jaspmotifs.gm       <- getJasparMotifs2(species = "Glycine max")
jaspmotifs <- c(jaspmotifs.at, jaspmotifs.zm, jaspmotifs.os, jaspmotifs.gm)

motif            <- matchMotifs(jaspmotifs, filtered_counts, genome = BSgenome.Gmax.a4.v1)
dev.motif        <- computeDeviations(object = filtered_counts, annotations = motif)
dev.motif.scores <- deviationScores(dev.motif)
motif.devs       <- deviations(dev.motif)
saveRDS(motif, file=paste0(prefix,".motif_matches.rds"))
write.table(t(dev.motif.scores), file=paste0(prefix,".motif.scores.txt"), quote=F, row.names=T, col.names=T, sep="\t")
write.table(t(motif.devs), file=paste0(prefix,".motif.deviations.txt"), quote=F, row.names=T, col.names=T, sep="\t")

## background peaks
bbpeaks <- getBackgroundPeaks(filtered_counts)
write.table(bbpeaks, file=paste0(prefix,".backgroundPeaks.mat.txt"), quote=F, row.names=T, col.names=T, sep="\t")
message("--Finished--")
