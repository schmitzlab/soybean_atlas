library("here")
library(devtools)
library(Seurat)
library(Socrates)


###############################################################
#                  Input parameters
###############################################################
#bed <- system.file("extdata", "test.tn5.bed.gz", package = "Socrates")
#ann <- system.file("extdata", "gencode.v19.annotation.gff3.gz", package = "Socrates")
#chr <- system.file("extdata", "hg19.txt", package = "Socrates")

args=commandArgs(T)
bed <- args[1]
ann <- args[2]
chr <- args[3]
min_tn5 <- as.numeric(args[4])
pre1 <- args[5]
gsize <- as.numeric(args[6])


###############################################################
#                  run QC
###############################################################

# load data
obj <- loadBEDandGenomeData(bed, ann, chr, attribute="Parent")

# count organellar reads
obj <- countRemoveOrganelle(obj, org_scaffolds=c("ChrCp", "ChrMt"), remove_reads=T)

# call ACR
obj <- callACRs(obj, genomesize=gsize, 
                shift= -75, 
                extsize=150,
                fdr=0.1,
                output=paste0(pre1,"_peaks"), 
                tempdir=paste0("./",pre1,"_macs2_temp"),
                verbose=T)

# build metadata
obj <- buildMetaData(obj, tss.window=1000, verbose=TRUE)

# generate sparse matrix
obj <- generateMatrix(obj,
                      filtered=F, 
                      windows=500, 
                      peaks=F, 
                      verbose=T)

# convert to Socrates format for downstream analysis. 
soc.obj <- convertSparseData(obj, verbose=T)

# save QC object
saveRDS(soc.obj, file=paste0(pre1,".raw.soc.rds"))


# filter cell

obj <- findCells(obj, 
                 doplot=T,
                 max.cells=15000,
                 set.tn5.cutoff=min_tn5,
                 filt.tss=T, 
                 filt.frip=T,
				 frip.min.freq=0.2,
                 prefix = pre1)



write.table(obj$meta.v3,file=paste0(pre1, ".filtered.metadata.txt"),quote=F,sep = "\t")
