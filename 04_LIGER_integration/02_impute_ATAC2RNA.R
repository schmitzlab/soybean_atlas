## impute ATAC data onto RNA cells ##
#module load R/4.3.1-foss-2022a
# load libraries
library(rliger)

#reference.
#https://rdrr.io/github/MacoskoLab/liger/f/vignettes/Integrating_scRNA_and_scATAC_data.Rmd#google_vignette

liger <- readRDS("Gm_atlas_All_endosperm.iNMF_liger_obj.rds")
atacp <- readRDS("./Gm_atlas_All_endosperm_cpm4_ACR_sparse.rds")
motif <- read.table("./05_Motif_deviation/Gm_atlas_All_endosperm_cpm4.ACR.sorted.bed.motif.deviations.txt")
prefix <- "Gm_atlas_All_endosperm"


# impute ATAC on RNA
ids <- colnames(liger@raw.data$atac)
liger@raw.data$atac <- atacp[,ids]
liger.atac.imp <- imputeKNN(liger, reference='atac', queries=list('rna'))

# add motif to raw data
liger@raw.data$atac <- Matrix(t(motif[ids,]), sparse=T)
liger.motif.imp <- imputeKNN(liger, reference='atac', queries='rna')

# retrieve imputed
impATAC <- liger.atac.imp@norm.data[['rna']]
impMotif <- liger.motif.imp@norm.data[['rna']]
obsRNA <- liger@norm.data[['rna']]

# save output
saveRDS(impATAC, file=paste0(prefix,"_liger_ACRs_mapto_RNA.rds"))
saveRDS(impMotif, file=paste0(prefix,"_liger_Motif_deviation_mapto_RNA.rds"))
saveRDS(obsRNA, file=paste0(prefix,"_liger_RNA.rds"))