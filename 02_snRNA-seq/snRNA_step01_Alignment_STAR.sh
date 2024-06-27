#!/bin/bash
#SBATCH --partition=batch
#SBATCH --output=starsolo..%j.out
#SBATCH --job-name=SCcount
#SBATCH --ntasks=1
#SBATCH --time=60:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=30
#SBATCH --mail-user 
#SBATCH --mail-type=FAIL,END

date;hostname;pwd

read1=Gm_hp2_S8_L003_R1_001.fastq.gz
read2=Gm_hp2_S8_L003_R2_001.fastq.gz

ml STAR

STAR --runThreadN 30 --readFilesCommand zcat \
    --genomeDir ./genome/glycine_max_v4_MtCp/wm82v4_1/gm_star_index/ \
    --outFileNamePrefix hp2 --readFilesIn $read2 $read1 \
    --sjdbGTFfile ./genome/glycine_max_v4_MtCp/wm82v4_1/Gmax_v4.0_MtCp.gtf \
    --soloType Droplet --soloCBwhitelist ./genome/glycine_max_v4_MtCp/wm82v4_1/3M-february-2018.txt \
    --soloUMIlen 12 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloCellFilter EmptyDrops_CR \
    --soloFeatures GeneFull \
    --soloMultiMappers PropUnique \
    --soloBarcodeReadLength 151 \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 200000000000 \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
