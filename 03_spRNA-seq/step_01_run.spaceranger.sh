#!/bin/bash
#SBATCH --partition=
#SBATCH --output=spaceranger.%j.out
#SBATCH --job-name=space
#SBATCH --ntasks=1
#SBATCH --time=90:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=30
#SBATCH --mail-user 

date;hostname;pwd

ml SpaceRanger
spaceranger count --id=gma_CS2_A \
               --transcriptome=Gm_scRNA-seq \
               --fastqs=./rawdata/Ziliang_Gm_spatial \
               --sample=sp_Gm_seed_CS2_A \
               --image=./images/CS2ndBatch/A1.jpg \
               --slide=V12N18-405 \
               --area=A1 \
               --loupe-alignment=./images/CS2ndBatch/V12N18-405-A1.json \
               --localcores=20 \
               --localmem=120 \
