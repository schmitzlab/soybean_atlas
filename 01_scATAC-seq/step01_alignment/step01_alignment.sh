#!/bin/bash
#SBATCH --job-name=scATAC-seq-pipe               #Job name (testBowtie2)
#SBATCH --partition=batch         #Queue name (batch)
#SBATCH --nodes=1                 # Run all processes on a single node
#SBATCH --ntasks=1                #Run in a single task on a single node
#SBATCH --cpus-per-task=32         # Number of CPU cores per task (8)
#SBATCH --mem=250G                 # Job memory limit (10 GB)
#SBATCH --export=ALL              # Do not load any users<U+0092> explicit environment variables
#SBATCH --output=%x_%j.out        # Standard output log, e.g., testBowtie2_1234.out
#SBATCH --error=%x_%j.err         # Standard error log, e.g., testBowtie2_1234.err
#SBATCH --mail-type=END,FAIL      # Mail events (BEGIN, END, FAIL, ALL)

set -e
set -u
set -o pipefail

cd $SLURM_SUBMIT_DIR
#cellranger-atac mkref cellRanger --config=./Arabidopsis_TAIR10_config.txt
module load CellRanger-ATAC/2.1.0
module load ucsc/443
module load picard/2.27.5-Java-15
module load R/4.3.1-foss-2022a
#module load Python/3.7.4-GCCcore-8.3.0 
module load MACS2/2.2.7.1-foss-2021b
module load HarfBuzz/5.3.1-GCCcore-12.2.0
module load FriBidi/1.0.12-GCCcore-12.2.0
module load SAMtools/1.16.1-GCC-11.3.0
export LC_ALL=C

###################################
#parameter for cellranger-atac count
cellranger_ref=/scratch/xz24199/GenomeRef/Glycine_max_v4_MtCp/Gm_cellranger_v2
fastq_path=/scratch/xz24199/00libs/AH_1/scATAC-seq/AH1_Gm_hp3
base=AH1_Gm_Hypocotyl-rep3

###################################
#parameter for process bam
chrom_info=Glycine_max_v4_MtCp_chrom.sizes
threads=32
memory=250
qual=30

if [ ! -d "./01cellranger" ]; then
  mkdir ./01cellranger
fi
cd ./01cellranger

cellranger-atac count --id=$base \
                   --reference=$cellranger_ref \
                   --fastqs=$fastq_path \
                   --localcores=$threads \
                   --localmem=$memory
cd ..
#process bam
if [ ! -d "./02process_bam" ]; then
  mkdir ./02process_bam
fi
cd 02process_bam

        echo "filtering bam based on mapQ..."
        echo "retaining only mapped reads ..."
	    samtools view -@ $threads -bhq $qual -f 3 ../01cellranger/${base}/outs/possorted_bam.bam > ${base}.mq${qual}.bam
		
        # run picard
        echo "removing dups - ${base} ..."
        java -Xmx120g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        INPUT=${base}.mq${qual}.bam \
        OUTPUT=${base}.mq${qual}.rmdup.bam \
        METRICS_FILE=${base}.rmdup.metrics \
        BARCODE_TAG=CR \
        ASSUME_SORT_ORDER=coordinate \
        REMOVE_DUPLICATES=true \
        USE_JDK_DEFLATER=true \
        USE_JDK_INFLATER=true

        samtools stats ${base}.mq${qual}.rmdup.bam > ${base}.mq${qual}.rmdup.bam.stats

        echo "fixing barcode and filter Cp Mt reads..."
        perl ../fixBC.pl ${base}.mq${qual}.rmdup.bam | samtools view -bhS - > ${base}.BC.mq${qual}.rmdup.bam
        samtools stats ${base}.BC.mq${qual}.rmdup.bam > ${base}.BC.mq${qual}.rmdup.bam.stats

        # make Tn5 bed files
        echo "making Tn5 bed files ..."
        perl ../makeTn5bed.pl ${base}.BC.mq${qual}.rmdup.bam | sort -k1,1 -k2,2n - > ${base}.tn5.mq${qual}.bed
		
		cd ../
