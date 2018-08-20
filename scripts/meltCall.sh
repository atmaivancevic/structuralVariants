#!/bin/bash  

## Script for calling MEIs with MELT
## Date: 17 Aug 2017
##
## Example usage:
## INDIR=/data/neurogenetics/alignments/Illumina/genomes/allGenomes sbatch --array 0-1 meltCall.sh

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=04:00:00
#SBATCH --mem=24GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# define key variables
export FASTDIR=/fast/users/$USER

OUTDIR=$FASTDIR/outputs
list=/data/neurogenetics/RefSeq/MELT/mei_list.txt
genes=/data/neurogenetics/RefSeq/MELT/hg19.genes.bed
ref=/data/neurogenetics/RefSeq/ucsc.hg19.fasta

# define query bam files
QUERIES=($(ls $INDIR/*.bam | xargs -n 1 basename))

# load modules 
module load Bowtie2/2.2.9-foss-2016b
module load Java/1.8.0_121
module load MELT/2.1.5-Java-1.8.0_121

# run the thing

### SV discovery/calling phase ###
echo $(date +"[%b %d %H:%M:%S] Starting MELT")
echo "Processing file: "${QUERIES[$SLURM_ARRAY_TASK_ID]}

java -jar $EBROOTMELT/MELT.jar \
Single -a \
-bamfile ${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]} \
-h $ref -w $OUTDIR/MEIcalling/MELT \
-t $list \
-n $genes \
-c 30 \
2> ${FASTDIR}/slurmLOG/meltCalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.log

echo $(date +"[%b %d %H:%M:%S] All done!")
