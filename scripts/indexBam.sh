#!/bin/bash

## Script for indexing bam files
## Date: 28 March 2018 
##
## Example usage:
## INDIR=/data/neurogenetics/alignments/Illumina/genomes/allGenomes_chr10 sbatch --array 0-125 indexBam.sh

#SBATCH -J indexBam
#SBATCH -o /fast/users/$USER/slurmOUT/indexBam/slurm-%A_%a.out

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=1-00:00
#SBATCH --mem=8GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# define query bam files
QUERIES=($(ls $INDIR/*.bam | xargs -n 1 basename))

# load modules
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25

# run the thing

echo $(date +"[%b %d %H:%M:%S] Get ready to start indexing")
echo "BAM file: "${QUERIES[$SLURM_ARRAY_TASK_ID]}

samtools index ${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]}

echo $(date +"[%b %d %H:%M:%S] Indexing completed successfully")

