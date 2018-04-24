#!/bin/bash

## Script for calling SVs with manta
## Use this script if you have trio data
## Date: 24 April 2018 
##
## Example usage:
## INBAM1=/data/neurogenetics/alignments/Illumina/genomes/CP/FR07959117/FR07959117.5.bam INBAM2=/data/neurogenetics/alignments/Illumina/genomes/CP/FR07959116/FR07959116.5.bam INBAM3=/data/neurogenetics/alignments/Illumina/genomes/CP/FR07959115/FR07959115.5.bam PREFIX=FR07959117_FR07959116_FR07959115 sbatch mantaTrioCall.sh

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1-00:00
#SBATCH --mem=5GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# define key variables
export FASTDIR=/fast/users/$USER

MANTACONFIG=/data/neurogenetics/executables/manta-1.1.1.centos5_x86_64/bin/configManta.py
OUTDIR=$FASTDIR/outputs
GENOMEDIR=$FASTDIR/references/genomes

# Run the thing

echo $(date +"[%b %d %H:%M:%S] Make output dir")
mkdir -p ${OUTDIR}/SVcalling/mantaOut/$PREFIX

echo $(date +"[%b %d %H:%M:%S] Set up manta workflow")
python $MANTACONFIG \
--bam $INBAM1 \
--bam $INBAM2 \
--bam $INBAM3 \
--referenceFasta $GENOMEDIR/ucsc.hg19.fasta \
--runDir ${OUTDIR}/SVcalling/mantaOut/$PREFIX

echo $(date +"[%b %d %H:%M:%S] Run manta workflow")
python ${OUTDIR}/SVcalling/mantaOut/$PREFIX/runWorkflow.py -m local -j 16

echo $(date +"[%b %d %H:%M:%S] All done!")

