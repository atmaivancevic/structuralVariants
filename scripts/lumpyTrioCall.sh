#!/bin/bash

## Script for calling SVs with lumpy
## Use this script if you have a trio of samples
## Date: 24 April 2018 
##
## Example usage:
## INBAM1=/data/neurogenetics/alignments/Illumina/genomes/CP/FR07959117/FR07959117.5.bam INBAM2=/data/neurogenetics/alignments/Illumina/genomes/CP/FR07959116/FR07959116.5.bam INBAM3=/data/neurogenetics/alignments/Illumina/genomes/CP/FR07959115/FR07959115.5.bam PREFIX=FR07959117_FR07959116_FR07959115 sbatch lumpyTrioCall.sh

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=3-00:00
#SBATCH --mem=24GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# define key variables
export FASTDIR=/fast/users/$USER

SPEEDSEQEXE=/data/neurogenetics/executables/speedseq/bin/speedseq
LUMPYDIR=/data/neurogenetics/executables/speedseq/src/lumpy-sv
OUTDIR=$FASTDIR/outputs
GENOMEDIR=$FASTDIR/references/genomes

# load modules
module load LUMPY/0.2.13-foss-2016b
module load SAMtools/1.3.1-foss-2016b 
module load Pysam/0.10.0-foss-2016uofa-Python-2.7.12

# run pipeline
# note: running lumpy via speedseq
# to auto-generate discordant/split bams
# and then genotype with svtyper

### Pre-processing ###

# samples have already been aligned to ref genome
# need to extract the discordant paired-end alignments for each sample
echo $(date +"[%b %d %H:%M:%S] Extract discordant paired-end alignments")
samtools view -b -F 1294 $INBAM1 > ${INBAM1%.bam}.discordants.unsorted.bam
samtools view -b -F 1294 $INBAM2 > ${INBAM2%.bam}.discordants.unsorted.bam
samtools view -b -F 1294 $INBAM3 > ${INBAM3%.bam}.discordants.unsorted.bam
echo $(date +"[%b %d %H:%M:%S] Done")

# and extract the split-read alignments
echo $(date +"[%b %d %H:%M:%S] Extract split-read alignments")
samtools view -h $INBAM1 \
    | ${LUMPYDIR}/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > ${INBAM1%.bam}.splitters.unsorted.bam

samtools view -h $INBAM2 \
    | ${LUMPYDIR}/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > ${INBAM2%.bam}.splitters.unsorted.bam

samtools view -h $INBAM3 \
    | ${LUMPYDIR}/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > ${INBAM3%.bam}.splitters.unsorted.bam
echo $(date +"[%b %d %H:%M:%S] Done") 

# sort both alignments for each sample
echo $(date +"[%b %d %H:%M:%S] Sort discordant and split alignments")
samtools sort ${INBAM1%.bam}.discordants.unsorted.bam -o ${INBAM1%.bam}.discordants
samtools sort ${INBAM1%.bam}.splitters.unsorted.bam -o ${INBAM1%.bam}.splitters

samtools sort ${INBAM2%.bam}.discordants.unsorted.bam -o ${INBAM2%.bam}.discordants
samtools sort ${INBAM2%.bam}.splitters.unsorted.bam -o ${INBAM2%.bam}.splitters

samtools sort ${INBAM3%.bam}.discordants.unsorted.bam -o ${INBAM3%.bam}.discordants
samtools sort ${INBAM3%.bam}.splitters.unsorted.bam -o ${INBAM3%.bam}.splitters
echo $(date +"[%b %d %H:%M:%S] Done")

### SV discovery/calling phase ###
echo $(date +"[%b %d %H:%M:%S] Start calling and genotyping")
$SPEEDSEQEXE sv \
-o $PREFIX \
-x ${FASTDIR}/references/genomes/wgEncodeDacMapabilityConsensusExcludable.bed \
-g \
-B $INBAM1,$INBAM2,$INBAM3 \
-D ${INBAM1%.bam}.discordants,${INBAM2%.bam}.discordants,${INBAM3%.bam}.discordants \
-S ${INBAM1%.bam}.splitters,${INBAM2%.bam}.splitters,${INBAM3%.bam}.splitters \
-R ${FASTDIR}/references/genomes/ucsc.hg19.fasta \
-t 8 \
-v \
2> ${FASTDIR}/slurmLOG/$PREFIX.lumpylog
echo $(date +"[%b %d %H:%M:%S] Done")
echo $(date +"[%b %d %H:%M:%S] All done!")
