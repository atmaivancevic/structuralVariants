#!/bin/bash

## Script for calling SVs with lumpy
## Use this script if you have non-trio samples (singletons or large families)
## Date: 24 April 2018 
##
## Example usage:
## INBAM1=/data/neurogenetics/alignments/Illumina/genomes/allGenomes_chr10/NA12878_chr10.bam PREFIX=NA12878_chr10 sbatch lumpySoloCalling.sh

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
# need to extract the discordant paired-end alignment
echo $(date +"[%b %d %H:%M:%S] Extract discordant paired-end alignments")
samtools view -b -F 1294 $INBAM1 > ${INBAM1%.bam}.discordants.unsorted.bam
echo $(date +"[%b %d %H:%M:%S] Done")

# and extract the split-read alignment
echo $(date +"[%b %d %H:%M:%S] Extract split-read alignments")
samtools view -h $INBAM1 \
    | ${LUMPYDIR}/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > ${INBAM1%.bam}.splitters.unsorted.bam
echo $(date +"[%b %d %H:%M:%S] Done") 

# sort both alignments for each sample
echo $(date +"[%b %d %H:%M:%S] Sort discordant and split alignments")
samtools sort ${INBAM1%.bam}.discordants.unsorted.bam -o ${INBAM1%.bam}.discordants
samtools sort ${INBAM1%.bam}.splitters.unsorted.bam -o ${INBAM1%.bam}.splitters
echo $(date +"[%b %d %H:%M:%S] Done")

### SV discovery/calling phase ###
echo $(date +"[%b %d %H:%M:%S] Starting calling and genotyping")
$SPEEDSEQEXE sv \
-o $PREFIX \
-x ${FASTDIR}/references/genomes/wgEncodeDacMapabilityConsensusExcludable.bed \
-g \
-B $INBAM1 \
-D ${INBAM1%.bam}.discordants \
-S ${INBAM1%.bam}.splitters \
-R ${FASTDIR}/references/genomes/ucsc.hg19.fasta \
-t 8 \
-v \
2> ${FASTDIR}/slurmLOG/$PREFIX.lumpylog
echo $(date +"[%b %d %H:%M:%S] Done")
echo $(date +"[%b %d %H:%M:%S] All done!")
