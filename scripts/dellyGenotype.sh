#!/bin/bash

## Script for genotyping SVs with delly
## Date: 24 April 2018 
##
## Example usage:
## INDIR=/data/neurogenetics/alignments/Illumina/genomes/allGenomes SITELIST=/fast/users/a1211880/outputs/SVcalling/dellyOut/sites.bcf sbatch --array 0-125 dellyGenotype.sh

#SBATCH -A robinson
#SBATCH --qos=long
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=7-00:00
#SBATCH --mem=8GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# define key variables
export FASTDIR=/fast/users/$USER

DELLYEXE=$FASTDIR/executables/delly-0.7.8/delly_v0.7.8_parallel_linux_x86_64bit
OUTDIR=$FASTDIR/outputs
GENOMEDIR=$FASTDIR/references/genomes

# define query bam files
QUERIES=($(ls $INDIR/*.bam | xargs -n 1 basename))

# run the thing

### SV discovery/calling phase ###
echo $(date +"[%b %d %H:%M:%S] Starting delly genotyping")
echo "Processing file: "${QUERIES[$SLURM_ARRAY_TASK_ID]}

$DELLYEXE call \
-g ${GENOMEDIR}/ucsc.hg19.fasta \
-v ${SITELIST} \
-o ${OUTDIR}/SVcalling/dellyOut/${QUERIES[$SLURM_ARRAY_TASK_ID]}.geno.bcf \
-x ${GENOMEDIR}/ucsc.hg19.excl.tsv \
${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]} \
2> ${FASTDIR}/slurmLOG/dellyCalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.log

echo $(date +"[%b %d %H:%M:%S] All done!")
