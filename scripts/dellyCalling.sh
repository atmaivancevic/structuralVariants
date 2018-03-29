#!/bin/bash

## Script for calling SVs with delly
## Date: 28 March 2018 
##
## Example usage:
## INDIR=/data/neurogenetics/alignments/Illumina/genomes/allGenomes_chr10 sbatch --array 0-125 dellyCalling.sh

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

# load modules
module load BCFtools/1.3.1-GCC-5.3.0-binutils-2.25

# run the thing

### SV discovery/calling phase ###
echo $(date +"[%b %d %H:%M:%S] Starting delly")
echo "Processing file: "${QUERIES[$SLURM_ARRAY_TASK_ID]}

$DELLYEXE call \
-g ${GENOMEDIR}/ucsc.hg19.fasta \
-o ${OUTDIR}/SVcalling/dellyOut/${QUERIES[$SLURM_ARRAY_TASK_ID]}.bcf \
-x ${GENOMEDIR}/ucsc.hg19.excl.tsv \
${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]} \
2> ${FASTDIR}/slurmLOG/dellyCalling/${QUERIES[$SLURM_ARRAY_TASK_ID]}.log

echo $(date +"[%b %d %H:%M:%S] All done!")
