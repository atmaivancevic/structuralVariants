#!/bin/bash

## Script for merging SVs with delly
## Date: 24 April 2018 
##
## Example usage:
## INDIR=/fast/users/a1211880/outputs/SVcalling/dellyOut sbatch dellyMerge.sh

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

# define key variables
export FASTDIR=/fast/users/$USER

DELLYEXE=$FASTDIR/executables/delly-0.7.8/delly_v0.7.8_parallel_linux_x86_64bit
OUTDIR=$FASTDIR/outputs
GENOMEDIR=$FASTDIR/references/genomes

# run the thing
echo $(date +"[%b %d %H:%M:%S] Got to dir")
cd $INDIR
pwd

echo $(date +"[%b %d %H:%M:%S] Set bcf list")
bcfList=$(find *.bcf)
echo $bcfList

echo $(date +"[%b %d %H:%M:%S] Merge all bcfs")
$DELLYEXE merge -o sites.bcf $bcfList

echo $(date +"[%b %d %H:%M:%S] All done!")
