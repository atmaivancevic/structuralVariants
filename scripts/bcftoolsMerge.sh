#!/bin/bash

## Script for merging genotyped SVs with bcftools
## Date: 24 April 2018 
##
## Example usage:
## INDIR=/fast/users/a1211880/outputs/SVcalling/dellyOut sbatch bcftoolsMerge.sh

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

# load modules
module load BCFtools/1.3.1-GCC-5.3.0-binutils-2.25

# run the thing
echo $(date +"[%b %d %H:%M:%S] Go to dir")
cd $INDIR
pwd

echo $(date +"[%b %d %H:%M:%S] Set bcf list")
genotypedList=$(find *.geno.bcf)
echo $genotypedList

echo $(date +"[%b %d %H:%M:%S] Merge all genotyped bcfs")
bcftools merge -m id -O b -o merged.bcf $genotypedList

echo $(date +"[%b %d %H:%M:%S] All done!")
