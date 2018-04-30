#!/bin/bash

## Script for cleaning up manta output
## Use this script after mantaSoloCall.sh
## Date: 30 April 2018 
##
## Example usage:
## INDIR=/fast/users/a1211880/outputs/SVcalling/mantaOut SAMPLE=FD00825684 sbatch mantaCleanup.sh

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=0:10:00
#SBATCH --mem=1GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Run the thing

echo $(date +"[%b %d %H:%M:%S] Go to dir")
cd $INDIR/$SAMPLE/results/variants
pwd

# unzip variant output file
echo $(date +"[%b %d %H:%M:%S] Unzip")
gunzip diploidSV.vcf.gz 

# rename output to include sample name
echo $(date +"[%b %d %H:%M:%S] Rename and move vcf")
cp diploidSV.vcf "$SAMPLE"_diploidSV.vcf
mv "$SAMPLE"_diploidSV.vcf $INDIR

echo $(date +"[%b %d %H:%M:%S] All done!")
