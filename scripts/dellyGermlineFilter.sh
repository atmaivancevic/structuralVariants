#!/bin/bash

## Script for filtering merged SVs with the delly germline filter
## Requires at least 20 unrelated samples
## Date: 24 April 2018 
##
## Example usage:
## INDIR=/fast/users/a1211880/outputs/SVcalling/dellyOut sbatch dellyGermlineFilter.sh

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

# set variables
DELLYEXE=$FASTDIR/executables/delly-0.7.8/delly_v0.7.8_parallel_linux_x86_64bit

# run the thing
echo $(date +"[%b %d %H:%M:%S] Go to dir")
cd $INDIR
pwd

echo $(date +"[%b %d %H:%M:%S] Check for merged bcf")
ls merged.bcf

echo $(date +"[%b %d %H:%M:%S] Apply delly germline filter")
$DELLYEXE filter -f germline -o germline.bcf merged.bcf

echo $(date +"[%b %d %H:%M:%S] All done!")
