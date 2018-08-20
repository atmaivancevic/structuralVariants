#!/bin/bash  

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time=10:00:00
#SBATCH --mem=32GB

# INDIR=/data/neurogenetics/alignments/Illumina/genomes/allGenomes OUTDIR=/fast/users/a1211880/outputs/SVcalling/cnvnatorOut sbatch --array 0-159 cnvnatorCall.sh

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# directories 
CHROMDIR=/data/neurogenetics/RefSeq/chromosomes/

# define query bam files
QUERIES=($(ls $INDIR/*.bam | xargs -n 1 basename))

# load modules
module load ROOT/v5.34.09-foss-2015b
module load GCC/4.9.3-binutils-2.25
module load CNVnator/0.3.2-GCC-4.9.3-binutils-2.25

# run the thing

### CNV discovery/calling phase ###
echo $(date +"[%b %d %H:%M:%S] Starting cnvnator")
echo "Processing file: "${QUERIES[$SLURM_ARRAY_TASK_ID]}

echo "Extracting read mapping from BAM files..."

cnvnator -root ${QUERIES[$SLURM_ARRAY_TASK_ID]}.root \
-tree ${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]}

echo "Generating a histogram..."

cnvnator -root ${QUERIES[$SLURM_ARRAY_TASK_ID]}.root \
-his 100 -d $CHROMDIR

echo "Calculating stats..."

cnvnator -root ${QUERIES[$SLURM_ARRAY_TASK_ID]}.root \
-stat 100

echo "Patitioning..."

cnvnator -root ${QUERIES[$SLURM_ARRAY_TASK_ID]}.root \
-partition 100

echo "Calling CNVs..."

cnvnator -root ${QUERIES[$SLURM_ARRAY_TASK_ID]}.root \
-call 100 \
-f ${INDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]} \
> ${OUTDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]}.out \
2> ${OUTDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]}.log

echo "Converting to VCF format..."

cnvnator2VCF.pl ${OUTDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]}.out $CHROMDIR \
> ${OUTDIR}/${QUERIES[$SLURM_ARRAY_TASK_ID]}.vcf

echo $(date +"[%b %d %H:%M:%S] All done!")
