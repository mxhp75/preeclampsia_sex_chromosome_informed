#!/bin/bash -l

# job requirements
#SBATCH --partition=general
#SBATCH --ntasks=1
# set --cpus-per-task to the number of threads we request.
#SBATCH --cpus-per-task=16
#SBATCH --time=6:00:00
#SBATCH --mem=16GB

# Notification configuration
#SBATCH --job-name=test_indexSingleBam
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melanie.smith@flinders.edu.au


## Melanie Smith ##
## 20230629 ##

## Script for running indexing of a single bam file ##
## To be run on DeepThought ##

##--------------------------------------------------------------------------------------------##
## Load modules
##--------------------------------------------------------------------------------------------##

module load shared
module load Miniconda3/4.9.2

##--------------------------------------------------------------------------------------------##
## Inside the 'featureCounts' conda environment
##--------------------------------------------------------------------------------------------##

# mamba
# subread
# multiqc
# samtools

##--------------------------------------------------------------------------------------------##
## Run samtools sort
##--------------------------------------------------------------------------------------------##

conda activate featureCounts

samtools index -b /scratch/user/smit1924/20230608_male_grch38/aligned_data/SCP_3954_T_S22_GRCh38_Aligned.sortedByCoord.out.bam

conda deactivate
