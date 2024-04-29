#!/bin/bash -l

# job requirements
#SBATCH --partition=general
#SBATCH --ntasks=1
# set --cpus-per-task to the number of threads we request.
#SBATCH --cpus-per-task=16
#SBATCH --time=6:00:00
#SBATCH --mem=16GB

# Notification configuration
#SBATCH --job-name=20240429_deduplicated_featureCounts
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melanie.smith@flinders.edu.au


## Melanie Smith ##
## 20240417 ##

## Script for running STAR alignment ##
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
## Set project root directory
##--------------------------------------------------------------------------------------------##

PROJROOT="/scratch/user/smit1924/preeclampsia_sex_chromosome_informed"
SCRIPTS="${PROJROOT}/deduplicated_clearBox/scripts"
##--------------------------------------------------------------------------------------------##
## Run samtools sort
##--------------------------------------------------------------------------------------------##

conda activate featureCounts

bash "${SCRIPTS}/03_featureCounts.sh"

conda deactivate
