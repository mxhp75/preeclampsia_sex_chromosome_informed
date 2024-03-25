#!/bin/bash -l

# job requirements
#SBATCH --partition=general
#SBATCH --ntasks=1
# set --cpus-per-task to the number of threads we request.
#SBATCH --cpus-per-task=4
#SBATCH --time=5:00:00
#SBATCH --mem=128GB

# Notification configuration
#SBATCH --job-name=preeclampsia_sex_informed_reduced_signatureMatrix
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melanie.smith@flinders.edu.au


## Melanie Smith ##
## 20240325 ##

## Script for running CIBERSORTx through apptainer -> docker ##
## To be run on DeepThought ##

##--------------------------------------------------------------------------------------------##
## Load modules
##--------------------------------------------------------------------------------------------##

module load shared
module load apptainer/1.1.9

##--------------------------------------------------------------------------------------------##
## Set environment variables
##--------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------##
## Set environment variables
##--------------------------------------------------------------------------------------------##

WORK_DIR="/scratch/user/smit1924/preeclampsia_sex_chromosome_informed/cibersortx_reduced_sigMatrix"
cd $WORK_DIR

##--------------------------------------------------------------------------------------------##
## Specify paths to data files and outputs
##--------------------------------------------------------------------------------------------##

OUT_DIR="${WORK_DIR}/outdir"

##--------------------------------------------------------------------------------------------##
## Run CIBERSORTx through apptainer
##--------------------------------------------------------------------------------------------##

apptainer run \
        -B ./data:/src/data \
        -B ${OUT_DIR}:/src/outdir \
        --pwd /src \
        fractions_latest.sif \
        --username melanie.smith@flinders.edu.au \
        --token e7d3e60a0c68bb01cbc06d543606066c \
        --verbose TRUE \
        --perm 100 \
        --fraction 0 \
        --single_cell TRUE \
        --refsample GSE182381_reduced_reference_sample.txt \
        --mixture bulkCountMatrix_all_gene_symbol.txt \
        --rmbatchSmode TRUE
