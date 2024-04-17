#!/bin/bash

# Script for assigning genomic features to `sortByCoordOut.bam` files using subread/featureCounts
# To be run on DeepThought

# Melanie Smith
# 20241704

set -euo pipefail

##--------------------------------------------------------------------------------------------##
## Set project root directory
##--------------------------------------------------------------------------------------------##

PROJROOT="/scratch/user/smit1924/preeclampsia_sex_chromosome_informed"

##--------------------------------------------------------------------------------------------##
## Set number of threads
##--------------------------------------------------------------------------------------------##

cores=16

##--------------------------------------------------------------------------------------------##
## Set project sub-directories
##--------------------------------------------------------------------------------------------##

ALIGN_OUTPUT="${PROJROOT}/deduplicated_clearBox/deduplicated_bams"   # The deduplicated data is here
READ_COUNTS="${PROJROOT}/deduplicated_clearBox/readCounts"      # Put the count files here

##--------------------------------------------------------------------------------------------##
## Set reference file
##--------------------------------------------------------------------------------------------##

GRCh38_GTF="/scratch/user/smit1924/refSeq/ref_annotations/gencode.v29.annotation.gtf"

##--------------------------------------------------------------------------------------------##
## Check the project root directory exists, else exit with error
##--------------------------------------------------------------------------------------------##

if [[ ! -d "${PROJROOT}" ]]; then
  echo "${PROJROOT} not found. Exiting with Error code 1." >&2
  exit 1
fi

echo "Found ${PROJROOT}"

##--------------------------------------------------------------------------------------------##
## Check the featureCounts directory exists, else make it
##--------------------------------------------------------------------------------------------##

if [[ ! -d "${READ_COUNTS}" ]]; then

mkdir -p "${READ_COUNTS}"
echo "Made ${READ_COUNTS}"

fi

##--------------------------------------------------------------------------------------------##
## Assign genomic features using subread/featureCounts
##--------------------------------------------------------------------------------------------##

# Feature Counts - obtaining all sorted bam files
SAMPLES=$(find "${ALIGN_OUTPUT}" -name "*_marked_duplicates.bam" | tr '\n' ' ')

# Running featureCounts on the *deduplicated* bam files
featureCounts -Q 10 -p --countReadPairs -g gene_id --fracOverlap 1 -T "${cores}" -a "${GRCh38_GTF}" -o "${READ_COUNTS}/counts.out" ${SAMPLES}

# Storing the output in a single file
cut -f1,7- "${READ_COUNTS}/counts.out" | sed 1d > "${READ_COUNTS}/deduplicated_readCounts.txt"

