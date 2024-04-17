#!/bin/bash

# Script to run Picard Tools Remove Duplicates
# Bams must be sorted
# Melanie Smith
# 20240412

# Set project root directory
PROJECT_ROOT_DIRECTORY=/scratch/user/smit1924/20230608_male_grch38
ALIGN_OUTPUT="${PROJECT_ROOT_DIRECTORY}/aligned_data"   # Aligned data directory
DEDUP_OUTPUT="${PROJECT_ROOT_DIRECTORY}/deduplicated_data" # Deduplicated bams directory
SAMPLE_BAM="${ALIGN_OUTPUT}/SCP4733_small.bam"

# Function to create directory if it doesn't exist
create_directory() {
  local directory=$1
  if [[ ! -d "${directory}" ]]; then
    echo "Directory ${directory} not found. Creating new directory."
    mkdir -p "${directory}" || { echo "Error: Failed to create directory ${directory}."; exit 1; }
  fi
}

# Check if output data directory exists, else create it
create_directory "${DEDUP_OUTPUT}"

# Check if aligned data directory exists, else exit with error
if [[ ! -d "${ALIGN_OUTPUT}" ]]; then
  echo "Error: Aligned data directory ${ALIGN_OUTPUT} not found. Exiting."
  exit 1
fi

echo "Found aligned data directory: ${ALIGN_OUTPUT}"

# Remove Duplicates
picard MarkDuplicates \
        -I "${SAMPLE_BAM}" \
        -O "${DEDUP_OUTPUT}/SCP4733_small_marked_duplicates.bam" \
        -M "${DEDUP_OUTPUT}/SCP4733_small_mark_dup_metrics.txt" \
        -REMOVE_DUPLICATES true
