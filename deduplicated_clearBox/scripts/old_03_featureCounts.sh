#!/bin/bash -l

## Melanie Smith ##
## 20230608 ##

## Script for assigning genomic features to `sortByCoordOut.bam` files using subread/featureCounts ##
## To be run on DeepThought ##

# nb: conda environment activated in the run script

##--------------------------------------------------------------------------------------------##
## set project root directory
##--------------------------------------------------------------------------------------------##

PROJROOT=/scratch/user/smit1924/20230608_male_grch38

##--------------------------------------------------------------------------------------------##
## set number of threads
##--------------------------------------------------------------------------------------------##

cores=16

##--------------------------------------------------------------------------------------------##
## set project sub-directories
##--------------------------------------------------------------------------------------------##

alignOutput=${PROJROOT}/aligned_data   # the aligned data is here
readCounts=${PROJROOT}/readCounts # put the count files here

##--------------------------------------------------------------------------------------------##
## set reference file
##--------------------------------------------------------------------------------------------##

GRCh38_gtf=/scratch/user/smit1924/refSeq/ref_annotations/gencode.v29.annotation.gtf

##--------------------------------------------------------------------------------------------##
## Check the project root directory exists, else exit with error
##--------------------------------------------------------------------------------------------##

if [[ ! -d ${PROJROOT} ]]
  then
    echo -e "${PROJROOT} not found.\nExiting with Error code 1\n"
    exit 1
fi
echo -e "Found ${PROJROOT}\n"

##--------------------------------------------------------------------------------------------##
## Check the featureCounts directory exists, else make it
##--------------------------------------------------------------------------------------------##

if [[ ! -d ${readCounts} ]]
  then
	mkdir -p ${PROJROOT}/readCounts
fi
echo -e "Found ${readCounts}\n"

##--------------------------------------------------------------------------------------------##
## Assign genomic features using subread/featureCounts
##--------------------------------------------------------------------------------------------##

## Feature Counts - obtaining all sorted bam files
SAMPLES=$(find ${alignOutput} -name "*Aligned.sortedByCoord.out.bam" | tr '\n' ' ')

## Running featureCounts on the *sorted* bam files
featureCounts -Q 10 -p --countReadPairs -g gene_id --fracOverlap 1 -T ${cores} -a ${GRCh38_gtf} -o ${readCounts}/counts.out ${SAMPLES}

## Storing the output in a single file
cut -f1,7- ${readCounts}/counts.out | \
   sed 1d > ${readCounts}/CDNOJANXX_male_readCounts.txt




