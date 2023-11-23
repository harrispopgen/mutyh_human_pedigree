#!/bin/bash

##### This script modifies output files to ensure that empty values (output from samtools depth) are replaced with zeros. #####

module load modules modules-init # Initialize modules environment

# Define the directory containing the output files to be processed.
OUTPUT_DIR=/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/output

# change empty rows to zeros in final output file with sed command
for file in "$OUTPUT_DIR"/*
do
 sed -i 's/\t$/\t0/g' "$file"
done
