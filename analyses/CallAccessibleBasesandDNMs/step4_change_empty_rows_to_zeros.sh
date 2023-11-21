#!/bin/bash
#$ -N add_zeroes
#$ -o /net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/errors
#$ -e /net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/errors
#$ -m bea

module load modules modules-init

OUTPUT_DIR=/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/output

# change empty rows to zeros in final output file
for file in "$OUTPUT_DIR"/*
do
 sed -i 's/\t$/\t0/g' "$file"
done
