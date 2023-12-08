#! /bin/bash

# Load required modules
module load modules modules-init modules-gs
module load python/3.9.13 numpy/1.23.1 pysam/0.19.0 cyvcf2/0.30.16 unfazed/08052022

# Define variables
BED_file_dir="/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/minimal_filtered_mutations/beds"
reference_VCF="/net/harris/vol1/sharing/MUTYH_pedigree/extraFiles/harris_grc_wgs_1.HF.final.vcf.gz"
pedigree_file_dir="/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/peds_for_unfazed"
BAM_file_directory="/net/harris/vol1/sharing/MUTYH_pedigree/bamFiles/rename"
output_dir="/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/minimal_filtered_mutations/phasing_output"

# Define an indexed array of BED/pedigree file pairs
file_pairs=(
    "trio1.ped:C11.bed"
    "trio1.ped:C12.bed"
    "family2_C22asfather_C21askid.ped:C21.bed"
    "family2_C23asfather_C21askid.ped:C21.bed"
    "family2_C21asfather_C22askid.ped:C22.bed"
    "family2_C23asfather_C22askid.ped:C22.bed"
    "family2_C21asfather_C23askid.ped:C23.bed"
    "family2_C22asfather_C23askid.ped:C23.bed"
    "trio3.ped:C31.bed"
    "trio3.ped:C32.bed"
    "family4_C42asfather_C41askid.ped:C41.bed"
    "family4_C41asfather_C42askid.ped:C42.bed"
)

# Loop through the array and run unfazed
for pair in "${file_pairs[@]}"; do
    IFS=":" read -ra parts <<< "$pair"
    ped_file=${parts[0]}
    bed_file=${parts[1]}

    # Remove the .ped and .bed file extensions
    ped_file_base=${ped_file%.ped}
    bed_file_base=${bed_file%.bed}

    # Generate the output text file name without extensions
    output_txt="${output_dir}/${ped_file_base}_${bed_file_base}_output.txt"

    # Run unfazed
    unfazed -d "${BED_file_dir}/${bed_file}" \
            -s $reference_VCF \
            -p "${pedigree_file_dir}/${ped_file}" \
            -g 38 \
            -b $BAM_file_directory \
            > $output_txt
done

