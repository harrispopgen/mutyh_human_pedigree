#!/bin/bash
#$ -o /net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/redo_hap_IBD_intersected_filtered_vcfs/errors
#$ -e /net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/redo_hap_IBD_intersected_filtered_vcfs/errors

module load modules modules-init
module load htslib/1.9-20 bcftools/1.12 bedtools/2.29.2

# Set the directory paths
vcf_dir="/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/redo_filtered_vcfs/output"
bed_dir="/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/surrogate_accessible_regions/hap_IBD_beds"
output_dir="/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/redo_hap_IBD_intersected_filtered_vcfs/output"
bed_mask_file="/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/redo_filtered_vcfs/output/averaged_accessible_bases_no_cent_no_telo_no_segdups_noSNPs_mappable.bed"
masked_bed_dir="/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/redo_hap_IBD_intersected_filtered_vcfs/masked_beds"

count_sites_removed_masked() {
    local masked_bed_file="$1"

    # Count the number of sites removed by each mask
    local masked_sites_removed=$(bedtools intersect -v -a "$bed_mask_file" -b "$masked_bed_file" | wc -l)

    # Output the counts for this masked BED file
    echo "For $(basename "$masked_bed_file"):"
    echo "Number of sites removed from denominator by mask: $masked_sites_removed"
}

# Define an associative array to map VCF files to their respective bed files
declare -A file_mapping
file_mapping["highqual_SNPs_family2_C21asfather_C22askid.PossibleDenovoAnnotation.accessible.vcf.gz"]="surrogate_parent_positive_mask_C21_C22.bed"
file_mapping["highqual_SNPs_P2_kid_P1_P3_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P1_P2_double_surrogates.bed P2_P3_double_surrogates.bed P2_v_surrogate_parentsP1_P3.bed"
file_mapping["highqual_SNPs_family2_C21asfather_C23askid.PossibleDenovoAnnotation.accessible.vcf.gz"]="surrogate_parent_positive_mask_C21_C23.bed"
file_mapping["highqual_SNPs_P2_kid_P1_P4_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P1_P2_double_surrogates.bed P2_P4_double_surrogates.bed P2_v_surrogate_parentsP1_P4.bed"
file_mapping["highqual_SNPs_family2_C22asfather_C21askid.PossibleDenovoAnnotation.accessible.vcf.gz"]="surrogate_parent_positive_mask_C21_C22.bed"
file_mapping["highqual_SNPs_P2_kid_P3_P4_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P2_P3_double_surrogates.bed P2_P4_double_surrogates.bed P2_v_surrogate_parentsP3_P4.bed"
file_mapping["highqual_SNPs_family2_C22asfather_C23askid.PossibleDenovoAnnotation.accessible.vcf.gz"]="surrogate_parent_positive_mask_C22_C23.bed"
file_mapping["highqual_SNPs_P3_kid_P1_P2_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P1_P3_double_surrogates.bed P2_P3_double_surrogates.bed P3_v_surrogate_parentsP1_P2.bed"
file_mapping["highqual_SNPs_family2_C23asfather_C21askid.PossibleDenovoAnnotation.accessible.vcf.gz"]="surrogate_parent_positive_mask_C21_C23.bed"
file_mapping["highqual_SNPs_P3_kid_P1_P4_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P1_P3_double_surrogates.bed P3_P4_double_surrogates.bed P3_v_surrogate_parentsP1_P4.bed"
file_mapping["highqual_SNPs_family2_C23asfather_C22askid.PossibleDenovoAnnotation.accessible.vcf.gz"]="surrogate_parent_positive_mask_C22_C23.bed"
file_mapping["highqual_SNPs_P3_kid_P2_P4_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P2_P3_double_surrogates.bed P3_P4_double_surrogates.bed P3_v_surrogate_parentsP2_P4.bed"
file_mapping["highqual_SNPs_family4_C41asfather_C42askid.PossibleDenovoAnnotation.accessible.vcf.gz"]="surrogate_parent_positive_mask_C41_C42.bed"
file_mapping["highqual_SNPs_P4_kid_P1_P2_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P1_P4_double_surrogates.bed P2_P4_double_surrogates.bed P4_v_surrogate_parentsP1_P2.bed"
file_mapping["highqual_SNPs_family4_C42asfather_C41askid.PossibleDenovoAnnotation.accessible.vcf.gz"]="surrogate_parent_positive_mask_C41_C42.bed"
file_mapping["highqual_SNPs_P4_kid_P1_P3_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P1_P4_double_surrogates.bed P3_P4_double_surrogates.bed P4_v_surrogate_parentsP1_P3.bed"
file_mapping["highqual_SNPs_P1_kid_P2_P3_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P1_P2_double_surrogates.bed P1_P3_double_surrogates.bed P1_v_surrogate_parentsP2_P3.bed"
file_mapping["highqual_SNPs_P4_kid_P2_P3_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P2_P4_double_surrogates.bed P3_P4_double_surrogates.bed P4_v_surrogate_parentsP2_P3.bed"
file_mapping["highqual_SNPs_P1_kid_P2_P4_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P1_P2_double_surrogates.bed P1_P4_double_surrogates.bed P1_v_surrogate_parentsP2_P4.bed"
file_mapping["highqual_SNPs_P1_kid_P3_P4_parents.PossibleDenovoAnnotation.accessible.vcf.gz"]="P1_P3_double_surrogates.bed P1_P4_double_surrogates.bed P1_v_surrogate_parentsP3_P4.bed"


# Loop through the VCF files
for vcf_file in "${!file_mapping[@]}"; do
  bed_files=(${file_mapping[$vcf_file]})

  base_name=$(basename "$vcf_file" .vcf.gz)

  vcf_input="$vcf_dir/$vcf_file"
  decompressed_vcf_input="$vcf_dir/${base_name}.vcf"
  output_file="$output_dir/${base_name}_hapIBD_masked.vcf.gz"

  tabix -p vcf "$vcf_input"
  gunzip -c "$vcf_input" > "$decompressed_vcf_input"

  # Create a temporary concatenated BED file for the current VCF file
  temp_bed_cat_file=$(mktemp)
  for bed in "${bed_files[@]}"; do
    cat "$bed_dir/$bed" >> "$temp_bed_cat_file"
  done

  # Create a masked BED file by applying the bed mask for the current VCF file
  masked_bed_file="$masked_bed_dir/surrogate_denominator_${base_name}_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed"
  bedtools intersect -u -a "$bed_mask_file" -b "$temp_bed_cat_file" > "$masked_bed_file"
  count_sites_removed_masked "$masked_bed_file"
  
  bcftools sort -Ov "$decompressed_vcf_input" -o "$decompressed_vcf_input.sorted"

  # Apply the masked bed file as a positive mask for the current VCF file
  bcftools view -T "$masked_bed_file" "$decompressed_vcf_input.sorted" | \
    bcftools sort -Oz -o "$output_file"

  tabix -p vcf "$output_file"

  rm "$temp_bed_cat_file"

  echo "Processed: $vcf_file"
  echo "Output: $output_file"
done


