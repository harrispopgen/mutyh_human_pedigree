#!/bin/bash

module load modules modules-init
module load htslib/1.9-20 bcftools/1.12 bedtools/2.29.2

# Set the directory paths
bed_mask_file="/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/redo_filtered_vcfs/output/averaged_accessible_bases_no_cent_no_telo_no_segdups_noSNPs_mappable.bed"
masked_bed_dir="/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/redo_hap_IBD_intersected_filtered_vcfs/masked_beds"
negative_mask_dir="/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/surrogate_accessible_regions/stringent_dense_cluster_beds"

mkdir -p "${masked_bed_dir}/stringent_final_denoms"

declare -A individual_negative_masks
individual_negative_masks["C21"]="C22asfather_C21askid_vs_C23asfather_C21askid_overlap.bed"
individual_negative_masks["C22"]="C21asfather_C22askid_vs_C23asfather_C22askid_overlap.bed"
individual_negative_masks["C23"]="C21asfather_C23askid_vs_C22asfather_C23askid_overlap.bed"
individual_negative_masks["P1"]="P2P3asparents_P1askid_vs_P2P4asparents_P1askid_vs_P3P4asparents_P1askid_overlap.bed"
individual_negative_masks["P2"]="P1P3asparents_P2askid_vs_P1P4asparents_P2askid_vs_P3P4asparents_P2askid_overlap.bed"
individual_negative_masks["P3"]="P1P4asparents_P3askid_vs_P1P2asparents_P3askid_vs_P2P4asparents_P3askid_overlap.bed"
individual_negative_masks["P4"]="P1P2asparents_P4askid_vs_P1P3asparents_P4askid_vs_P2P3asparents_P4askid_overlap.bed"

declare -A individual_bed_files
individual_bed_files["C21"]="surrogate_denominator_highqual_SNPs_family2_C22asfather_C21askid.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_family2_C23asfather_C21askid.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed"
individual_bed_files["C22"]="surrogate_denominator_highqual_SNPs_family2_C21asfather_C22askid.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_family2_C23asfather_C22askid.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed"
individual_bed_files["C23"]="surrogate_denominator_highqual_SNPs_family2_C21asfather_C23askid.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_family2_C22asfather_C23askid.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed"
individual_bed_files["P1"]="surrogate_denominator_highqual_SNPs_P1_kid_P2_P3_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_P1_kid_P2_P4_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_P1_kid_P3_P4_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed"
individual_bed_files["P2"]="surrogate_denominator_highqual_SNPs_P2_kid_P1_P3_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_P2_kid_P1_P4_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_P2_kid_P3_P4_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed"
individual_bed_files["P3"]="surrogate_denominator_highqual_SNPs_P3_kid_P1_P2_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_P3_kid_P1_P4_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_P3_kid_P2_P4_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed"
individual_bed_files["P4"]="surrogate_denominator_highqual_SNPs_P4_kid_P1_P2_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_P4_kid_P1_P3_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed surrogate_denominator_highqual_SNPs_P4_kid_P2_P3_parents.PossibleDenovoAnnotation.accessible_accessible_no_cent_no_telo_no_segdups_noSNPs_mappable.bed"

# ===================== BED File Output Creation =====================

# Process BED files for each individual
for individual in "${!individual_negative_masks[@]}"; do
  temp_merged_bed=$(mktemp)
  final_masked_bed_file="${masked_bed_dir}/stringent_final_denoms/final_${individual}_denominator.bed"

  # Merge the BED files for this individual
  IFS=' ' read -ra bed_files <<< "${individual_bed_files[$individual]}"
  for bed_file in "${bed_files[@]}"; do
    cat "$masked_bed_dir/$bed_file" >> "$temp_merged_bed"
  done

  sorted_temp=$(mktemp)
  bedtools sort -i "$temp_merged_bed" > "$sorted_temp" || { echo "Sorting failed for $individual"; exit 1; }

  merged_temp=$(mktemp)
  bedtools merge -i "$sorted_temp" > "$merged_temp" || { echo "Merge failed for $individual"; exit 1; }

  bp_before=$(awk -F'\t' 'BEGIN {sum=0} {sum += $3 - $2} END {print sum}' "$merged_temp")

  bedtools intersect -v -a "$merged_temp" -b "$negative_mask_dir/${individual_negative_masks[$individual]}" > "$final_masked_bed_file" || { echo "Intersection failed for $individual"; exit 1; }

  bp_after=$(awk -F'\t' 'BEGIN {sum=0} {sum += $3 - $2} END {print sum}' "$final_masked_bed_file")

  echo "Processed: $individual"
  echo "Output: $final_masked_bed_file"
  echo "Count before: $bp_before"
  echo "Count after: $bp_after"

  rm -f "$temp_merged_bed" "$sorted_temp" "$merged_temp"

done

echo "Finished creating BED files!"

