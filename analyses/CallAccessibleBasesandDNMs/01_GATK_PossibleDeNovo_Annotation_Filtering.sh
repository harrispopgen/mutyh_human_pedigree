#! /bin/bash

##### This script runs the GATK variant annotator for the PossibleDeNovo annotation on all individuals in an extended pedigree. #####
##### It is applied to a multi-sample VCF without prior masks. The script also filters high and low-quality de novo SNPs using bcftools. #####

module load modules modules-init modules-gs # initialize modules environment
module load GATK/4.2.6.1 # Load GATK version 4.2.6.1

# Define input VCF and output directory
vcf=/net/harris/vol1/sharing/MUTYH_pedigree/extraFiles/harris_grc_wgs_1.HF.final.vcf.gz
outdir=/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/unfiltered_vcfs/output

# # Define the list of pedigree files
ped_files=(
	# Array of file paths to PED files for different family configurations (all surrogate trio combinations considered)
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/trio1.ped"
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/trio3.ped"
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/onechildperped/family2_C21asfather_C22askid.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/onechildperped/family2_C22asfather_C21askid.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/onechildperped/family2_C23asfather_C21askid.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/onechildperped/family4_C41asfather_C42askid.ped"
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/onechildperped/family2_C21asfather_C23askid.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/onechildperped/family2_C22asfather_C23askid.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/onechildperped/family2_C23asfather_C22askid.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/onechildperped/family4_C42asfather_C41askid.ped"
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P1_kid_P2_P3_parents.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P1_kid_P3_P4_parents.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P2_kid_P1_P4_parents.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P3_kid_P1_P2_parents.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P3_kid_P2_P4_parents.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P4_kid_P1_P3_parents.ped"
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P1_kid_P2_P4_parents.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P2_kid_P1_P3_parents.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P2_kid_P3_P4_parents.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P3_kid_P1_P4_parents.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P4_kid_P1_P2_parents.ped"  
	"/net/harris/vol1/data/human_mutyh_pedigree/pedFiles/proband_combos/P4_kid_P2_P3_parents.ped"  
)

# Loop through each pedigree file
for pedfile in "${ped_files[@]}"; do
  # Generate the output file name based on the PED file name
  filename=$(basename "$pedfile" .ped)
  GATK_outfile="${filename}.HF.PossibleDenovoAnnotation.vcf.gz"

  # Run GATK VariantAnnotator
  gatk VariantAnnotator \
    -A PossibleDeNovo \
    -V "$vcf" \
    --pedigree "$pedfile" \
    -O "$outdir/$GATK_outfile"

  ##### Filter out high & low-quality de novo SNPs from GATK annotated vcf file using bcftools #####
  final_vcf="$outdir/$GATK_outfile"
  final_outfile="highqual_SNPs_${filename}.PossibleDenovoAnnotation.vcf.gz"
  
  # Filtering process
  zcat "$final_vcf" | grep -v "\./\." | bcftools view -m2 -M2 -v snps -i '%FILTER="PASS"' | bcftools filter -e 'INFO/hiConfDeNovo="." & INFO/loConfDeNovo="."' -Oz -o "$outdir/SNVs_only/$final_outfile"
done

