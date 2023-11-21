#! /bin/bash
#$ -l h_rt=150:00:00,mfree=10G
#$ -o /net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/unfiltered_vcfs/errors
#$ -e /net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/unfiltered_vcfs/errors
#$ -m bea
#$ -P sage

##### GATK variant annotator: PossibleDeNovo of all individuals in extended pedigree (run on vcf with no prior masks applied) #####

module load modules modules-init modules-gs # initialize modules
module load GATK/4.2.6.1

vcf=/net/harris/vol1/sharing/MUTYH_pedigree/extraFiles/harris_grc_wgs_1.HF.final.vcf.gz
outdir=/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/PossibleDeNovo_files/unfiltered_vcfs/output

# List of ped files
ped_files=(
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

for pedfile in "${ped_files[@]}"; do
  # Generate the output file name
  filename=$(basename "$pedfile" .ped)
  GATK_outfile="${filename}.HF.PossibleDenovoAnnotation.vcf.gz"

  # Run GATK VariantAnnotator
  gatk VariantAnnotator \
    -A PossibleDeNovo \
    -V "$vcf" \
    --pedigree "$pedfile" \
    -O "$outdir/$GATK_outfile"

  ##### bcftools: filtering out high & low quality de novo SNPs from GATK annotated vcf file #####

  final_vcf="$outdir/$GATK_outfile"
  final_outfile="highqual_SNPs_${filename}.PossibleDenovoAnnotation.vcf.gz"

  zcat "$final_vcf" | grep -v "\./\." | bcftools view -m2 -M2 -v snps -i '%FILTER="PASS"' | bcftools filter -e 'INFO/hiConfDeNovo="." & INFO/loConfDeNovo="."' -Oz -o "$outdir/SNVs_only/$final_outfile"
done

