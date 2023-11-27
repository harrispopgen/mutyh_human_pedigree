# Calling putative de novo mutations and accessible bases per pedigree individual

Code used for calling DNMs with GATK PossibleDeNovo, filtering out mutations by various parameters, and using samtools depth to determine the number of accessible sites per each individual's genome. The surrogate DNM calling technique is also used in this code.

### Step 1: 
`01_GATK_PossibleDeNovo_Annotation_Filtering.sh`

This script annotates and filters possible de novo mutations (DNMs) across multiple individuals in an extended human family pedigree. It utilizes the GATK VariantAnnotator for the PossibleDeNovo annotation and bcftools for filtering.

#### Primary Functions:
1. **Annotation**: Applies the GATK VariantAnnotator to a multi-sample VCF file to annotate possible de novo mutations.
2. **Filtering**: Filters the annotated VCF to identify high-quality de novo SNPs with bcftools, removing low-quality and non-SNP variants.

#### Necessary Inputs:
1. **Multi-sample VCF file**: A VCF file containing variant data for multiple individuals (e.g., in our case, `harris_grc_wgs_1.HF.final.vcf.gz`).
2. **PED files**: A series of PED files defining family structures and individuals to consider for DNM calling. These files include various family trio configurations (including all surrogate trio combinations).

#### Required Packages/Modules:
- GATK (Genome Analysis Toolkit) version 4.2.6.1
- bcftools for variant filtering

#### Output:
The script generates VCF files for each PED file analyzed, with two sets of output files:
1. **Annotated VCFs**: Files containing PossibleDeNovo annotations (e.g., `family1.HF.PossibleDenovoAnnotation.vcf.gz`). Optional to keep these intermediate files in the final output. 
2. **Filtered High-Quality SNPs VCFs**: Files containing only high-quality SNPs annotated as possible de novo mutations (e.g., `highqual_SNPs_family1.PossibleDenovoAnnotation.vcf.gz`). These files are the ones used in downstream analysis!

Each output VCF is named based on the corresponding PED file and is stored in the designated output directory.

### Step 2:

### Step 3:

### Step 4:

### Step 5:

### Step 6:

### Step 7:

### Step 8:

### Step 9:

### Step 10:

### Step 11:

