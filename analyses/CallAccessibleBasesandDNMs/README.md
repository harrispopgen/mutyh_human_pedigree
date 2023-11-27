# Calling putative de novo mutations and accessible bases per pedigree individual

Code used for calling DNMs with GATK PossibleDeNovo, filtering out mutations by various parameters, and using samtools depth to determine the number of accessible sites per each individual's genome. The surrogate DNM calling technique is also used in this code.

## Table of Contents
1. [Step 1: GATK Possible De Novo Annotation & Filtering](#step-1-gatk-possible-de-novo-annotation--filtering)
2. [Step 2: Parallel Coverage Computation](#step-2-coverage-computation)
3. [Step 3: Add Zeros to Empty Rows](#step-3-add-zeros-to-rows)
4. [Step 4: Merge Files, Calculate and Filter Average Depth Scores](#step-4-merge-files-and-calculate-average-depths)
5. [Step 5: Average Depths across Chromosomes across Individuals](#step-5-average-depths-across-chromosomes-across-individuals)
6. [Step 6: Mask BED and VCF files](#step-6-mask-bed-and-vcf-files)
7.
8.
9.
10.
11.

<hr>

## Step 1: GATK Possible De Novo Annotation & Filtering
`01_GATK_PossibleDeNovo_Annotation_Filtering.sh`

This script annotates and filters possible de novo mutations (DNMs) across multiple individuals in an extended human family pedigree. It utilizes the GATK VariantAnnotator for the PossibleDeNovo annotation and bcftools for filtering.

#### Primary Functions:
1. **Annotation**: Applies the GATK VariantAnnotator to a multi-sample VCF file to annotate possible de novo mutations.
2. **Filtering**: Filters the annotated VCF to identify high-quality de novo SNPs with bcftools, removing low-quality and non-SNP variants.

#### Necessary Inputs:
1. **Multi-sample VCF file**: A VCF file containing variant data for multiple individuals (e.g., in our case, `harris_grc_wgs_1.HF.final.vcf.gz`).
2. **PED files**: A series of PED files defining family structures and individuals to consider for DNM calling. These files include various family trio configurations (including all surrogate trio combinations).

#### Required Packages:
- GATK (Genome Analysis Toolkit) version 4.2.6.1
- bcftools for variant filtering

#### Output:
The script generates VCF files for each PED file analyzed, with two sets of output files:
1. **Annotated VCFs**: Files containing PossibleDeNovo annotations (e.g., `family1.HF.PossibleDenovoAnnotation.vcf.gz`). Optional to keep these intermediate files in the final output. 
2. **Filtered High-Quality SNPs VCFs**: Files containing only high-quality SNPs annotated as possible de novo mutations (e.g., `highqual_SNPs_family1.PossibleDenovoAnnotation.vcf.gz`). These files are the ones used in downstream analysis!

Each output VCF is named based on the corresponding PED file and is stored in the designated output directory.

<hr>

## Step 2: Coverage Computation
`02_ParallelCoverageComputation_15Individuals23Chromosomes.sh`

This script computes coverage scores for 15 individuals across 23 chromosomes, processing each combination in parallel. It divides the genome into 10K base pair (bp) chunks and calculates coverage using the samtools depth command.

#### Primary Functions:
1. **Parallel Computation**: Utilizes Grid Engine job arrays to run 345 jobs in parallel, covering all individual-chromosome combinations.
2. **Coverage Calculation**: Computes coverage in 10K bp chunks for each individual and chromosome, outputting the data to individual files.

#### Necessary Inputs:
1. **BAM Files**: BAM files for each individual.

#### Required Packages:
- `samtools` version 1.9

#### Process Details:
1. **Initialization**: Sets up the modules environment and loads required versions of samtools.
2. **Arrays Definition**: Defines arrays of individuals and chromosomes to be processed.
3. **Coverage Calculation**: For each job, selects an individual and chromosome based on job array ID, calculates coverage in 10K bp chunks, and outputs the data to a text file. We use a minimum mapping quality score of 20 for the depth calculation.

#### Outputs:
- Coverage txt files for each individual and chromosome combination, named in the format: `individualID_chromosome.txt` (e.g., `C11.1079256_chr1.txt`).
- Each file contains the chromosome, start and end positions of each chunk, and the calculated coverage scores across each 10K bp chunk.

<hr>

### Step 3: Add Zeros to Rows
`03_ChangeEmptyRowstoZeros.sh`

This script is designed to process output files from the coverage calculation step, specifically modifying them to replace any empty values (which are an output artifact from the samtools depth command) with zeros.

#### Primary Function:
- **Data Cleaning**: Modifies output files by replacing empty values in the data with zeros, which is necessary for downstream analyses. 

#### Necessary Inputs:
- **Output Files Directory**: The script targets txt files in the output directory of `Step 2`, which contain the coverage data.

#### Required Packages:
- NA

#### Outputs:
- The script modifies the files in-place; therefore, the output is the updated txt files in the same directory with empty values replaced by zeros. These can be named appropriately when modifying this code.

<hr>

### Step 4: Merge Files and Calculate Average Depths
`04_MergeDepthFiles_and_CalculateAvgDepths.py`

This script is designed to merge depth files across different samples and chromosomes, and then calculate the average depth for each 10KB window (averaged across the 15 individuals in our dataset). We also exclude inaccessible regions from the final output by filtering out averaged depth scores < 12 and >120 (as calculated by samtools depth in step 2).  

#### Primary Functions:
- **merge_data**: Merges data across samples and chromosomes.
- **log_error**: Logs errors to a specified file.

#### Necessary Inputs:
- Directory path to directory containing txt files of samtools depth output per chromosome per individual (with zeroes added from step 3)
- Directory path for error logging.

#### Required Packages:
- os
- pandas
- random
- numpy

#### Output:
- Merged depth files per chromosome across all individuals
- Average depth calculations per 10KB window for each chromosome, **filtered by accessibility**
- Log files for errors

This script consolidates and averages depth information across multiple individuals, and then filters out inaccessible regions (depth < 12 or depth >120). The sampling and validation steps included in the script are optional but are useful for verifying the integrity of the merged data.

#### Optional Sampling and Validation:
- **Sampling Percentage**: Set to a small percentage (e.g., 0.001) for validation purposes.
- This process involves randomly sampling a set number of sites and comparing the averaged depth values in the output file with those in the original individual depth files to ensure the data from the merged files matches the original input. This step confirms the validity of the merging and averaging processes.

<hr>

### Step 5: Average Depths across Chromosomes across Individuals
`05_Calculate_Depth_Summary_Statistics.py`

This script calculates the average (mean) depth across 10KB chunks per chromosome, already averaged across individuals from the previous code. While it also computes various summary statistics, these were mainly optional for our purposes. The primary goal is to generate BED files with and without the average depth scores.

#### Primary Functions:
- **Calculate Mean Depth**: Computes the mean depth for each chromosome.
- **Generate Summary Statistics**: Optionally calculates median, minimum, maximum depths, and 2.5 times the mean depth.
- **Count # of Accessible Bases**: Accumulates the number of accessible bases across chromosomes (filtered by average coverage in step 4).
- **Output Summary Table**: Generates a summary table of statistics.
- **Convert BED Files**: Processes and outputs BED files, with and without depth scores.

#### Necessary Inputs:
- Directory path to output of step 4 (averaged depth files per chromosome across individuals)

#### Required Packages:
- os
- pandas

#### Output:
- A summary table in CSV format containing the mean, median, range, and 2.5 times the mean depth for each chromosome.
- Total count of the average number of accessible bases per individual (previously filtered for accessibility in step 4!).
- BED files with and without the average depth scores.

The script processes depth data to create output BED files of the average depth scores across all chromosomal 10kb chunks and individuals. The optional summary statistics offer additional depth insights.

<hr>

### Step 6: Mask BED and VCF files
`06_Grab_AccessibleSites_fromVCFs_andBEDs.sh`

This script filters original VCF files for regions that are mappable and have an average accessible depth (calculated in steps 4 & 5). It excludes regions like telomeres, centromeres, segmental duplications, and common SNPs.

#### Primary Functions:
- Filter VCF files for uniquely mappable and accessible regions.
- Exclude specific genomic regions (telomeres, centromeres, etc.) from analysis.
- Calculate the number of sites removed by each genomic mask (important for checking which masks are significantly affecting your the number of bases in the final callsets)

#### Necessary Inputs:
- PossibleDeNovo annotated VCF files (from Step 1).
- BED files for accessible regions, mappable regions, telomeres, centromeres, segmental duplications, and common SNPs.

#### Required Packages:
- `htslib/1.9-20`
- `bcftools/1.12`
- `bedtools/2.29.2`

#### Output:
- Masked VCF files based on the specified genomic regions.
- Masked BED file of accessible bases after applying the same masks.
- Count of variants and sites removed at each step, recorded in an output file.

#### Download Links for Masking Files:
- Telomere regions ("Gap Locations"): [https://genome.ucsc.edu/cgi-bin/hgTables](#)
- Centromere regions: [https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=centromeres](#)
- Mappable regions (UMAP S100): [https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1667695194_a5i1edTDwMtMljwsHOLN6iODnGrM&c=chrX&g=hub_213889_Umap_100](#)
- Segmental duplication regions: [https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_doSchema=describe+table+schema](#)
- Common SNP regions (dbSNP 153): [https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=varRep&hgta_track=dbSnp153Composite&hgta_table=dbSnp153Common&hgta_doSchema=describe+table+schema](#)

*NOTE: accessible regions are calculated by averaging the depth scores of 10kb chromosomal chunks across all individuals in the sequencing set up. This BED file should be generated in step 5!*

<hr>

### Step 7:

<hr>

### Step 8:

<hr>

### Step 9:

<hr>

### Step 10:

<hr>

### Step 11:

