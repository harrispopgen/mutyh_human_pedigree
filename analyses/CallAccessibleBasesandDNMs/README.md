# Calling putative de novo mutations and accessible bases per pedigree individual

Code used for calling DNMs with GATK PossibleDeNovo, filtering out mutations by various parameters, and using samtools depth to determine the number of accessible sites per each individual's genome. The surrogate DNM calling technique is also used in this code.

## Table of Contents
1. [Step 1: GATK Possible De Novo Annotation & Filtering](#step-1-gatk-possible-de-novo-annotation--filtering)
2. [Step 2: Parallel Coverage Computation](#step-2-coverage-computation)
3. [Step 3: Add Zeros to Empty Rows](#step-3-add-zeros-to-rows)
4. [Step 4: Merge Files, Calculate and Filter Average Depth Scores](#step-4-merge-files-and-calculate-average-depths)
5. [Step 5: Average Depths across Chromosomes across Individuals](#step-5-average-depths-across-chromosomes-across-individuals)
6. [Step 6: Mask BED and VCF files](#step-6-mask-bed-and-vcf-files)
7. [Step 7: Generate IBD Segments for Pedigree Individuals](#step-7-Identity-by-Descent-calling)
8. [Step 8: Identify Shared hap-IBD Tracts](#step-8-identify-shared-ibd-tracts)
9. [Step 9: Grab Shared IBD Regions from (surrogate) VCF and BED files](#step-9-grab-shared-ibd-regions-from-vcf-and-bed-files)
10. [Step 10: Filter Masked VCFs by Site and Sample Specific Metrics](#step-10-filter-masked-vcfs-by-site-and-sample-specific-metrics)
11. [Step 11: Apply sparsity filter, interactively visualize mutations, and generate final mutation counts callset](#step-11-apply-sparsity-filter-and-generate-final-mutation-counts-callset)
12. [Step 12: Generate final callset of accessible bases for surrogate individuals](#step-12-generate-final-denominators)
13. [Step 13: Calculate and plot number of accessible bases and rates per individual](#step-13-generate-final-counts-and-rates)

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
- `GATK (Genome Analysis Toolkit) version 4.2.6.1`
- `bcftools/1.12`

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

## Step 3: Add Zeros to Rows
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

## Step 4: Merge Files and Calculate Average Depths
`04_MergeDepthFiles_and_CalculateAvgDepths.py`

This script is designed to merge depth files across different samples and chromosomes, and then calculate the average depth for each 10KB window (averaged across the 15 individuals in our dataset). We also exclude inaccessible regions from the final output by filtering out averaged depth scores < 12 and >120 (as calculated by samtools depth in step 2).  

#### Primary Functions:
- **merge_data**: Merges data across samples and chromosomes.
- **log_error**: Logs errors to a specified file.

#### Necessary Inputs:
- Directory path to directory containing txt files of samtools depth output per chromosome per individual (with zeroes added from step 3)
- Directory path for error logging.

#### Required Packages:
- `os`
- `pandas`
- `random`
- `numpy`

#### Output:
- Merged depth files per chromosome across all individuals
- Average depth calculations per 10KB window for each chromosome, **filtered by accessibility**
- Log files for errors

This script consolidates and averages depth information across multiple individuals, and then filters out inaccessible regions (depth < 12 or depth >120). The sampling and validation steps included in the script are optional but are useful for verifying the integrity of the merged data.

#### Optional Sampling and Validation:
- **Sampling Percentage**: Set to a small percentage (e.g., 0.001) for validation purposes.
- This process involves randomly sampling a set number of sites and comparing the averaged depth values in the output file with those in the original individual depth files to ensure the data from the merged files matches the original input. This step confirms the validity of the merging and averaging processes.

<hr>

## Step 5: Average Depths across Chromosomes across Individuals
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
- `os`
- `pandas`

#### Output:
- A summary table in CSV format containing the mean, median, range, and 2.5 times the mean depth for each chromosome.
- Total count of the average number of accessible bases per individual (previously filtered for accessibility in step 4!).
- BED files with and without the average depth scores.

The script processes depth data to create output BED files of the average depth scores across all chromosomal 10kb chunks and individuals. The optional summary statistics offer additional depth insights.

<hr>

## Step 6: Mask BED and VCF files
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

## Step 7: Identity by Descent calling

#### Necessary Inputs:
- Phased 15-individual WGS VCF from Beagle
- 1000 Genomes Whole Genome Sequencing data
- Recombination maps

#### Required Packages:
- `htslib/1.12`
- `bcftools/1.12`
- `hapibd`

This nextflow pipeline merges the the 15 individual WGS dataset with the high-coverage 1000 Genomes dataset, keeps common variants that are found at a minor allele frequency of >10%, and removes all remaining sites that are un-phased or contain missing data (see processes `merge_with_onekg`, `filter_vcf`, `remove_unphased_and_missing_site` for details)
performs IBD calling using hapibd (process `infer_ibd`)
converts the hapibd output format (see table 1) into a custom format (table 2) for downstream processing outlined in Step __ using the Python script `reformat_ibd_output.py` (see process `convert_dataformat`)

<b>Table 1:</b>
|:----------:|
| sample 1 id |
| haplotype index of sample 1 (1 or 2) |
| sample 2 id |
| haplotype index of sample 2 (1 or 2) |
| chromosome |
| start position of IBD segment |
| end position of the IBD segment |

<b>Table 2:</b>
|:----------:|
| chromosome |
| start position of IBD segment |
| end position of IBD segment |
| sample_haplotype 1 |
| sample_haplotype 2 |

The sample haplotype in the custom format is obtained by first mapping the haplotype index (either 1 or 2) in the hapibd output file to `a` or `b` respectively, and then concatenating it to the hapibd sample id. For example, for a given IBD segment in hapibd, if sample 1 id = S11, and haplotype index = 1, the resulting sample_haplotype in the custom format will be S11a.
Relevant path definitions and resource requirements for individual processes are found in the `nextflow.config` file, although note that the cluster settings are specific to the sun grid engine (SGE) and all files are located on a local hard drive, therefore the pipeline is not expected to be reproducible out of the box.

<hr>

## Step 8: Identify Shared IBD Tracts
`08_GetGrandparent_hapIBD_segs.py`

This script is designed to analyze genomic segments shared among the first generation of individuals in our pedigree. It identifies segments inherited from the grandparents based on Identity by Descent (IBD) data.

#### Primary Functions:
- Determine if two hap-IBD segments overlap
- Read and process centromere and telomere tract data
- Merge overlapping segments and remove small gaps between them (considering centromere and telomere regions)
- Read and process IBD data for identifying shared segments among first generation individuals
- Export data structures for duo and trio segment overlaps within the familial structures of our pedigree

#### Necessary Inputs:
- Centromere tract data: [https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=centromeres](#)
- Telomere tract data: [https://genome.ucsc.edu/cgi-bin/hgTables](#)
- IBD files for each chromosome (generated in step 8)

#### Required Packages:
- `os`
- `matplotlib` (configured for non-GUI environments)

#### Output:
The script generates a detailed mapping of shared hap-IBD segments among the first generation (grandparent) individuals of our pedigree, considering duo and trio combinations. The output is organized into dictionaries (shared, duo_segs, trio_segs) for downstream analysis.

<hr>

## Step 9: Grab Shared IBD Regions from VCF and BED files
`09_GrabSurrogate_hapIBDregions.sh`

This script processes VCF files in conjunction with corresponding BED files to apply haplotype-based IBD regions as postive masks. It uses `bedtools` and `bcftools` for applying these masks to the BED files and VCF files, respectively.

#### Primary Functions:
- Maps VCF files to their corresponding hap-IBD BED files, applying them as a positive mask.
- Counts sites removed in the VCF by each IBD mask using `bedtools`.
- Also masks the original accessible bases bed file (generated in step 5) to each hap-IBD bed file, corresponding to the relevant individual used in the surrogate method.

#### Necessary Inputs:
- VCF and BED files generated from step 6.
- BED file outputs generated from step 7 (containing shared IBD regions between different combinations of individuals).

#### Required Packages:
- `htslib/1.9-20`
- `bcftools/1.12`
- `bedtools/2.29.2`

#### Output:
- hap-IBD masked VCF files and BED files per each individual where the surrogate method was applied
- Counts of sites removed from the numerator (VCF) and denominator (modified accessible bases BED file) by the mask.

The script integrates VCF and BED files to apply specific IBD regions as a positive mask, necessary for calling true de novo mutations in individuals where siblings are used as surrogate parents.

<hr>

## Step 10: Filter Masked VCFs by Site and Sample Specific Metrics
`10_postmasking_filter_site_and_sample_specific_scores.Rmd`

This script takes the output VCFs generated from step 6 (individuals where the surrogate method was not applied) and step 9 (individuals where the surrogate method was applied. It filters VCFs (converted to CSVs for easier processing) by standard quality metrics, then converts the filtered output into renamed csvs for downstream processing in python. 

NOTE: to convert VCFs to CSV files, we used the tool pdbio, with the command: `pdbio vcf2csv --expand-info --expand-samples $sample.vcf.gz > $sample.csv`

#### Primary Functions:
- process data frames to clean up unecessary columns
- apply the following site specific filters:
	- QD > 2.0
	- FS < 60.0
	- MQRankSum > -12.5
        - ReadPosRankSum > -8.0
        - SOR < 3.0
- apply the following sample specific filters (for each individual, the child and parents all must pass these filters at each variant site)
	- 12 > DP < 120
	- GQ >= 20
- generate VAF column ("AD_frac") for downstream filtering
- optional: merge, count, and plot mutations for individuals where the surrogate method was not applied
- save filtered csv files

#### Necessary Inputs:
- csv files converted from the VCF file outputs of step 6 and step 9 

#### Required Packages:
- `dplyr`
- `tidyr`
- `ggplot2` (optional)

#### Output:
- filtered csv files which can be plugged immediately into step 11 for further processing

<hr>

## Step 11: Apply sparsity filter and generate final mutation counts callset
`11_TagSparseMutations_GenerateNegMasksforSurrogates.ipynb`

Various functions to filter out dense mutations in surrogate individuals using a sliding window approach; create and apply BED files where there are dense regions of mutations; visually inspect mutations unique and shared across different parents per surrogate individual; and grab and count "sparse mutations".

NOTE: we use the Jupyter notebook environment to run the Dash app for interactive visualizations. More information on configuring this specific environment can be found here: https://dash.plotly.com/dash-in-jupyter
 
#### Primary Functions:
- identify "dense" and "sparse" mutations per csv callset using a sliding window approach
- generate BED files of "dense" clusters to be used as negative masks on the final denominator for surrogate individuals
- visually inspect dense and sparse mutations with an interactive graph
- visually inspect mutations shared across different callsets for the same individual (for example, C21 with C22 as father vs. C21 with C23 as father)
- grab all "shared" and "unique" mutations across callsets per individual, count and plot
- visually inspect VAF distribution of mutations, apply filter (0.3 < VAF > 0.7)
- create a BED file containing the positions of each individual's mutations

#### Necessary Inputs:
- filtered csv files generated from step 10

#### Required Packages:
- `os`
- `numpy`
- `pandas`
- `zipfile`
- `re`
- `collections`
- `pybedtools`
- `matplotlib`
- `dash`
- `jupyter_dash` (if using Jupyter Notebook GUI)
- `plotly`

#### Output:
- html files of Dash app visualizations, zipped
- BED files of "dense" clusters per surrogate individual callset
- final tally of mutations per individual in csv format
- BED files of all identified mutations per individual
- optional: plot of VAF distribution of mutations per individual

<hr>

## Step 12: Generate Final Denominators
`12_GenerateFinalBEDs_forSurrogates.sh`

Apply BED files generated from step 11 to mask indivdual denominator regions (generated in step 9) where mutations cannot be called in surrogate individuals (AKA "dense" clusters). 

#### Primary Functions:
- Merge BED files generated from step 9 for each individual.
- Sort combined BED files.
- Intersect the combined BED files with negative mask files generated in step 11 (per individual) to finalize denominator BED files.
- Calculates base pair counts before and after the intersection.
 
#### Necessary Inputs:
- Masked BED files generated from step 9.
- negative mask files generated from step 11.
- An associative array of individual-specific negative masks.
- An associative array of BED files for each individual.

#### Required Packages:
- `htslib/1.9-20`
- `bcftools/1.12`
- `bedtools/2.29.2`

#### Output:
- Final masked BED files for each individual stored (which are used as the "denominators" per each surrogate individual)
- Console output displaying the base pair counts before and after processing.

<hr>

## Step 13: Generate Final Counts and Rates
`13_generate_final_counts_rates.sh`

Calculate and vizualize final mutation counts and rates generated from the final VCF and BED files for all individuals. 

#### Primary Functions:
-  **Counting Accessible Bases**: Iterates through BED files (generated in step 12 for surrogate individuals, and step 5 for the non-surrogate individuals) & calculate the total number of accessible bases.
-  **Data Visualization**: Create bar plots to visually represent mutation counts, denominators, and normalized mutation rates across all individuals of the pedigree.

#### Necessary Inputs:
- Directory path containing BED files (`final_denoms`), generated from step 12 for surrogate individuals and step 5 for non-surrogate individuals (ie the accessible bases BED file). 

#### Required Packages:
- `os`
- `pandas`
- `matplotlib.pyplot`
- `matplotlib.ticker.ScalarFormatter`

#### Output:
- **Console Output**: Displays the number of accessible bases processed from each file.
- **Bar Plots**:
  1. Mutation Counts by Data Frame.
  2. Denominators multiplied by 2 for each data frame.
  3. Normalized mutation counts (mutation counts divided by denominators times 2), ie, the number of DNMs per base pair per generation for each pedigree individual. 
