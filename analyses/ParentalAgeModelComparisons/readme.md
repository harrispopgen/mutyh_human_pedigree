# Comparing observed mutation counts to expectations under a parental age model

These scripts were used to compare observed DNM mutation counts and spectra (from our study and Sherwood et al. 2023) to those expected under a parental age model based on the Poisson regressions in J—nsson et al. (2017).

### Step 1: 
`step_1_processUnfazedResults.R`
This is a little helper script to process files from the Unfazed pipeline into input format needed for Step 2. 
It takes as input: 
1. `total_counts_with_parents_IGV_filtered.txt`, a file which contains per-individual information about total DNMs called, total phased, parental ages, and total accessible bases. Note that these counts reflect sites that have been through the "minimal" VAF filter and went through IGV inspection. 

2. `phased_mutations_per_individual/[ID]_mutations_IGV_filtered.txt`, a set of files (one per individual) containing the coordinates of phased mutations, which parental haplotype they were phased to, and the 1-mer mutation type. Note that the parent generation wasn't phased, so there are no files for individuals P1-P4.

This script outputs:
1. `ALLCHILDREN.TotalCounts.PhasedToEachParent.WithPhasingSucessRate.txt`, a summary of mutation counts and phasing success rates in a format that is easier for Step 2 to work with 

2. `ALLCHILDREN.todate.SpectrumCounts.PhasedToEachParent.long.txt`, a file containing the *phased* 1-mer spectrum for each individual (excluding P1-P4 which weren't phased), in 'long' format


### Step 2: 
`step_2_CalculateExpectedSpectraBasedonParentalAge.JonssonTableS9.minimalVAF.bugfix.IGVFiltered.tidied.R`

This is a big multi-step script that carries out many analyses on data from this study and Sherwood et al.'s. It generates all the figures related to the Jonsson regression, Poisson probabilities, heatmaps, and minimum effect size estimates. It is best run in RStudio bit by bit to keep track of what it's doing and look at plots as they are generated.

It takes many files as input:


1. `Jonsson.Table9.ForMutyhProject.Wide.ForR.txt`: a wide version of Jonsson et al.'s Table S9 showing their per-mutation-type Poisson regression's maternal and paternal slopes and intercepts.

2. `parental_ages.txt`: a table listing the parental ages at birth for each individual in the pedigree

3. `mut_type_counts_noX_IGV_filtered.txt`: a table of *unphased* 1-mer mutation spectrum counts for each individual in 'long' format. These counts went through minimal VAF filtering and IGV inspection.

4. `total_counts_with_parents_IGV_filtered.txt`, a file which contains per-individual information about total DNMs called, total phased, parental ages, and total accessible bases. Note that these counts reflect sites that have been through the "minimal" VAF filter and went through IGV inspection. In the context of this script, we use this table to get information about the accessible genome size for each individual.

5. `ALLCHILDREN.TotalCounts.PhasedToEachParent.WithPhasingSucessRate.txt`:  a summary of mutation counts and phasing success rates in a format that is easier for Step 2 to work with 

6. `ALLCHILDREN.todate.SpectrumCounts.PhasedToEachParent.long.txt`: a file generated in Step 1, containing the *phased* 1-mer spectrum for each individual (excluding P1-P4 which weren't phased), in 'long' format

7. `SherwoodMutyhPaper2023/TableOfDNMsandPhasing.txt`: a file containing mutation count and phasing success information from Sherwood et al. 2023 (from their Table 1) for comparison with our paper. Their reported fractions were converted to mutation counts by multiplying them by the reported DNM SNV load. 

8. `SherwoodMutyhPaper2023/Sherwood.SummedUpSpectrum.TableS2.forR.txt`: a file containing the mutation spectra from Sherwood et al. (reported in their Table S2) that were summed up per group (individuals grouped by parent carrier status, e.g. mutyh biallelic)

The script outputs a great many text files and figures into the `scriptoutput/` directory. Final figures for the manuscript are generated in Step 3, below.


### Step 3:
`step_3_makeManuscriptReadyFigures.R`

This script takes in the many dataframe text files generated in Step 2 (`scriptoutput/`) and replots main text and SI figure panels.

