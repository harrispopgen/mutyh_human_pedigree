# README

This folder contains the code to generate [igv reports](https://github.com/igvteam/igv-reports) out of the
trio sequencing de novo variants.

## Required coniguration

An `.env` file is required containing the path to
the required files:

- Read data stored in bams (`BAMPATH`), folder containing bam files.
- Final callset in bed files (`BEDPATH_MINIMAL`), containing the folder with bed files that correspond to the final DNM call.
- A pedigree file (`PEDPATH`) that includes the relationship among samples in each pedigree. 
- The main VCF resulted from the HaplotypeCaller step (`MAINVCF`).

See the file `.env.sample` for an example.

The pipeline will also need an extra file containing the chromosome
sizes of the given assembly (GRCh38 here). To do so I have used a
faidx file that corresponds to the indexed version of [the fasta file available here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/).

## Software install

To run the pipeline we need to install and activate a conda
environtment for Snakemake. To do so follow the guide
[here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

The run the main steps and the IGV reports, the pipeline uses
modules. The cluster should have thus the following installed
modules.

* module load python/3.7.7
* module load IGV/2.9.4
* module load igv-reports/1.7.0
* module load GATK/4.2.6.1

Although in our execution the following packages were also loaded 
modules, the pipeline will need them just to be available in `$PATH`.

* module load python/3.9.13
* module load bedtools/2.29.2
* module load samtools/1.17
* module load bcftools/1.17
* module load tabix/0.2.6

## Usage

Once the files above are generated and
referenced into the .env file and the software is
installed
we can run the pipeline from the Make file:

```bash
source .env
make all
```
