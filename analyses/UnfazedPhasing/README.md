## Description

This Bash script is designed to automate the execution of the `unfazed` tool for phasing genetic variants in a set of pedigree and BED files. `unfazed` is a command-line utility used for phasing genetic variants, and this script streamlines the process for multiple input files. 

## Prerequisites

Before running this script, ensure you have the following:

- **Python**: Python 3.9 is required. 

- **Python Packages**:
  - `numpy` version 1.23.1
  - `pysam` version 0.19.0
  - `cyvcf2` version 0.30.16
  - [unfazed GitHub Repository](https://github.com/jbelyeu/unfazed) version 08052022

- **Input Files**:
  - BED Files: These files should be located in the directory specified by `BED_file_dir` and are used for variant phasing.
  - Reference VCF File: The reference VCF file is specified by `reference_VCF` and is used for genetic variant processing.
  - Pedigree Files: Pedigree files should be located in the directory specified by `pedigree_file_dir`.
  - BAM Files: BAM files are specified by `BAM_file_directory` and are used for genetic variant processing.

## Usage

1. Modify the script's variables:
   - `BED_file_dir`: Path to the directory containing BED files.
   - `reference_VCF`: Path to the reference VCF file.
   - `pedigree_file_dir`: Path to the directory containing pedigree files.
   - `BAM_file_directory`: Path to the directory containing BAM files.
   - `output_dir`: Directory where the output files will be stored.

2. Define the array of BED/pedigree file pairs in the `file_pairs` variable. Each entry should be in the format "pedigree_file:BED_file".

3. Run the script using the following command:
   ```bash
   bash script_name.sh

