# Mis-polarized sites analysis

This folder generates the file `myhped_ms6.csv` that is used in
the mutational signature fitting analysis
[here](../../analyses/MutationalSignatureFitting).
There renamed as `3merSpectrum.David.ContainsRedCalls.20230920.csv`.

## Required files

We need to download the files for the ancestral human genome
from
[here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2)
and de-compress it.

We also need the WholeGenome Sequence from the UCSC version of
`hg38` that can be downloaded from 
[iGenomes](http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz)
.

Store the path of the decompressed directory in a `.env` file follwing
the example in `.env.sample`

We also need to copy the mutation data in the folder
`raw.data.nobackup`

```bash
mkdir -p raw.data.nobackup
mkdir -p raw.data.nobackup/igvreports
mkdir -p raw.data.nobackup/unisamplevcf
```

Use the folder `unisamplevcf` for the generated
unisample vcfs in the
[igv reports folder](https://github.com/harrispopgen/MUTYH_germline_effects/tree/main/data/IgvReports)
and the `igvreports` folder for the
csv files extracted from the manually scored
sheets.

The name of the files needs to be like (below) and must
match the name of the sheets tab.

```
{offspring}_{parent1}_{parent2}_{scorer}_MINIMAL.csv
```  

## Required software

Python requirements can be found in `requirements.txt` and can be
installed through pip.

```bash
python3 -m pip install -r requirements.txt
```

This project uses GNU sed version 4.9, note that 
the scripts in
this folder are not compatible with the default
installation of macOS.

## Usage

To generate the output files, run:

```bash
source .env
make all
```
