# Mutational signature analysis
Code used for formatting 3-mer spectra for sigprofilerextractor (or sigfit), and carrying out signature extraction in sigprofilerextractor

### Step 1: 
`formatForSigfitAndSPE.R`
This is a helper script that converts a 3-mer mutation spectrum into formats that can be read by SigprofilerExtractor (or Sigfit; they use different formats)

It generated spectra per-individual, per-family and per-group

The script outputs:

1. `AllSpectra.FormattedForSigProfilerExtractor.PERINDIVIDUAL.txt`: a file containing the *per-individual* 3-mer spectra formatted for sigprofilerextractor.

2. `AllSpectra.FormattedForSigProfilerExtractor.PERFAMILY.txt`: *per-family* 3-mer spectra (summed up across all children of a family)

### Step 2:

Run sigprofilerextractor using R wrapper (SigProfilerExtractorR) based on spectra that are per-individual or summed per-family

The script takes in the the AllSpectra.*.txt files generated above. 