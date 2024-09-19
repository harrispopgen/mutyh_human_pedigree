
require(reticulate)
# works if you use /opt/anaconda3/bin/pip install SigProfilerExtractor 
# and then use /opt/anaconda3/bin/python
#use_python("/usr/local/bin/python3")
#py_config() # NOTE: had to go to Rstudio >> Preferences >> Python and select /usr/bin/python3 for 
#library(SigProfilerAssignmentR) #updated 9/21/23 # package ‘SigProfilerAssignmentR’ is not available for this version of R
library(SigProfilerExtractorR) #updated 9/21/23
todaysdate=format(Sys.Date(),"%Y%m%d")
inputdir="/Users/annabelbeichman/Documents/UW/Human_MUTYH/results/signatureFitting/"
# uses default cosmic version 3.3
# you wrote out data ready for SPE/SPA in previous script
############## run sig profiler extractor #########
spe_outdir=paste0(inputdir,"sigprofilerextractor/",todaysdate,"_Results_DSCorrection_COSMICv3.3/")
dir.create(spe_outdir)

# run per family:
sigprofilerextractor(input_type = "table", 
                     paste0(spe_outdir,"perFamily/"), 
                     input_data=paste0(inputdir,"AllSpectra.FormattedForSigProfilerExtractor.PERFAMILY.txt"), 
                     reference_genome="GRCh38",
                     minimum_signatures=1,
                     maximum_signatures=3,
                     nmf_replicates=100) # 
# above 3 got error "There is no signature over the thresh-hold stability. We are selecting the lowest possible number of signatures."
## note: max sigs must be < total samples
# There is no signature over the thresh-hold stability. We are selecting the lowest possible number of signatures.
# Error: IndexError: list index out of range

# run per individual:
sigprofilerextractor(input_type = "table", 
                     paste0(spe_outdir,"perIndividual/"), 
                     input_data=paste0(inputdir,"AllSpectra.FormattedForSigProfilerExtractor.PERINDIVIDUAL.txt"), 
                     reference_genome="GRCh38",
                     minimum_signatures=1,
                     maximum_signatures=10,
                     nmf_replicates=100) # 

### make sure that rows are in same order as example data and that signatures look good (they do.)