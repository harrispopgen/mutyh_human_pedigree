require(ggplot2)
require(ggrepel)
require(tidyr)
require(dplyr)
require(reshape2)
### stylistic choices:
# replace 'proband' with 'parent'
# call them DNMs not SNVs (unlike Sherwood)
# no hyphenation of mono- and bi-allelic
# MAKE SURE YOU'RE ALWAYS USING EXPECTED COUNTS THAT HAVE BEEN CORRECTED FOR ACC. GENOME SIZE (variables will end with something like *_TimesAccGenomeRatio)
########### Replotting figures for Manuscript ########### 
wd="/Users/annabelbeichman/Documents/UW/Human_MUTYH/results/jonssonRegressions/2024-09-09_Results_MinVAF_bugfix_IGVFiltered_MinEffectSizeAnalysis_DSCorrection_BUGFIX_FORREVISIONS/" # working directory -- 20240822 now has DS and shared allele correction; post bug-fix 20240909
indir=paste0(wd,"/script_output/") # get from Google drive

# main figure outdir:
outdir=paste0(wd,"manuscript_quality_figures/")
dir.create(outdir,showWarnings = F)

# for SI figures:
SIoutdir=paste0(outdir,"SIFigures/")
dir.create(SIoutdir,showWarnings = F)
########### Total mutation count vs Jonsson scatterplot ########### 
totalPerChild_EXPECTED_and_OBSERVED <- read.table(paste0(indir,"MutationCounts.totalPerChild_EXPECTED_and_OBSERVED.ContainsPoissonProbs.txt"),header=T,sep = "\t")

# have "Parent generation" be "monoallelic mother+father"
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$label=="Parent generation",]$label <- "monoallelic mother+father"

compareTotalCounts_ToJonsson <- ggplot(totalPerChild_EXPECTED_and_OBSERVED,aes(x=EXPECTED_totalPerChild_TimesAccGenomeRatio,y=Observed_totalPerChild,color=label))+
  geom_point(size=2)+
  geom_abline()+
  geom_text_repel(aes(label=child_id),show.legend = F)+
  theme_bw()+
  ylab("Observed total DNM burden")+
  xlab("Expected total DNM burden based on parental ages\n and corrected for accessible genome size")+
  scale_color_brewer(palette="Dark2")+
  theme(legend.title = element_blank())+
  theme(text = element_text(size=14))

compareTotalCounts_ToJonsson
ggsave(paste0(outdir,"ourData.compareTotalCounts_ToJonsson.pdf"),compareTotalCounts_ToJonsson,height=5,width=7)
ggsave(paste0(outdir,"ourData.compareTotalCounts_ToJonsson.png"),compareTotalCounts_ToJonsson,height=5,width=7)



########### SPECTRA ########### 


##### OBSERVED AND EXPECTED SPECTRA visualizations ##### 
# Faceted scatter plot (main text)
# column plot (counts and fractions in SI)
expectedAndObservedSpectra <- read.table(paste0(indir,"expectedAndObservedSpectra_PerIndividual.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),header=T,sep="\t")

# remove P2:
expectedAndObservedSpectra <- expectedAndObservedSpectra[expectedAndObservedSpectra$child_id!="P2",]

####### spectra: faceted scatter plot #######
expectedObservedSpectra_facetedScatterplot <- ggplot(expectedAndObservedSpectra,aes(x=Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio,y=Observed_Count,color=group))+
  geom_point()+
  facet_wrap(~gsub("\\.",">",Mut_Class_DontSeparateCpG),scales="free")+
  geom_abline()+
  geom_text_repel(aes(label=child_id),show.legend=F)+
  theme_bw()+
  ylab("Observed total DNM burden")+
  xlab("Expected total DNM burden based on parental ages\n and corrected for accessible genome size")+
  scale_color_brewer(palette="Dark2")+
  theme(legend.title = element_blank(),text=element_text(size=14))
expectedObservedSpectra_facetedScatterplot
ggsave(paste0(outdir,"Spectra.observed.expected.facetedscatterplot.pdf"),expectedObservedSpectra_facetedScatterplot,height=5,width=8)
ggsave(paste0(outdir,"Spectra.observed.expected.facetedscatterplot.png"),expectedObservedSpectra_facetedScatterplot,height=5,width=8)


##### spectra: heatmap (our data per family + sherwood + BXD mice) ######
# per family: expectedAndObservedSpectra_SummedOverFamilies_melt #
expectedAndObservedSpectra_SummedOverFamilies_melt_FORHEATMAP <- read.table(paste0(indir,"expectedAndObservedSpectra_SummedOverFamilies_melt.FORHEATMAP.txt"),header=T,sep="\t") # this variable: expectedAndObservedSpectra_SummedOverFamilies_melt

# sherwood: sherwoodSummedSpectrumObserved_PlusExpectations_melt; this variable: Observed_Over_Expected_Counts
sherwoodSummedSpectrumObserved_PlusExpectations_melt_FORHEATMAP <- read.table(paste0(indir,"sherwoodSummedSpectrumObserved_PlusExpectations_melt.FORHEATMAP.txt"),header=T,sep="\t")

# first need to do a bit of processing for BXD:
# these are counts of each mutation type from Kelley ;
# Add CpGs to C>Ts 


BXD_Counts <- read.table("/Users/annabelbeichman/Documents/UW/Human_MUTYH/information/BXD_counts/BXD.Counts.FromKelley.20231103.txt",header=T,sep="\t") # update this path as needed
head(BXD_Counts)


# add cpgs to c>t category: 
BXD_Counts$C.T_plusCpG.TpG <- BXD_Counts$C.T + BXD_Counts$CpG.TpG
# get rid of old categories:
BXD_Counts <- select(BXD_Counts,-c("C.T","CpG.TpG"))
# and rename to C.T_plusCpG.TpG to C.T 
BXD_Counts <- BXD_Counts %>%
  rename("C.T" = "C.T_plusCpG.TpG")

# melt it:
BXD_Counts_melt <- melt(BXD_Counts,id.vars = c("strain_type","label","total_generations_inbreeding"),variable.name = "Mut_Class_DontSeparateCpG",value.name = "mutation_count")

head(BXD_Counts_melt)

BXD_Counts_melt$mutation_count_DividedByGensInbreeding <- BXD_Counts_melt$mutation_count/BXD_Counts_melt$total_generations_inbreeding

## no longer making data compositional
# get fractions:
#BXD_Counts_melt <- BXD_Counts_melt %>%
#  group_by(strain_type,label) %>% 
#  mutate(mutation_fraction = mutation_count/sum(mutation_count)) %>%
#  ungroup()

# instead, want to divide D strain by B and BXD68 by B
# make it wide again:
BXD_Counts_melt_wide <- pivot_wider(BXD_Counts_melt,names_from = "strain_type",id_cols = "Mut_Class_DontSeparateCpG",values_from = "mutation_count_DividedByGensInbreeding")
head(BXD_Counts_melt_wide)

BXD_Counts_melt_wide$D_over_B <- BXD_Counts_melt_wide$D/BXD_Counts_melt_wide$B
BXD_Counts_melt_wide$BXD68_over_B <- BXD_Counts_melt_wide$BXD68/BXD_Counts_melt_wide$B
BXD_Counts_melt_wide$variable <- "mutation_count_DividedByGensInbreeding"
head(BXD_Counts_melt_wide)



########## heatmap data wrangling: combine bxd, sherwood and our data for heatmap ########

sherwoodSummedSpectrumObserved_PlusExpectations_melt_FORHEATMAP$dataset <- "Sherwood et al.\n(Y179C/G396D)"
# just select obs/expected:
sherwoodSummedSpectrumObserved_PlusExpectations_melt_FORHEATMAP <- sherwoodSummedSpectrumObserved_PlusExpectations_melt_FORHEATMAP[sherwoodSummedSpectrumObserved_PlusExpectations_melt_FORHEATMAP$variable=="Observed_Over_Expected_Counts",]
## make same columns in each dataset so we can bind: 

sherwoodSummedSpectrumObserved_PlusExpectations_melt_FORHEATMAP$label <- sherwoodSummedSpectrumObserved_PlusExpectations_melt_FORHEATMAP$sherwood_spectrum_label


expectedAndObservedSpectra_SummedOverFamilies_melt_FORHEATMAP$dataset <- "This study\n(Y179C/V344M)"
# select just obs/exp
expectedAndObservedSpectra_SummedOverFamilies_melt_FORHEATMAP <- expectedAndObservedSpectra_SummedOverFamilies_melt_FORHEATMAP[expectedAndObservedSpectra_SummedOverFamilies_melt_FORHEATMAP$variable=="Observed_Over_Expected_Count_Total_SUMMEDPERFAMILY",]
expectedAndObservedSpectra_SummedOverFamilies_melt_FORHEATMAP$label <- expectedAndObservedSpectra_SummedOverFamilies_melt_FORHEATMAP$family

######## heatmap data wrangling: subset dfs for heatmaps and combine #########
# info we need for heatmap: label (yaxis), obs/exp (or D/B BXD/D), mutation type
# subset it for heatmap and melt it, resulting in Mut_Class_DontSeparateCpG variable (D_over_B,) value
BXD_subsetForHeatmap <- melt(select(BXD_Counts_melt_wide,c("Mut_Class_DontSeparateCpG","D_over_B","BXD68_over_B")),variable.name = "label",id.vars = c("Mut_Class_DontSeparateCpG"))
BXD_subsetForHeatmap$dataset <- "Sasani et al.\n(BXD mice)"
# subset our data: 
ourData_perFamily_subsetForHeatmap <- select(expectedAndObservedSpectra_SummedOverFamilies_melt_FORHEATMAP,c("Mut_Class_DontSeparateCpG","variable","value","label","dataset"))

# subset sherwood: only select their MUTYH families?
sherwood_subsetForHeatmap <- select(sherwoodSummedSpectrumObserved_PlusExpectations_melt_FORHEATMAP,c("Mut_Class","variable","value","label","dataset"))
# rename Mut_Class to match other datasets
sherwood_subsetForHeatmap <- sherwood_subsetForHeatmap %>%
  rename("Mut_Class_DontSeparateCpG"="Mut_Class") %>%
  filter(!label %in% c("POLE","POLD1"))

### bind together for heatmap: ###
allDataCombined_ForHeatmap <- bind_rows(ourData_perFamily_subsetForHeatmap,sherwood_subsetForHeatmap,BXD_subsetForHeatmap)  # note 'variable' is NA for bxd and that's ok

# order datasets:
allDataCombined_ForHeatmap$dataset <- factor(allDataCombined_ForHeatmap$dataset,levels=c("This study\n(Y179C/V344M)","Sherwood et al.\n(Y179C/G396D)","Sasani et al.\n(BXD mice)"))

# add some nicer labels:
allDataCombined_ForHeatmap$niceLabel <- allDataCombined_ForHeatmap$label
# nicely label sherwood:
allDataCombined_ForHeatmap[allDataCombined_ForHeatmap$label=="MUTYH_monoallelic",]$niceLabel <- "monoallelic mother+father"
allDataCombined_ForHeatmap[allDataCombined_ForHeatmap$label=="MUTYH_biallelic",]$niceLabel <- "biallelic mother"
# label BXD:
allDataCombined_ForHeatmap[allDataCombined_ForHeatmap$label=="D_over_B",]$niceLabel <- "D strain / B strain"
allDataCombined_ForHeatmap[allDataCombined_ForHeatmap$label=="BXD68_over_B",]$niceLabel <- "BXD68 outlier / B strain"

# organize nice labels:
# need to go in reverse order within each group because heatmap fills in from bottom for some reason
allDataCombined_ForHeatmap$niceLabel <- factor(allDataCombined_ForHeatmap$niceLabel,levels=c("Family 4\n(monoallelic mother)","Family 3\n(biallelic father)","Family 2\n(biallelic mother)","Family 1\n(biallelic mother)","Parent generation\n(monoallelic mother+father)","monoallelic mother+father","biallelic mother","control","BXD68 outlier / B strain","D strain / B strain"))


########### PLOT: heatmap #########
HEATMAP_combinedDatasets <- ggplot(allDataCombined_ForHeatmap,aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=niceLabel,fill=value))+
  geom_tile()+
  facet_grid(dataset~.,scales ="free_y",space="free")+
  geom_text(aes(label=round(value,2)))+
  ylab("")+xlab("")+
  scale_fill_gradient2(midpoint=1,high="red",low="blue",mid="white",na.value = "grey",breaks=c(0.5,1,2,3,4,5,6))+
  theme_classic()+theme(legend.title = element_blank(),text=element_text(size=12),strip.placement = "top")+
  ggtitle("C>A mutation counts appear consistently elevated\nrelative to expectations across datasets")
HEATMAP_combinedDatasets
ggsave(paste0(outdir,"Spectra.observed.expected.HEATMAP.includingBXD.Sherwood.pdf"),HEATMAP_combinedDatasets,height=7,width=8)
ggsave(paste0(outdir,"Spectra.observed.expected.HEATMAP.includingBXD.Sherwood.png"),HEATMAP_combinedDatasets,height=7,width=8)


######### SPECTRA SIGNIFICANCE: per individual ##########


expectedAndObservedSpectra$sigColor <- "non-significant"
expectedAndObservedSpectra[expectedAndObservedSpectra$Poisson_Probability_OfObservingGreaterThanOrEqualToCount<0.05,]$sigColor <- "significant"


poissonProbPlot_PerIndividual <- ggplot(expectedAndObservedSpectra,aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount),color=sigColor))+
  geom_point(size=4,shape=1)+
  facet_wrap(~paste0(child_id,"\n",group))+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Per-individual: mutation types that are significantly elevated\nabove expectations of parental age model")+
  xlab("")+
  ylab("-log10(p-value)")+
  theme(text=element_text(size=14),legend.title=element_blank(),legend.position = "none")+
  scale_color_manual(values=c("black","tomato"))

poissonProbPlot_PerIndividual

ggsave(paste0(outdir,"OurData.PoissonProbs.Full.Spectra.PerIndividual.GreaterThanEqualTo.pdf"),poissonProbPlot_PerIndividual,height=6,width=9)
ggsave(paste0(outdir,"OurData.PoissonProbs.Full.Spectra.PerIndividual.GreaterThanEqualTo.png"),poissonProbPlot_PerIndividual,height=6,width=9)

######### SPECTRA SIGNIFICANCE: per family ##########
expectedAndObservedSpectra_SummedOverFamilies <- read.table(paste0(indir,"expectedAndObservedSpectra_SummedOverFamilies.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),header=T,sep="\t")

expectedAndObservedSpectra_SummedOverFamilies$sigColor <- "non-significant"
expectedAndObservedSpectra_SummedOverFamilies[expectedAndObservedSpectra_SummedOverFamilies$Poisson_Probability_OfObservingGreaterThanOrEqualToCount<0.05,]$sigColor <- "significant"


poissonProbPlot_PerFamily <- ggplot(expectedAndObservedSpectra_SummedOverFamilies,aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount),color=sigColor))+
  geom_point(size=4,shape=1)+
  facet_wrap(~paste0(family))+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Per-family: mutation types that are significantly elevated\nabove expectations of parental age model")+
  xlab("")+
  ylab("-log10(p-value)")+
  theme(text=element_text(size=14),legend.title=element_blank(),legend.position = "none")+
  scale_color_manual(values=c("black","tomato"))

poissonProbPlot_PerFamily
ggsave(paste0(outdir,"OurData.PoissonProbs.Full.Spectra.PerFamily.GreaterThanEqualTo.pdf"),poissonProbPlot_PerFamily,height=6,width=9)
ggsave(paste0(outdir,"OurData.PoissonProbs.Full.Spectra.PerFamily.GreaterThanEqualTo.png"),poissonProbPlot_PerFamily,height=6,width=9)


########### Minimum Effect size (power) analysis plots ########### 
## looking at the minimum C>A counts (and implied mutyh effect sizes) needed to get significant results


##### Min effect sizes per individual (children from Families 1-3 only) ##### 
ExpectedCounts_ForFindingCACountThreshold_plusObserved <- read.table(paste0(indir,"MinimumEffectSize.ExpectedCounts_ForFindingCACountThreshold.Fams1-3Only.txt"),header=T,sep="\t")

minEffectSizePlot_onOverallCACount_Fams1_3 <- ggplot(ExpectedCounts_ForFindingCACountThreshold_plusObserved,aes(x=child_id,y=minCACountForSignficance_qpois_0.05Threshold_UpperTail_plus1))+
  facet_wrap(~family,nrow=1,scales="free_x")+
  geom_point(size=10,shape=95,color="black")+
  geom_point(aes(x=child_id,y=Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio,color="expected C>A count"),size=3)+
  geom_point(aes(x=child_id,y=Observed_Count,color="observed C>A count"),size=3)+
  theme_bw()+
  theme(legend.title=element_blank(),legend.position = "bottom")+
  ggtitle("Observed and expected C>A DNM counts\nBlack lines denote mutator detection thresholds")+
  xlab("")+
  ylab("Overall C>A mutation count\ninherited from both parents")+
  scale_color_manual(values=c("#A6CEE3","#1F78B4"))+
  theme(text=element_text(size=12))

minEffectSizePlot_onOverallCACount_Fams1_3
ggsave(paste0(outdir,"minimumEffectSizesForSignificance.PerIndividual.EffectSizeOnOverallCACount.pdf"),minEffectSizePlot_onOverallCACount_Fams1_3,height=3.5,width=5.5)

# ratios that this impact on overall C>A implies in carrier parent: 
minEffectSizePlot_Fams_1_3_Ratios <- ggplot(ExpectedCounts_ForFindingCACountThreshold_plusObserved,aes(x=child_id,y=minEffectSizeForSignificance_CarrierParent,color="mutator detection threshold"))+
  geom_point(shape=95,size=10,show.legend = F)+
  facet_wrap(~family,nrow=1,scales="free_x")+
  geom_point(aes(x=child_id,y=EstimateOfObservedCarrierParentRatio_assumingNonCarrierEqualsJonsson,color="estimate of carrier parent contribution (empirical data)"),shape=16,size=3)+
  theme_bw()+
  theme(legend.title=element_blank())+
  ggtitle("Est. enrichment of C>A DNMs from carrier parent\nBlack lines denote mutator detection thresholds")+
  xlab("")+
  ylab("ratio of estimated observed/expected\nC>A count from carrier parent")+
  scale_color_manual(values=c("darkorange3",'black'))+
  guides(shape="none")+
  theme(text=element_text(size=12),legend.position="none")
minEffectSizePlot_Fams_1_3_Ratios
ggsave(paste0(outdir,"minimumEffectSizesForSignificance.PerIndividual.AllMutTypes.Ratios.pdf"),minEffectSizePlot_Fams_1_3_Ratios,height=3,width=5)



##### Min effect sizes per family (Families 1-3 only) + Sherwood ##### 
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved <- read.table(paste0(indir,"ExpectedCounts_ForFindingCACountThreshold.Fams1-3Only.PERFAMILY.PlusSherwood.txt"),header=T,sep="\t")


# plot showing overall obs and expected C>A count with mutator detection thresholds (PER FAMILY)
minEffectSizePlot_onOverallCACount_Fams1_3_PERFAMILY <- ggplot(ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved,aes(x=family,y=minCACountForSignficance_qpois_0.05Threshold_UpperTail_plus1_PERFAMILY))+
  geom_point(size=10,shape=95,color="black")+
  geom_point(aes(x=family,y=Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY,color="expected C>A count"),size=3)+
  geom_point(aes(x=family,y=Observed_Count_SUMMEDPERFAMILY,color="observed C>A count"),size=3)+
  theme_bw()+
  theme(legend.title=element_blank(),legend.position = "bottom")+
  ggtitle("Observed and expected C>A DNM counts per family\nBlack lines denote mutator detection thresholds")+
  xlab("")+
  ylab("Overall C>A mutation count per family")+
  scale_color_manual(values=c("#A6CEE3","#1F78B4"))+
  theme(text=element_text(size=12))

minEffectSizePlot_onOverallCACount_Fams1_3_PERFAMILY
ggsave(paste0(outdir,"minimumEffectSizesForSignificance.PerFamily.EffectSizeOnOverallCACount.PlusSherwood.pdf"),minEffectSizePlot_onOverallCACount_Fams1_3_PERFAMILY,height=3.5,width=5)

# plot showing effect sizes in carrier parent that the mutator detector threshold implies
minEffectSizePlot_Fams_1_3_Ratios_PERFAMILY <- ggplot(ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved,aes(x=family,y=minEffectSizeForSignificance_CarrierParent_PERFAMILY,color="mutator detection threshold"))+
  geom_point(shape=95,size=10,show.legend = F)+
  geom_point(aes(x=family,y=EstimateOfObservedCarrierParentRatio_assumingNonCarrierEqualsJonsson_PERFAMILY,color="estimate of carrier parent contribution (empirical data)"),shape=16,size=3)+
  theme_bw()+
  theme(legend.title=element_blank())+
  ggtitle("Est. enrichment of C>A DNMs from carrier parent\nBlack lines denote mutator detection thresholds")+
  xlab("")+
  ylab("ratio of estimated observed/expected\nC>A count from carrier parent")+
  scale_color_manual(values=c("darkorange3",'black'))+
  guides(shape="none")+
  theme(text=element_text(size=12),legend.position="none")
minEffectSizePlot_Fams_1_3_Ratios_PERFAMILY
ggsave(paste0(outdir,"minimumEffectSizesForSignificance.PerFamily.AllMutTypes.Ratios.PlusSherwood.pdf"),minEffectSizePlot_Fams_1_3_Ratios_PERFAMILY,height=3.5,width=5.5)



########### SI Figures ########### 
## Supplementary figures 

##### SI: accessible genome ratios #####
denominatorInfo <- read.table(paste0(indir,"denominatorInfo.txt"),header=T,sep="\t")
accBasesPlot <- ggplot(denominatorInfo,aes(x=child_ID,y=Ratio_IndividualAcc_over_Jonsson))+
  geom_col()+
  theme_bw()+
  geom_hline(yintercept=1,linetype="dashed")+
  ylab("accessible genome ratio")+
  xlab("")
accBasesPlot
ggsave(paste0(SIoutdir,"accessibleGenomeRatioPlot.pdf"),accBasesPlot)
ggsave(paste0(SIoutdir,"accessibleGenomeRatioPlot.png"),accBasesPlot)


##### SI: Overall mutation count scatter plots with Sherwood ##### 
colorsForScatterplot=list("control"="#666666","POLE"="#E6AB02","POLD1"="#D95F02","Parent generation"="#A6761D", "MUTYH_Mother monoallelic"=  "#1B9E77","MUTYH_Both monoallelic"=    "#66A61E","MUTYH_Mother biallelic"="#E7298A","MUTYH_Father biallelic"="#7570B3")

# add labels to our data to match sherwood:
totalPerChild_EXPECTED_and_OBSERVED$label2 <- totalPerChild_EXPECTED_and_OBSERVED$label
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$label=="monoallelic mother",]$label2 <- "MUTYH_Mother monoallelic"
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$label=="biallelic father",]$label2 <- "MUTYH_Father biallelic"
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$label=="biallelic mother",]$label2 <- "MUTYH_Mother biallelic"
# parent generation:
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$label=="monoallelic mother+father",]$label2 <- "MUTYH_Both monoallelic" # both parents are monoallelic for probands

sherwoodData_withExpected <- read.table(paste0(indir,"sherwoodData_withExpected.txt"),header=T,sep="\t")
# exclude: MUTYH_C:II.1 the Sherwood somatic bleedthrough individual
sherwoodData_withExpected <- sherwoodData_withExpected[sherwoodData_withExpected$Child.ID!="MUTYH_C:II.1",]

# all sherwood points EXCLUDING MUTYH_C:II.1 as in Sherwood 
compareTotalCounts_OurDataSherwood <- ggplot(sherwoodData_withExpected,aes(x=EXPECTED_TOTAL_MAT_PLUS_PAT,y=De.novo.SNV.burden,color=label2,shape="Sherwood"))+
  geom_point(size=4)+
  geom_point(data=totalPerChild_EXPECTED_and_OBSERVED,aes(x=EXPECTED_totalPerChild_TimesAccGenomeRatio,y=Observed_totalPerChild,color=label2,shape="our dataset"),size=4)+
  geom_text_repel(data=totalPerChild_EXPECTED_and_OBSERVED,aes(x=EXPECTED_totalPerChild_TimesAccGenomeRatio,y=Observed_totalPerChild,label=child_id))+
  geom_abline()+
  xlim(20,85)+ylim(20,250)+
  theme_bw()+
  ylab("Observed total DNM burden")+
  xlab("Expected total DNM burden based on parental ages\n and corrected for accessible genome size")+
  scale_color_manual(values=colorsForScatterplot)+
  theme(legend.title = element_blank())+
  scale_shape_manual(values=c(17,16))+
  theme(text=element_text(size=14))

compareTotalCounts_OurDataSherwood
ggsave(paste0(SIoutdir,"compareObsExpectedCounts.OurData.PlusSherwood.pdf"),compareTotalCounts_OurDataSherwood,height=5,width=7)
ggsave(paste0(SIoutdir,"compareObsExpectedCounts.OurData.PlusSherwood.png"),compareTotalCounts_OurDataSherwood,height=5,width=7)


##### SI: MAYBE CUTTING: Sherwood poisson probabilities: total mutation counts ##### 

# add color if significant
sherwoodData_withExpected$sigColor_TotalCount <- "non-significant"
sherwoodData_withExpected[sherwoodData_withExpected$Poisson_Probability_OfObservingGreaterThanEqualToCount<0.05,]$sigColor_TotalCount <- "significant"

# exclude somatic bleed-through individual MUTYH_C:II.1

SherwoodPoissonTotalCountPlot <- ggplot(sherwoodData_withExpected[sherwoodData_withExpected$Child.ID!="MUTYH_C:II.1",],aes(x= Child.ID,y=-log10(Poisson_Probability_OfObservingGreaterThanEqualToCount),color=sigColor_TotalCount))+
  geom_point(size=4,shape=1)+
  facet_wrap(~label2,scales="free_x")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1))+
  xlab("")+
  ylab("-log10(p-value)")+
  #ggtitle("Poisson probaility of observing a mutation count >= what we observed under the Jonsson null model")+
  scale_color_manual(values=c("black","tomato"))+
  theme(legend.title = element_blank(),legend.position="none")
SherwoodPoissonTotalCountPlot

ggsave(paste0(SIoutdir,"Sherwood.PoissonProb.TotalCounts.pdf"),SherwoodPoissonTotalCountPlot,height=5,width=7)
ggsave(paste0(SIoutdir,"Sherwood.PoissonProb.TotalCounts.png"),SherwoodPoissonTotalCountPlot,height=5,width=7)

##### SI : Poisson Significance testing of total mutation counts >= Jonsson ##### 

# add color if significant
totalPerChild_EXPECTED_and_OBSERVED$sigColor <- "non-significant"
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$Poisson_Probability_OfObservingGreaterThanOrEqualToCount<0.05,]$sigColor <- "significant"

totalCountPoissonProbsPlot <- ggplot(totalPerChild_EXPECTED_and_OBSERVED,aes(x=child_id,y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount),color=sigColor))+
  geom_point(size=4,shape=1)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  xlab("")+
  ylab("-log10(p-value)")+
  #ggtitle("Poisson probaility of observing a mutation count >=\nwhat we observed under the Jonsson null model")+
  scale_color_manual(values=c("black","tomato"))+
  theme(legend.position = "none")

totalCountPoissonProbsPlot

ggsave(paste0(SIoutdir,"OurData.PoissonProbs.mutationtotals.pdf"),totalCountPoissonProbsPlot,height=3,width=4)
ggsave(paste0(SIoutdir,"OurData.PoissonProbs.mutationtotals.png"),totalCountPoissonProbsPlot,height=3,width=4)

##### SI: observed and expected spectra: column plot: counts (SI) ##### 
# use melted spectrum table: 
expectedAndObservedSpectra_melt <- read.table(paste0(indir,"expectedAndObservedSpectra_melt.txt"),sep="\t",header=T)
# exclude P2:
expectedAndObservedSpectra_melt <- expectedAndObservedSpectra_melt[expectedAndObservedSpectra_melt$child_id!="P2",]

head(expectedAndObservedSpectra_melt)

### make nicer variable labels:
expectedAndObservedSpectra_melt$variable_label <- expectedAndObservedSpectra_melt$variable

expectedAndObservedSpectra_melt[expectedAndObservedSpectra_melt$variable=="Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio",]$variable_label <- "Expected count (corrected for parental age and acc. genome size)"

expectedAndObservedSpectra_melt[expectedAndObservedSpectra_melt$variable=="Observed_Count",]$variable_label <- "Observed count"

expectedAndObservedSpectra_melt[expectedAndObservedSpectra_melt$variable=="Expected_FRACTION_total",]$variable_label <- "Expected fraction (corrected for parental age and acc. genome size)"

expectedAndObservedSpectra_melt[expectedAndObservedSpectra_melt$variable=="Observed_Fraction",]$variable_label <- "Observed fraction"


## column plot showing spectrum counts:
spectrum_column_plot_counts <- ggplot(expectedAndObservedSpectra_melt[expectedAndObservedSpectra_melt$variable %in% c("Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio","Observed_Count"),],aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=value,fill=variable_label))+
  geom_col(position="dodge")+
  facet_wrap(~paste0(child_id,"\n",group),scales="free_y",ncol=4)+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Comparing expected mutation counts based on parental age\nand accessible genome size to observed counts")+
  ylab("mutation count")+xlab("mutation type")+
  theme(legend.position="bottom",legend.title = element_blank(),text=element_text(size=14))

spectrum_column_plot_counts

ggsave(paste0(SIoutdir,"Spectra.observed.expected.columnplot.counts.pdf"),spectrum_column_plot_counts,height=6,width=10)
ggsave(paste0(SIoutdir,"Spectra.observed.expected.columnplot.counts.png"),spectrum_column_plot_counts,height=6,width=10)

##### SI: observed and expected spectra: column plot: fractions (SI) ##### 

spectrum_column_plot_fractions <- ggplot(expectedAndObservedSpectra_melt[expectedAndObservedSpectra_melt$variable %in% c("Expected_FRACTION_total","Observed_Fraction"),],aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=value,fill=variable_label))+
  geom_col(position="dodge")+
  facet_wrap(~paste0(child_id,"\n",group),ncol=4)+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  theme(legend.position="bottom",legend.title = element_blank(),text=element_text(size=14))+
  ggtitle("Comparing expected mutation fractions based on parental age\nand accessible genome size to observed fractions")+
  ylab("mutation fraction")+xlab("mutation type")

spectrum_column_plot_fractions

ggsave(paste0(SIoutdir,"Spectra.observed.expected.columnplot.fractions.pdf"),spectrum_column_plot_fractions,height=6,width=10)
ggsave(paste0(SIoutdir,"Spectra.observed.expected.columnplot.fractions.png"),spectrum_column_plot_fractions,height=6,width=10)

##### SI: Sherwood poisson probabilities of elevated spectrum counts #####
sherwoodSummedSpectrumObserved_PlusExpectations <- read.table(paste0(indir,"Sherwood.expectedAndObservedSpectra_SummedOverGroups.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),header=T,sep="\t")

head(sherwoodSummedSpectrumObserved_PlusExpectations)


# add significance coloring: 
sherwoodSummedSpectrumObserved_PlusExpectations$sigColor <- "non-significant"

sherwoodSummedSpectrumObserved_PlusExpectations[sherwoodSummedSpectrumObserved_PlusExpectations$Poisson_Probability_OfObservingGreaterThanEqualToCount<0.05,]$sigColor <- "significant"


SHERWOOD_poissonProbPlot_PerGroup <- ggplot(sherwoodSummedSpectrumObserved_PlusExpectations,aes(x=gsub("\\.",">",Mut_Class),y=-log10(Poisson_Probability_OfObservingGreaterThanEqualToCount),color=sigColor))+
  geom_point(size=4,shape=1)+
  facet_wrap(~sherwood_spectrum_label)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Sherwood et al.: Mutation types that are significantly elevated\nabove expectations of parental age model (spectra summed per group)")+
  xlab("")+ylab("-log10(p-value)")+
  theme(text=element_text(size=14),legend.position = "none") +
  scale_color_manual(values=c("black","red"))
SHERWOOD_poissonProbPlot_PerGroup
ggsave(paste0(SIoutdir,"Sherwood.Spectrum.PoissonProbs.Full.PerGroup.pdf"),SHERWOOD_poissonProbPlot_PerGroup,height=6,width=8)
ggsave(paste0(SIoutdir,"Sherwood.Spectrum.PoissonProbs.Full.PerGroup.png"),SHERWOOD_poissonProbPlot_PerGroup,height=6,width=8)


#####  SI: Total mutations phased to each parent ##### 

# read in data:
# column plots and poisson probabilities
expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY = read.table(paste0(indir,"OurData.PhasingCounts.ObservedExpected.txt"),header=T,sep="\t") # obs and expected counts phased to each parent (overall counts, not spectra)

# going to restrict just to families 1-3 where we have some a prior expectations of more mutations phasing to either parent:
expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY <- expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY[!expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$child_id %in% c("C41","C42"),] # exclude monoallelics


# phasing fractions for column plot:
expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY = read.table(paste0(indir,"OurData.PhasingFractions.ObservedExpected.txt"),header=T,sep="\t") # # obs and expected fractions phased to each parent (overall fractions, not spectra)
# going to restrict just to families 1-3 where we have some a prior expectations of more mutations phasing to either parent: 
expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY <- expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY[!expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$child_id %in% c("C41","C42"),] # exclude monoallelics
 

# overall phasing data with poisson probabilities:
expectedAndObservedCounts_PlusPhasingInfo = read.table(paste0(indir,"OurData.expectedAndObservedCounts_PlusPhasingInfo.plusPoissonProbs.txt"),header=T,sep="\t") # for poisson probabilities plot 
# restrict to just families 1-3:
expectedAndObservedCounts_PlusPhasingInfo <- expectedAndObservedCounts_PlusPhasingInfo[!expectedAndObservedCounts_PlusPhasingInfo$child_id %in% c("C41","C42"),]

##### SI: Phasing mutations to each parent: counts column plot: #####

# make nicer labels for carrier parent status:
expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$niceVariableLabel <- ""

expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$variable=="Observed_count_maternally_phased",]$niceVariableLabel <- "Observed maternally phased count"

expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$variable=="Observed_count_paternally_phased",]$niceVariableLabel <- "Observed paternally phased count"

expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$variable=="Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes",]$niceVariableLabel <- "Expected paternally phased count"

expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$variable=="Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes",]$niceVariableLabel <- "Expected maternally phased count"

# need to order variables so the right ones are next to each other in the plot: 
expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$niceVariableLabel <- factor(expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$niceVariableLabel,levels=c("Expected maternally phased count","Observed maternally phased count","Expected paternally phased count","Observed paternally phased count"))

# nicer labels for carrier parent: 
expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$carrierParentLabel <- ""

expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$child_id %in% c("C11","C12","C21","C22","C23"),]$carrierParentLabel <- "Biallelic mother"

expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$child_id %in% c("C31","C32"),]$carrierParentLabel <- "Biallelic father"

# order carrierParentLabels:
expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$carrierParentLabel <- factor(expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$carrierParentLabel,levels=c("Biallelic mother","Biallelic father"))

phasingCountPlot_ourData <- ggplot(na.omit(expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY),aes(x=child_id,y=value,fill=niceVariableLabel))+
  geom_col(position = "dodge")+
  facet_wrap(~carrierParentLabel,scales="free_x")+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  theme(text = element_text(size=14),legend.title=element_blank())+
  ggtitle("Comparing observed and expected counts of mutations phased\nto each parent. Expectations from parental age model * phasing rate\n(all mutation types combined)")+
  ylab("count of phased mutations\nassigned to each parent")+xlab("")+
  theme(legend.position = "bottom",text=element_text(size=14))+
  guides(fill = guide_legend(ncol = 2)) # make legend two columns
phasingCountPlot_ourData

ggsave(paste0(SIoutdir,"OurData.PhasedToEachParent.Counts.Families1-3Only.pdf"),phasingCountPlot_ourData,height=4,width=7)

##### SI: Phasing mutations to each parent: fractions column plot #####

# make nicer labels for carrier parent status:
expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$niceVariableLabel <- ""

expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$variable=="Observed_fraction_maternally_phased",]$niceVariableLabel <- "Observed maternally phased fraction"

expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$variable=="Observed_fraction_paternally_phased",]$niceVariableLabel <- "Observed paternally phased fraction"

expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$variable=="Expected_Paternal_Fraction_Phasing_AllTypes",]$niceVariableLabel <- "Expected paternally phased fraction"

expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$variable=="Expected_Maternal_Fraction_Phasing_AllTypes",]$niceVariableLabel <- "Expected maternally phased fraction"

# need to order variables so the right ones are next to each other in the plot: 
expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$niceVariableLabel <- factor(expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$niceVariableLabel,levels=c("Expected maternally phased fraction","Observed maternally phased fraction","Expected paternally phased fraction","Observed paternally phased fraction"))

# nicer labels for carrier parent: 
expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$carrierParentLabel <- ""

expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$child_id %in% c("C11","C12","C21","C22","C23"),]$carrierParentLabel <- "Biallelic mother"

expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$child_id %in% c("C31","C32"),]$carrierParentLabel <- "Biallelic father"

# order carrierParentLabels:
expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$carrierParentLabel <- factor(expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$carrierParentLabel,levels=c("Biallelic mother","Biallelic father"))

phasingFractionPlot_ourData <- ggplot(na.omit(expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY),aes(x=child_id,y=value,fill=niceVariableLabel))+
  geom_col(position = "dodge")+
  facet_wrap(~carrierParentLabel,scales="free_x")+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  theme(text = element_text(size=14),legend.title=element_blank())+
  ggtitle("Comparing observed and expected proportions of mutations phased\nto each parent. Expectations from parental age model * phasing rate\n(all mutation types combined)")+
  ylab("fraction of phased mutations\nassigned to each parent")+xlab("")+
  theme(legend.position = "bottom",text=element_text(size=14))+
  guides(fill = guide_legend(ncol = 2)) # make legend two columns
phasingFractionPlot_ourData

ggsave(paste0(SIoutdir,"OurData.PhasedToEachParent.Fractions.Families1-3Only.pdf"),phasingFractionPlot_ourData,height=4,width=7)

##### SI: Phasing mutations to each parent: Poisson probabilites #####
poissonPlot_TotalPhasedToEachParent <- 
  ggplot(expectedAndObservedCounts_PlusPhasingInfo,aes(x=child_id,y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal)))+
  geom_point(aes(color="maternally phased"),size=4,shape=1)+
  geom_point(aes(x=child_id,y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal),color="paternally phased"),size=4,shape=1)+
  theme_bw()+
  theme(legend.title=element_blank(),text=element_text(size=14))+
  ylab("-log10(p-value)")+xlab("")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  ggtitle("Testing for significant deviations from expectations\nin the number of mutations phased to each parent (all mutation types)")

poissonPlot_TotalPhasedToEachParent

ggsave(paste0(SIoutdir,"OurData.PhasedToEachParent.PoissonProbs.AllMutationTypes.pdf"),poissonPlot_TotalPhasedToEachParent,height=4, width=7)

##### SI: Sherwood counts phased to each parent (MUTYH_Biallelic only) #####

sherwoodData_withExpected # dataframe for plotting Sherwood poisson probabilties (already read in above)

sherwoodData_withExpected_melt <- read.table(paste0(indir,"sherwoodData_withExpected_melt.ForPlottingPhasedColumnPlots.txt"),header=T,sep="\t") # dataframe for plotting Sherwood phasing fractions column plot 
# restrict just to mutyh biallelic:
sherwoodData_withExpected_melt <- sherwoodData_withExpected_melt[sherwoodData_withExpected_melt$label=="MUTYH_Mother biallelic",]

# make nicer variable labels:

#### YOU ARE HERE! 
# make nicer labels for carrier parent status:
sherwoodData_withExpected_melt$niceVariableLabel <- ""

sherwoodData_withExpected_melt[sherwoodData_withExpected_melt$variable=="Observed_Maternal_Count_phased",]$niceVariableLabel <- "Observed maternally phased count"

sherwoodData_withExpected_melt[sherwoodData_withExpected_melt$variable=="Observed_Paternal_Count_phased",]$niceVariableLabel <- "Observed paternally phased count"

sherwoodData_withExpected_melt[sherwoodData_withExpected_melt$variable=="Expected_Paternal_Count_TOTAL_timesProportionPhasedInSherwood",]$niceVariableLabel <- "Expected paternally phased count"

sherwoodData_withExpected_melt[sherwoodData_withExpected_melt$variable=="Expected_Maternal_Count_TOTAL_timesProportionPhasedInSherwood",]$niceVariableLabel <- "Expected maternally phased count"

# need to order variables so the right ones are next to each other in the plot: 
sherwoodData_withExpected_melt$niceVariableLabel <- factor(sherwoodData_withExpected_melt$niceVariableLabel,levels=c("Expected maternally phased count","Observed maternally phased count","Expected paternally phased count","Observed paternally phased count"))


sherwoodPhasedToEachParent_obsExpectedPlot <- ggplot(sherwoodData_withExpected_melt,aes(x=Child.ID,y=value,fill=niceVariableLabel))+
  geom_col(position='dodge')+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Sherwood et al. observed counts phased to each parent,\nwith expectations under parental age model")+
  facet_wrap(~label,scales="free" )+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  xlab("")+ylab("count of mutations phased to each parent")+
  theme(legend.title = element_blank(),text=element_text(size=14),legend.position = "bottom")+
  guides(fill = guide_legend(ncol = 2)) # make legend two columns

sherwoodPhasedToEachParent_obsExpectedPlot
ggsave(paste0(SIoutdir,"Sherwood.CountsPhasedToEachParent.ColumnPlot.MutyhBiallelicOnly.pdf"),sherwoodPhasedToEachParent_obsExpectedPlot,height = 5, width=7)
ggsave(paste0(SIoutdir,"Sherwood.CountsPhasedToEachParent.ColumnPlot.MutyhBiallelicOnly.png"),sherwoodPhasedToEachParent_obsExpectedPlot,height = 5, width=7)

### add fractions: 
####### YOU ARE HERE: show fractions! (are in Sherwood table) # need to calc for Jonsson too if we want expectations
# get expected fractions phased to each parent (already have expected counts)
sherwoodData_withExpected$Expected_Maternal_Phasing_Fraction <- sherwoodData_withExpected$Expected_Maternal_Count_TOTAL_timesProportionPhasedInSherwood / (sherwoodData_withExpected$Expected_Maternal_Count_TOTAL_timesProportionPhasedInSherwood + sherwoodData_withExpected$Expected_Paternal_Count_TOTAL_timesProportionPhasedInSherwood)


sherwoodData_withExpected$Expected_Paternal_Phasing_Fraction <- sherwoodData_withExpected$Expected_Paternal_Count_TOTAL_timesProportionPhasedInSherwood / (sherwoodData_withExpected$Expected_Maternal_Count_TOTAL_timesProportionPhasedInSherwood + sherwoodData_withExpected$Expected_Paternal_Count_TOTAL_timesProportionPhasedInSherwood)

sherwoodData_withExpected_MeltForPhasingFractionPlot <- melt(select(sherwoodData_withExpected,"Child.ID","label","Proportion.phased.to.mother","Proportion.phased.to.father","Expected_Maternal_Phasing_Fraction"  ,"Expected_Paternal_Phasing_Fraction"),id.vars = c("Child.ID","label"))

# make nice variable names:
sherwoodData_withExpected_MeltForPhasingFractionPlot$niceVariableLabel <- ""
sherwoodData_withExpected_MeltForPhasingFractionPlot[sherwoodData_withExpected_MeltForPhasingFractionPlot$variable=="Proportion.phased.to.mother",]$niceVariableLabel <- "Observed maternally phased fraction"
sherwoodData_withExpected_MeltForPhasingFractionPlot[sherwoodData_withExpected_MeltForPhasingFractionPlot$variable=="Proportion.phased.to.father",]$niceVariableLabel <- "Observed paternally phased fraction"
sherwoodData_withExpected_MeltForPhasingFractionPlot[sherwoodData_withExpected_MeltForPhasingFractionPlot$variable=="Expected_Maternal_Phasing_Fraction",]$niceVariableLabel <- "Expected maternally phased fraction"
sherwoodData_withExpected_MeltForPhasingFractionPlot[sherwoodData_withExpected_MeltForPhasingFractionPlot$variable=="Expected_Paternal_Phasing_Fraction",]$niceVariableLabel <- "Expected paternally phased fraction"
# order variables:
sherwoodData_withExpected_MeltForPhasingFractionPlot$niceVariableLabel <- factor(sherwoodData_withExpected_MeltForPhasingFractionPlot$niceVariableLabel, levels=c("Expected maternally phased fraction","Observed maternally phased fraction","Expected paternally phased fraction","Observed paternally phased fraction"))

# mutyh only: 
sherwoodPhasedToEachParent_obsExpectedPlot_FRACTIONS <- ggplot(sherwoodData_withExpected_MeltForPhasingFractionPlot[sherwoodData_withExpected_MeltForPhasingFractionPlot$Child.ID %in% c("MUTYH_B:II.1" ,"MUTYH_B:II.2"),],aes(x=Child.ID,y=value,fill=niceVariableLabel))+
  geom_col(position='dodge')+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Sherwood et al. observed fractions phased to each parent,\nwith expectations under parental age model")+
  facet_wrap(~label,scales="free" )+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  xlab("")+ylab("fraction of mutations phased to each parent")+
  theme(legend.title = element_blank(),text=element_text(size=14),legend.position = "bottom")+
  guides(fill = guide_legend(ncol = 2)) # make legend two columns

sherwoodPhasedToEachParent_obsExpectedPlot_FRACTIONS
ggsave(paste0(SIoutdir,"Sherwood.FractionsPhasedToEachParent.ColumnPlot.MutyhBiallelicOnly.pdf"),sherwoodPhasedToEachParent_obsExpectedPlot_FRACTIONS,height = 5, width=7)
ggsave(paste0(SIoutdir,"Sherwood.FractionsPhasedToEachParent.ColumnPlot.MutyhBiallelicOnly.png"),sherwoodPhasedToEachParent_obsExpectedPlot_FRACTIONS,height = 5, width=7)

##### SI: Sherwood Poisson probabilities overall mutations phased to each parent ######
## restrict to just MUTYH_biallelic 
Sherwood_TotalPhasedTOEachParent_PoissonPlot_MutyhBiallelicOnly <- ggplot(sherwoodData_withExpected[sherwoodData_withExpected$sherwood_spectrum_label=="MUTYH_biallelic",],aes(x=Child.ID,y=-log10(Poisson_Probability_OfObservingGreaterThanEqualToCount_Maternal)))+
  geom_point(aes(color="maternally phased"),size=4,shape=1)+
  geom_point(aes(x=Child.ID,y=-log10(Poisson_Probability_OfObservingGreaterThanEqualToCount_Paternal),color="paternally phased"),size=4,shape=1)+
  theme_bw()+
  facet_wrap(~label,scales="free")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme(axis.text.x = element_text(angle = 45,hjust=1),legend.title = element_blank(),text=element_text(size=14))+
  ggtitle("Testing for significant deviations from expectations\nin the number of mutations phased to each parent (all mutation types)\n(Sherwood et al. biallelic MUTYH family)")+
  xlab("")+ylab("-log10(p-value)")
Sherwood_TotalPhasedTOEachParent_PoissonPlot_MutyhBiallelicOnly

ggsave(paste0(SIoutdir,"Sherwood.PoissonProbs.TotalPhasedToEachParent.AllTypes.pdf"),Sherwood_TotalPhasedTOEachParent_PoissonPlot_MutyhBiallelicOnly,height = 4, width=5)
ggsave(paste0(SIoutdir,"Sherwood.PoissonProbs.TotalPhasedToEachParent.AllTypes.png"),Sherwood_TotalPhasedTOEachParent_PoissonPlot_MutyhBiallelicOnly,height = 4, width=5)

##### SI: Plotting per-family phased spectra as heatmap ##### 
# just families 1-3
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies <- read.table(paste0(indir,"PHASED.expectedAndObservedSpectra_PerFamily.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),header=T,sep='\t')

# want to melt so can plot mother and father phased plots side by side
# fams 1-3 only
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies_melt_filterFORHEATMAP <-
  expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies  %>%
  filter(family %in% c("Family 1\n(biallelic mother)" ,  "Family 2\n(biallelic mother)"  , "Family 3\n(biallelic father)" )) %>%
  select(family,Mut_Class_DontSeparateCpG,Observed_Over_Expected_Count_Maternal_SUMMEDPERFAMILY,Observed_Over_Expected_Count_Paternal_SUMMEDPERFAMILY ) %>%
  melt()
  
# add nice labels:
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies_melt_filterFORHEATMAP$niceVariableLabel <- ""

expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies_melt_filterFORHEATMAP[expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies_melt_filterFORHEATMAP$variable=="Observed_Over_Expected_Count_Maternal_SUMMEDPERFAMILY",]$niceVariableLabel <- "Maternally phased: observed/expected"

expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies_melt_filterFORHEATMAP[expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies_melt_filterFORHEATMAP$variable=="Observed_Over_Expected_Count_Paternal_SUMMEDPERFAMILY",]$niceVariableLabel <- "Paternally phased: observed/expected"

phasedSpectraPerFamilyPlot <- ggplot(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies_melt_filterFORHEATMAP,aes(x=gsub("\\.", ">",Mut_Class_DontSeparateCpG),y=family,fill=value))+
  facet_wrap(~niceVariableLabel)+
  geom_tile()+
  geom_text(aes(label=round(value,2)))+
  xlab("")+ylab("")+
  scale_fill_gradient2(midpoint=1,high="red",low="blue",mid="white",na.value = "grey",breaks=c(0,0.5,1,2,3))+
  theme_classic()+theme(legend.title = element_blank(),text=element_text(size=12),strip.placement = "top")+
  ggtitle("Per-family: Comparing observed and expected\nspectra of muations phased to each parent")
  
phasedSpectraPerFamilyPlot
ggsave(paste0(SIoutdir,"OurData.PhasedSpectra.Heatmap.pdf"),phasedSpectraPerFamilyPlot,height = 4, width=8)
ggsave(paste0(SIoutdir,"OurData.PhasedSpectra.Heatmap.png"),phasedSpectraPerFamilyPlot,height = 4, width=8)

##### SI: Poisson probabilities of phased spectra #########
# just fams 1-3
poissonPlot_PhasedSpectra_PerFamily <- ggplot(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies[expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$family %in% c("Family 1\n(biallelic mother)" ,  "Family 2\n(biallelic mother)"  , "Family 3\n(biallelic father)"),],aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal),color="maternally phased",shape=mat_significant))+
  geom_point(size=4)+
  geom_point(aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal),color="paternally phased",shape=pat_significant),size=4)+
  facet_wrap(~family)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("")+
  xlab("")+
  ylab("-log10(p-value)")+
  theme(text=element_text(size=14),legend.title = element_blank(),legend.position = "bottom")+
  scale_shape_manual(values=c(1,8))+guides(shape="none") #exclude shape from legend

poissonPlot_PhasedSpectra_PerFamily


ggsave(paste0(SIoutdir,"OurData.PhasedSpectra.PoissonProbabilities.pdf"),poissonPlot_PhasedSpectra_PerFamily,height = 4, width=8)
ggsave(paste0(SIoutdir,"OurData.PhasedSpectra.PoissonProbabilities.png"),poissonPlot_PhasedSpectra_PerFamily,height = 4, width=8)

