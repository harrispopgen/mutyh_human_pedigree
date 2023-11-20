##################### Code to generate expected counts of mutations of different types based on parental age ##############
require(dplyr)
# uses the Poisson count regressions from Jonsson et al's Table 9
require(ggplot2)
require(reshape2)
require(ggrepel)
require(RColorBrewer)
require(tidyr)

# quick note about rounding: you don't round expected counts from jonsson model. the only place rounding happens is when expectations become 'observations' which should be integers. but lambda remains real numbers.
# Sherwood observations for phased data (not overall spectra) should be rounded since they are based on multiplying %s phased and so aren't integers (their overall DNM count is an integer, as are their summed up spectra)
todaysdate=Sys.Date()
outdir=paste0("/Users/annabelbeichman/Documents/UW/Human_MUTYH/results/jonssonRegressions/",todaysdate,"_Results_MinVAF_bugfix_IGVFiltered_MinEffectSizeAnalysis/script_output/")

dir.create(outdir,showWarnings = F)
dataIndir="/Users/annabelbeichman/Documents/UW/Human_MUTYH/results/Unfazed/230912_minimal_phased_mutations_IGV_filtered/" # results from Unfazed, minimal VAF, bugfix, have been IGV filtered
jonssonTableS9 <- read.table("/Users/annabelbeichman/Documents/UW/Human_MUTYH/information/Jonsson_nature24018-s2/Jonsson.Table9.ForMutyhProject.Wide.ForR.txt",header=T) # note this has had the T>N nucleotides rev-comped. was manually typed up in excel from Jonsson's table S9 and checked visually. (T>A became A>T; T>C became A>G; T>G became A>C)
# includes CpGs for now!
jonssonTableS9

# carrier parent status: 
biallelicMother=c("C11","C12","C23","C21","C22") # two families combined
biallelicFather=c("C31","C32")
monoallelicMother=c("C41","C42")
monoallelic_BothParents=c("P1","P2","P3","P4")

# families:
family1=c("C11","C12")
family2=c("C23","C21","C22")
family3=c("C31","C32")
family4=c("C41","C42")
parentGeneration=c("P1","P2","P3","P4") # note P2 excluded from most analyses below due to somatic load

# function that takes in this table, maternal age and paternal age (at conception of child) and returns mutation counts
# this is if you have individual mat/pat ages:
expectedDNMSpectrumBasedOnParentalAgesFunc <- function(slopeInterceptTable=jonssonTableS9,mat_age,pat_age) {
  
  # y = mx + b --> b = intercept; m = slope; x = parental age
  slopeInterceptTable$Expected_Maternal_Count <- (slopeInterceptTable$Maternal_Slope*mat_age) + slopeInterceptTable$Maternal_Intercept 
  
  slopeInterceptTable$Expected_Paternal_Count <- (slopeInterceptTable$Paternal_Slope*pat_age) + slopeInterceptTable$Paternal_Intercept 
  
  # get proportions: 
    slopeInterceptTable$Expected_Maternal_Fraction <- slopeInterceptTable$Expected_Maternal_Count/sum(slopeInterceptTable$Expected_Maternal_Count)
  
    slopeInterceptTable$Expected_Paternal_Fraction <- slopeInterceptTable$Expected_Paternal_Count/sum(slopeInterceptTable$Expected_Paternal_Count)

    # for book-keeping, add a record of the mat and pat ages
    slopeInterceptTable$mat_age <- mat_age
    slopeInterceptTable$pat_age <- pat_age

  return(select(slopeInterceptTable,c(Mut_Class,Expected_Maternal_Count,Expected_Maternal_Fraction,Expected_Paternal_Count,Expected_Paternal_Fraction,mat_age,pat_age)))
  
}

# test the function: 
#test = expectedDNMSpectrumBasedOnParentalAgesFunc(jonssonTableS9,mat_age=33,pat_age=37)
#test
#sum(test$Expected_Maternal_Count)
#sum(test$Expected_Paternal_Count)
# okay, father has 2.9x more mutations, that seems about right. 
## note that this includes CpG>TpG as its own category (get added to C>T later in the script)

#  for manuscript: calculate ~min/max mut loads (min: 15 age, max: 50 age)
#test15 = expectedDNMSpectrumBasedOnParentalAgesFunc(jonssonTableS9,mat_age=15,pat_age=15)
#sum(test15$Expected_Maternal_Count)+sum(test15$Expected_Paternal_Count) # 37.2

#test50 = expectedDNMSpectrumBasedOnParentalAgesFunc(jonssonTableS9,mat_age=50,pat_age=50)
#sum(test50$Expected_Maternal_Count)+sum(test50$Expected_Paternal_Count) # 99.5



######## read in empirical mat and pat ages and apply function: ######
mat_pat_agesEmpirical <- read.table("/Users/annabelbeichman/Documents/UW/Human_MUTYH/information/ParentalAgesAtConception/parental_ages.txt",header = T)
mat_pat_agesEmpirical


child_IDs = unique(mat_pat_agesEmpirical$child_ID)

expectedMutationCountsList=list()
for(child in child_IDs){
  mat_age=mat_pat_agesEmpirical[mat_pat_agesEmpirical$child_ID==child,]$mat_age
  pat_age=mat_pat_agesEmpirical[mat_pat_agesEmpirical$child_ID==child,]$pat_age
  
  childResult <- expectedDNMSpectrumBasedOnParentalAgesFunc(jonssonTableS9,mat_age=mat_age,pat_age=pat_age)
  
  expectedMutationCountsList[[child]] <- childResult
}

expectedMutationCounts_df <- bind_rows(expectedMutationCountsList,.id="child_id")

########### sum mat + pat counts for unphased data ##########
expectedMutationCounts_df <- expectedMutationCounts_df %>%
  mutate(Expected_Maternal_Count_PLUS_Expected_Paternal_Count=(Expected_Maternal_Count+Expected_Paternal_Count))

# get the fraction of each mutation type (NOT by summing mat+pat fractions, but by dividing mat+pat counts of each type (mat+pat) by total counts)

expectedMutationCounts_df <- expectedMutationCounts_df %>%
  group_by(child_id) %>%
  mutate(Expected_FRACTION_total = Expected_Maternal_Count_PLUS_Expected_Paternal_Count/sum(Expected_Maternal_Count_PLUS_Expected_Paternal_Count)) %>%
  ungroup()

# total per child: (across all mutation types)
totalPerChild_EXPECTED <- expectedMutationCounts_df %>%
  group_by(child_id) %>%
  summarise(EXPECTED_totalPerChild=sum(Expected_Maternal_Count_PLUS_Expected_Paternal_Count)) %>%
  ungroup()

### NOTE these expected counts get corrected below to deal with differences in total acc. genome in our study vs jonsson #######
####### read in empirical mutations per child #######
empiricalMutationCounts <- read.table(paste0(dataIndir,"/mut_type_counts_noX_IGV_filtered.txt"),header=T,sep="\t")


empiricalMutationCounts <- empiricalMutationCounts %>%
  group_by(child_ID) %>%
  mutate(fraction=count/sum(count)) %>%
  ungroup()

empiricalMutationCounts_TOTALPERCHILD <- empiricalMutationCounts %>%
  group_by(child_ID) %>%
  summarise(Observed_totalPerChild=sum(count)) %>%
  ungroup()

empiricalMutationCounts_TOTALPERCHILD
# note P2 is excluded 
# these counts are green + yellow 
# these match the google sheet
###### fractions from empirical data  ####

# total per individual
##### add CpG>TpG to overall C>T count ####
expectedMutationCounts_df <- expectedMutationCounts_df %>%
  mutate(Mut_Class_DontSeparateCpG = Mut_Class)

expectedMutationCounts_df[expectedMutationCounts_df$Mut_Class=="CpG.TpG",]$Mut_Class_DontSeparateCpG <- "C.T"

expectedMutationCounts_df_CpGsAddedIn <- expectedMutationCounts_df %>%
  group_by(child_id,Mut_Class_DontSeparateCpG,mat_age,pat_age) %>%
  summarise(Expected_Maternal_Count=sum(Expected_Maternal_Count),Expected_Maternal_Fraction=sum(Expected_Maternal_Fraction),Expected_Paternal_Count=sum(Expected_Paternal_Count),Expected_Paternal_Fraction=sum(Expected_Paternal_Fraction),Expected_Maternal_Count_PLUS_Expected_Paternal_Count=sum(Expected_Maternal_Count_PLUS_Expected_Paternal_Count),Expected_FRACTION_total=sum(Expected_FRACTION_total)) %>%
  ungroup()

########### Correct Expectations by differences in denominator between our study and Jonsson ##############
# this table has acc. genome info
denominatorInfo <- read.table(paste0(dataIndir,"total_counts_with_parents_IGV_filtered.txt"),sep="\t",header=T)

JonssonAccessibleBases=2682890000 # 2.6 Gb, from Jonsson et al.

# note don't need factor of 2 bc ratio will cancel out 
denominatorInfo$JonssonAccessibleBases <- JonssonAccessibleBases

denominatorInfo$Ratio_IndividualAcc_over_Jonsson <- denominatorInfo$total_accessible_bases/denominatorInfo$JonssonAccessibleBases # so if our denom is bigger, this ratio will be >1 and count will go up, if our denom is smaller, ratio will be <1 and count will go down

accBasesPlot <- ggplot(denominatorInfo,aes(x=child_ID,y=Ratio_IndividualAcc_over_Jonsson))+
  geom_col()+
  theme_bw()+
  geom_hline(yintercept=1,linetype="dashed")+
  ylab("accessible genome ratio")+
  xlab("")
accBasesPlot
ggsave(paste0(outdir,"accessibleGenomeRatioPlot.pdf"),accBasesPlot)
ggsave(paste0(outdir,"accessibleGenomeRatioPlot.png"),accBasesPlot)

write.table(denominatorInfo,paste0(outdir,"denominatorInfo.txt"),quote=F,sep="\t",row.names=F)

###### add accessible genome info to expected df:
expectedMutationCounts_df_CpGsAddedIn <- merge(expectedMutationCounts_df_CpGsAddedIn,select(denominatorInfo,c(child_ID,Ratio_IndividualAcc_over_Jonsson)),by.x="child_id",by.y="child_ID")

# mat count:
expectedMutationCounts_df_CpGsAddedIn$Expected_Maternal_Count_TimesAccGenomeRatio <- expectedMutationCounts_df_CpGsAddedIn$Expected_Maternal_Count*expectedMutationCounts_df_CpGsAddedIn$Ratio_IndividualAcc_over_Jonsson

# pat count:
expectedMutationCounts_df_CpGsAddedIn$Expected_Paternal_Count_TimesAccGenomeRatio <- expectedMutationCounts_df_CpGsAddedIn$Expected_Paternal_Count*expectedMutationCounts_df_CpGsAddedIn$Ratio_IndividualAcc_over_Jonsson

# also adjust the total mat+pat count:
expectedMutationCounts_df_CpGsAddedIn$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio <- expectedMutationCounts_df_CpGsAddedIn$Expected_Maternal_Count_PLUS_Expected_Paternal_Count*expectedMutationCounts_df_CpGsAddedIn$Ratio_IndividualAcc_over_Jonsson

# also adjust total snp counts:
totalPerChild_EXPECTED <- merge(totalPerChild_EXPECTED,select(denominatorInfo,c(child_ID,Ratio_IndividualAcc_over_Jonsson)),by.x="child_id",by.y="child_ID")

totalPerChild_EXPECTED$EXPECTED_totalPerChild_TimesAccGenomeRatio <- totalPerChild_EXPECTED$EXPECTED_totalPerChild * totalPerChild_EXPECTED$Ratio_IndividualAcc_over_Jonsson


# and add in totals:
totalPerChild_EXPECTED_and_OBSERVED <- merge(totalPerChild_EXPECTED,empiricalMutationCounts_TOTALPERCHILD,by.x="child_id",by.y="child_ID")

############# plot SNP totals vs Jonsson as scatterplot ##########
totalPerChild_EXPECTED_and_OBSERVED$label <- ""
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$child_id %in% c("P1","P2","P3","P4"),]$label <- "Parent generation"
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$child_id %in% biallelicMother,]$label <- "biallelic mother"
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$child_id %in% biallelicFather,]$label <- "biallelic father"
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$child_id %in% monoallelicMother,]$label <- "monoallelic mother"

compareTotalCounts_ToJonsson <- ggplot(totalPerChild_EXPECTED_and_OBSERVED,aes(x=EXPECTED_totalPerChild_TimesAccGenomeRatio,y=Observed_totalPerChild,color=label))+
  geom_point(size=2)+
  geom_abline()+
  geom_text_repel(aes(label=child_id),show.legend = F)+
  theme_bw()+
  ylab("Observed total DNM burden")+
  xlab("Expected total DNM burden based on parental ages\n and corrected for accessible genome size")+
  scale_color_brewer(palette="Dark2")+
  theme(legend.title = element_blank())+
  ggtitle("Comparing total DNM burdens\nto expectations based on parental age\nIGV 'red' calls excluded")+
  theme(text = element_text(size=14))

compareTotalCounts_ToJonsson
ggsave(paste0(outdir,"ourData.compareTotalCounts_ToJonsson.pdf"),compareTotalCounts_ToJonsson,height=5,width=7)


### Get Poisson probabilities on total counts ##########
#P2 is NA so have to remove and make as numeric
totalPerChild_EXPECTED_and_OBSERVED <- totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$child_id!="P2",] # exclude P2 as it is NA

# lower.tail=F means calculating P(X>k) instead of P(X<=k)
# k is observed total -1  ; subtracted  1 to get P(X>=k)=P(X>(k-1)) instead of P(X>k) (works because are integers)
# lambda=expected count from Jonsson (corrected for genome acc)
totalPerChild_EXPECTED_and_OBSERVED$Poisson_Probability_OfObservingGreaterThanOrEqualToCount <- ppois(q=as.numeric(totalPerChild_EXPECTED_and_OBSERVED$Observed_totalPerChild)-1,lambda=totalPerChild_EXPECTED_and_OBSERVED$EXPECTED_totalPerChild_TimesAccGenomeRatio,lower.tail=F) # uses acc genome ratio! 

totalCountPoissonProbsPlot <- ggplot(totalPerChild_EXPECTED_and_OBSERVED,aes(x=child_id,y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount)))+
  geom_point(size=4,shape=1)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Testing whether overall mutation counts are incompatible with Jonsson predictions\nAcc. genome correction applied")+
  xlab("")
totalCountPoissonProbsPlot

ggsave(paste0(outdir,"OurData.PoissonProbs.mutationtotals.pdf"),totalCountPoissonProbsPlot,height=5,width=7)

write.table(totalPerChild_EXPECTED_and_OBSERVED,paste0(outdir,"MutationCounts.totalPerChild_EXPECTED_and_OBSERVED.ContainsPoissonProbs.txt"),quote=F,row.names=F,sep="\t")
##### merge with observed ########

expectedAndObservedSpectra <- merge(expectedMutationCounts_df_CpGsAddedIn,empiricalMutationCounts,by.x=c("child_id","Mut_Class_DontSeparateCpG"),by.y=c("child_ID","mut_type"),all=T)

expectedAndObservedSpectra <- expectedAndObservedSpectra %>%
  rename(Observed_Count=count,Observed_Fraction=fraction)

# add obs/expected counts
expectedAndObservedSpectra$Observed_Over_Expected_Count_Total <- expectedAndObservedSpectra$Observed_Count/expectedAndObservedSpectra$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio

# add group info:
head(expectedAndObservedSpectra)
expectedAndObservedSpectra$group <- ""
expectedAndObservedSpectra[expectedAndObservedSpectra$child_id %in% biallelicMother,]$group <- "biallelic mother"
expectedAndObservedSpectra[expectedAndObservedSpectra$child_id %in% biallelicFather,]$group <- "biallelic father"
expectedAndObservedSpectra[expectedAndObservedSpectra$child_id %in% monoallelicMother,]$group <- "monoallelic mother"
expectedAndObservedSpectra[expectedAndObservedSpectra$child_id %in% monoallelic_BothParents,]$group <- "monoallelic mother+father" # adding P1-P4 as part of mono mom+dad group

###### melt it #########
expectedAndObservedSpectra_melt <- melt(expectedAndObservedSpectra)

#################### Plot observed and expected spectra from our data (NO PHASING) ########

#### PLOT COUNTS ########
### plot counts adjusted by denominator 
countPlot <- ggplot(expectedAndObservedSpectra_melt[expectedAndObservedSpectra_melt$variable %in% c("Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio","Observed_Count") & expectedAndObservedSpectra_melt$child_id!="P2",],aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=value,fill=variable))+
  geom_col(position="dodge")+
  facet_wrap(~paste0(child_id,"\n",group),scales="free_y")+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  ggtitle("Comparing expected mutation counts based on parental age at\nconception (Jonsson S9) to observed counts (*unphased*)\n[expectations now corrected for diffs in acc. genome]")+
  ylab("mutation count")+xlab("mutation type (CpG>TpGs added to C>Ts)")+
  theme(legend.position="bottom",legend.title = element_blank(),text=element_text(size=14))

countPlot
ggsave(paste0(outdir,"Expected.Observed.Counts.pdf"),countPlot,height=7,width=11)

facetedcountplot <- ggplot(expectedAndObservedSpectra[expectedAndObservedSpectra$child_id!="P2",],aes(x=Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio,y=Observed_Count,color=group))+
  geom_point()+
  facet_wrap(~gsub("\\.",">",Mut_Class_DontSeparateCpG),scales="free")+
  geom_abline()+
  geom_text_repel(aes(label=child_id))+
  theme_bw()+
  ylab("Observed total SNV burden\n(before phasing)")+
  xlab("Expected total SNV burden based on parental age\n(Jonsson S9; corrected for diffs in acc genome)")+
  scale_color_brewer(palette="Dark2")+
  theme(legend.title = element_blank(),text=element_text(size=14))
facetedcountplot
ggsave(paste0(outdir,"faceted.counts.pdf"),facetedcountplot,height=8,width=12)
##### PLOT FRACTIONS #######
fractionPlot <- ggplot(expectedAndObservedSpectra_melt[expectedAndObservedSpectra_melt$variable %in% c("Expected_FRACTION_total","Observed_Fraction") & expectedAndObservedSpectra$child_id!="P2",],aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=value,fill=variable))+
  geom_col(position="dodge")+
  facet_wrap(~paste0(child_id,"\n",group),ncol=3)+
  theme_bw()+
  scale_fill_brewer(palette="Paired")+
  theme(legend.position="bottom",legend.title = element_blank(),text=element_text(size=14))+
  ggtitle("Comparing expected mutation fractions based on parental age at\nconception (Jonsson S9) to observed fractions")+
  ylab("mutation fraction")+xlab("mutation type (CpG>TpGs added to C>Ts)")
fractionPlot
ggsave(paste0(outdir,"Expected.Observed.Fractions.pdf"),fractionPlot,height=7,width=11)


# exclude P2 
ObsOverExpectedPlot_HEATMAP <- ggplot(expectedAndObservedSpectra_melt[expectedAndObservedSpectra_melt$variable=="Observed_Over_Expected_Count_Total" & expectedAndObservedSpectra_melt$child_id!="P2",],aes(y=child_id,x=gsub(pattern = "\\.",replacement = ">",Mut_Class_DontSeparateCpG),fill=value))+
  geom_tile()+
  facet_wrap(~group,strip.position = "left",scales="free_y")+
  scale_fill_gradient2(midpoint=1,high="red",low="blue",mid="white",na.value = "grey")+  
  geom_text(aes(label=round(value,2)))+
  theme_bw()+
  theme(text=element_text(size=14))+
  xlab("")+
  ylab("")+
  ggtitle("Observed count / Expected count\n(genome acc correction) (unphased)")
ObsOverExpectedPlot_HEATMAP  
ggsave(paste0(outdir,"OurData.Observed_Over_Expected.Counts.Unphased.HEATMAP.noP2.usethis.pdf"),ObsOverExpectedPlot_HEATMAP,height=5,width=10)

write.table(expectedAndObservedSpectra_melt,paste0(outdir,"expectedAndObservedSpectra_melt.txt"),sep="\t",row.names=F,quote=F)
######## Poisson Probabilities: per individual for spectra #########
# https://www.geeksforgeeks.org/a-guide-to-dpois-ppois-qpois-and-rpois-in-r/

# so the function we want to use is ppois(q=observed count,lambda=jonsson expected count,lower.tail=F) # lower.tail=F is VERY IMPORTANT (makes it the equivalent of 1-ppois(q,lambda)). The probability with lower.tail=T is that we observe a q <= our observed q, given lambda, but the probability we want is prob of observing a q > than the q we observed, given lambda so we use lower.tail=F

# because we want X>=x not, X>x, we are going to do cdf(x-1|lambda) instead of cdf(x|lambda)
# should work because are integers ; see here: https://stats.stackexchange.com/questions/495116/how-to-calculate-greater-than-or-equal-to-in-r-using-ppois
expectedAndObservedSpectra$Poisson_Probability_OfObservingGreaterThanOrEqualToCount <- ppois(q=(expectedAndObservedSpectra$Observed_Count-1),lambda=expectedAndObservedSpectra$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio,lower.tail=F)

# no P2
poissonProbPlot_PerIndividual <- ggplot(expectedAndObservedSpectra[expectedAndObservedSpectra$child_id!="P2",],aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount)))+
  geom_point(size=4,shape=1)+
  facet_wrap(~paste0(child_id,"\n",group))+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Poisson prob of observing rand. var > our (observed count-1)\naka X>=x (integers)\ngiven expected Jonsson count\ngenome acc correction applied")+
  xlab("")+
  theme(text=element_text(size=14))

poissonProbPlot_PerIndividual
ggsave(paste0(outdir,"OurData.PoissonProbs.Full.Spectra.PerIndividual.GreaterThanEqualTo.pdf"),poissonProbPlot_PerIndividual,height=8,width=10)

write.table(expectedAndObservedSpectra,paste0(outdir,"expectedAndObservedSpectra_PerIndividual.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),quote=F,row.names=F,sep="\t")

############# Sum up individuals per group and plot spectra ##############

### exclude P2 from sums:
expectedAndObservedSpectra_SummedOverGroups <- expectedAndObservedSpectra %>%
  filter(child_id!="P2") %>%
  group_by(group,Mut_Class_DontSeparateCpG) %>%
  summarise(Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERGROUP=sum(Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio),Observed_Count_SUMMEDPERGROUP=sum(Observed_Count)) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(Expected_Fraction_PerGroup=Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERGROUP/sum(Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERGROUP),Observed_Fraction_PerGroup=Observed_Count_SUMMEDPERGROUP/sum(Observed_Count_SUMMEDPERGROUP),Observed_Over_Expected_Count_Total_SUMMEDPERGROUP=Observed_Count_SUMMEDPERGROUP/Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERGROUP) %>% # add obs/expected
  ungroup()

head(expectedAndObservedSpectra_SummedOverGroups)
### NOTE these are UNPHASED counts so not sept into mat/pat, added exp mat/pat together. no phasing fraction correction. but yes to denominator corretion.
head(expectedAndObservedSpectra_SummedOverGroups)

expectedAndObservedSpectra_SummedOverGroups_melt <- melt(expectedAndObservedSpectra_SummedOverGroups,id.vars = c("group","Mut_Class_DontSeparateCpG"))



head(expectedAndObservedSpectra_SummedOverGroups_melt)
#### plot observed and expected spectra summed over groups ####

unique(expectedAndObservedSpectra_SummedOverGroups_melt$variable)

# order vars:
expectedAndObservedSpectra_SummedOverGroups_melt$variable <- factor(expectedAndObservedSpectra_SummedOverGroups_melt$variable,levels=c("Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERGROUP","Observed_Count_SUMMEDPERGROUP","Expected_Fraction_PerGroup","Observed_Fraction_PerGroup","Observed_Over_Expected_Count_Total_SUMMEDPERGROUP"))



groupSpectraPlot_COUNTS <- ggplot(expectedAndObservedSpectra_SummedOverGroups_melt[expectedAndObservedSpectra_SummedOverGroups_melt$variable %in% c("Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERGROUP","Observed_Count_SUMMEDPERGROUP"),],aes(x=Mut_Class_DontSeparateCpG,y=value,fill=variable))+
  geom_col(position="dodge")+
  theme_bw()+
  facet_wrap(~group)+
  scale_fill_brewer(palette = "Paired")+
  ylab("count (summed over multiple individuals)")+
  theme(legend.position = "bottom",legend.title = element_blank())+
  ggtitle("P2 excluded from counts")
groupSpectraPlot_COUNTS
ggsave(paste0(outdir,"OurData.spectra.unphased.SUMMEDpergroup.Counts.pdf"),groupSpectraPlot_COUNTS,height=7,width=14)


groupSpectraPlot_FRACTIONS <- ggplot(expectedAndObservedSpectra_SummedOverGroups_melt[expectedAndObservedSpectra_SummedOverGroups_melt$variable %in% c("Expected_Fraction_PerGroup","Observed_Fraction_PerGroup"),],aes(x=Mut_Class_DontSeparateCpG,y=value,fill=variable))+
  geom_col(position="dodge")+
  theme_bw()+
  facet_wrap(~group)+
  scale_fill_brewer(palette = "Paired")+
  ylab("fraction (per group)")+
  theme(legend.position = "bottom",legend.title = element_blank())+  ggtitle("P2 excluded from counts")

groupSpectraPlot_FRACTIONS
ggsave(paste0(outdir,"spectra.unphased.SUMMEDpergroup.Fractions.pdf"),groupSpectraPlot_FRACTIONS,height=7,width=14)


groupObsOverExpectedPlot_HEATMAP <- ggplot(expectedAndObservedSpectra_SummedOverGroups_melt[expectedAndObservedSpectra_SummedOverGroups_melt$variable=="Observed_Over_Expected_Count_Total_SUMMEDPERGROUP",],aes(y=group,x=gsub(pattern = "\\.",replacement = ">",Mut_Class_DontSeparateCpG),fill=value))+
  geom_tile()+
  scale_fill_gradient2(midpoint=1,high="red",low="blue",mid="white",na.value = "grey")+  
  geom_text(aes(label=round(value,2)))+
  theme_bw()+
  theme(text=element_text(size=14))+
  xlab("")+
  ylab("")+
  ggtitle("GROUPS (spectra summed over groups)\nObserved count / Expected count\n(genome acc correction) (unphased)\nP2 excluded from counts")
groupObsOverExpectedPlot_HEATMAP  
ggsave(paste0(outdir,"OurData.GROUPS.Observed_Over_Expected.Counts.Unphased.HEATMAP.usethis.pdf"),groupObsOverExpectedPlot_HEATMAP,height=5,width=7)

##################### Poisson probabilities: per GROUP #####################

head(expectedAndObservedSpectra_SummedOverGroups)

# Going to have it be observed because 1-cdf(x|lambda) is X>x but we want X>=x
# since everything are integers, X>=x will be 1-cdf(x-1|lambda)
# see how this changes things
expectedAndObservedSpectra_SummedOverGroups$Poisson_Probability_OfObservingGreaterThanOrEqualToCount <- ppois(q=(expectedAndObservedSpectra_SummedOverGroups$Observed_Count_SUMMEDPERGROUP-1),lambda=expectedAndObservedSpectra_SummedOverGroups$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERGROUP,lower.tail=F)


poissonProbPlot_PerGroup <- ggplot(expectedAndObservedSpectra_SummedOverGroups,aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount)))+
  geom_point(size=4,shape=1)+
  facet_wrap(~group)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Poisson prob of observing X > (our observed count-1)\naka X>= our observed count (integers) \ngiven expected Jonsson count\ngenome acc correction applied\nSPECTRA SUMMED OVER GROUPS")+
  xlab("")+
  theme(text=element_text(size=14))
poissonProbPlot_PerGroup
ggsave(paste0(outdir,"OurData.PoissonProbs.Full.Spectra.PerGroup.GreaterThanEqualTo.pdf"),poissonProbPlot_PerGroup,height=6,width=8)

write.table(expectedAndObservedSpectra_SummedOverGroups,paste0(outdir,"expectedAndObservedSpectra_SummedOverGroups.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),quote=F,row.names=F,sep="\t")

######## group by family ####################

expectedAndObservedSpectra$family <- ""
expectedAndObservedSpectra[expectedAndObservedSpectra$child_id %in% c("C11","C12"),]$family <- "Family 1\n(biallelic mother)"
expectedAndObservedSpectra[expectedAndObservedSpectra$child_id %in% c("C21","C22","C23"),]$family <- "Family 2\n(biallelic mother)"
expectedAndObservedSpectra[expectedAndObservedSpectra$child_id %in% c("C31","C32"),]$family <- "Family 3\n(biallelic father)"
expectedAndObservedSpectra[expectedAndObservedSpectra$child_id %in% c("C41","C42"),]$family <- "Family 4\n(monoallelic mother)"
expectedAndObservedSpectra[expectedAndObservedSpectra$child_id %in% c("P1","P2","P3","P4"),]$family <- "Parent generation\n(monoallelic mother+father)"
# excluding P2 
# sum over families: 
expectedAndObservedSpectra_SummedOverFamilies <- expectedAndObservedSpectra %>%
  filter(child_id!="P2") %>%
  group_by(family,Mut_Class_DontSeparateCpG) %>%
  summarise(Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY=sum(Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio),Observed_Count_SUMMEDPERFAMILY=sum(Observed_Count)) %>%
  ungroup() %>%
  group_by(family) %>%
  mutate(Expected_Fraction_PerFamily=Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY/sum(Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY),Observed_Fraction_PerFamily=Observed_Count_SUMMEDPERFAMILY/sum(Observed_Count_SUMMEDPERFAMILY),Observed_Over_Expected_Count_Total_SUMMEDPERFAMILY=Observed_Count_SUMMEDPERFAMILY/Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY) %>% # add obs/expected
  ungroup()


head(expectedAndObservedSpectra_SummedOverFamilies)

expectedAndObservedSpectra_SummedOverFamilies_melt <- melt(expectedAndObservedSpectra_SummedOverFamilies,id.vars = c("family","Mut_Class_DontSeparateCpG"))



head(expectedAndObservedSpectra_SummedOverFamilies_melt)
#### plot observed and expected spectra summed over groups ####

unique(expectedAndObservedSpectra_SummedOverFamilies_melt$variable)

# order vars:
expectedAndObservedSpectra_SummedOverFamilies_melt$variable <- factor(expectedAndObservedSpectra_SummedOverFamilies_melt$variable,levels=c("Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY","Observed_Count_SUMMEDPERFAMILY","Expected_Fraction_PerFamily","Observed_Fraction_Perfamily","Observed_Over_Expected_Count_Total_SUMMEDPERFAMILY"))

expectedAndObservedSpectra_SummedOverFamilies$Poisson_Probability_OfObservingGreaterThanOrEqualToCount <- ppois(q=(expectedAndObservedSpectra_SummedOverFamilies$Observed_Count_SUMMEDPERFAMILY-1),lambda=expectedAndObservedSpectra_SummedOverFamilies$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY,lower.tail=F)

poissonProbPlot_PerFamily <- ggplot(expectedAndObservedSpectra_SummedOverFamilies,aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount)))+
  geom_point(size=4,shape=1)+
  facet_wrap(~family)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Poisson prob of observing X > (our observed count-1)\naka X>= our observed count (integers) \ngiven expected Jonsson count\ngenome acc correction applied\nSPECTRA SUMMED OVER FAMILIES")+
  xlab("")+
  theme(text=element_text(size=14))
poissonProbPlot_PerFamily
ggsave(paste0(outdir,"OurData.PoissonProbs.Full.Spectra.PerFamily.GreaterThanEqualTo.pdf"),poissonProbPlot_PerFamily,height=6,width=8)

write.table(expectedAndObservedSpectra_SummedOverFamilies,paste0(outdir,"expectedAndObservedSpectra_SummedOverFamilies.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),quote=T,row.names=F,sep="\t") # have quote =T to deal with family labels

####### per family heatmap ########

expectedAndObservedSpectra_SummedOverFamilies_melt <- melt(expectedAndObservedSpectra_SummedOverFamilies,id.vars = c("family","Mut_Class_DontSeparateCpG"))

familyObsOverExpectedPlot_HEATMAP <- ggplot(expectedAndObservedSpectra_SummedOverFamilies_melt[expectedAndObservedSpectra_SummedOverFamilies_melt$variable=="Observed_Over_Expected_Count_Total_SUMMEDPERFAMILY",],aes(y=family,x=gsub(pattern = "\\.",replacement = ">",Mut_Class_DontSeparateCpG),fill=value))+
  geom_tile()+
  scale_fill_gradient2(midpoint=1,high="red",low="blue",mid="white",na.value = "grey",breaks=c(0.5,1,1.5))+  
  geom_text(aes(label=round(value,2)))+
  theme_bw()+
  theme(text=element_text(size=14),legend.title = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("FAMILY (spectra summed over families)\nObserved count / Expected count\n(genome acc correction) (unphased)\nP2 excluded from counts")
familyObsOverExpectedPlot_HEATMAP  
ggsave(paste0(outdir,"OurData.FAMILIES.Observed_Over_Expected.Counts.Unphased.HEATMAP.usethis.pdf"),familyObsOverExpectedPlot_HEATMAP,height=5,width=7)

write.table(expectedAndObservedSpectra_SummedOverFamilies_melt,paste0(outdir,"expectedAndObservedSpectra_SummedOverFamilies_melt.FORHEATMAP.txt"),quote=T,sep="\t",row.names=F) # quote = T to deal with family label
####### compare maternal and paternal expected counts after PHASING with Unfazed (Candice) #######
# 20230724 - Candice added mutation types to tables 
empiricalMutations_PHASED <- read.table(paste0(dataIndir,"/ALLCHILDREN.TotalCounts.PhasedToEachParent.WithPhasingSucessRate.txt"),sep="\t",header=T) # contains info on phasing success etc

# make nicer names of columns:
empiricalMutations_PHASED <- empiricalMutations_PHASED %>%
  rename(Observed_fraction_paternally_phased=fraction_paternally_phased,Observed_fraction_maternally_phased=fraction_maternally_phased,Observed_count_maternally_phased=count_maternally_phased,Observed_count_paternally_phased=count_paternally_phased)
# sum up expected spectra and then merge with phasing success:
head(expectedAndObservedSpectra)

# this is a bit redundant, but is making a version of the expectedAndObservedSpectra df that is ready to be merged with the phasing information: 
expectedAndObservedCounts <- expectedAndObservedSpectra %>%
  mutate(Observed_Count = ifelse(is.na(Observed_Count), 0, Observed_Count)) %>%   # need to replace count NA as count 0 for C21
  group_by(child_id) %>%
  summarise(Observed_Count_BothParents_NoPhasing_AllTypes=sum(Observed_Count),Expected_Maternal_Count_TimesAccGenomeRatio_NoPhasing_AllTypes=sum(Expected_Maternal_Count_TimesAccGenomeRatio),Expected_Paternal_Count_TimesAccGenomeRatio_NoPhasing_AllTypes=sum(Expected_Paternal_Count_TimesAccGenomeRatio),Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_NoPhasing_AllTypes=sum(Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio)) %>%
  ungroup()  ### THIS NOW ACCOUNTS FOR DENOMINATORS


#  exclude P1-P4 here (since not phased:)
expectedAndObservedCounts_PlusPhasingInfo <- merge(expectedAndObservedCounts[!expectedAndObservedCounts$child_id %in% c("P1","P2","P3","P4"),],empiricalMutations_PHASED[!empiricalMutations_PHASED$child_ID %in% c("P1","P2","P3","P4"),],by.x="child_id",by.y="child_ID",all=T) # exclude P1-P4 because not phased

# calculate expected counts based on phasing success percentage
# multiplying expected mother and father counts by the fraction of phasing success to reduce them down
# mult mother by phasing frac
expectedAndObservedCounts_PlusPhasingInfo$Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes <- expectedAndObservedCounts_PlusPhasingInfo$Expected_Maternal_Count_TimesAccGenomeRatio_NoPhasing_AllTypes * expectedAndObservedCounts_PlusPhasingInfo$fraction_phased

# mult father by phasing frac
expectedAndObservedCounts_PlusPhasingInfo$Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes <- expectedAndObservedCounts_PlusPhasingInfo$Expected_Paternal_Count_TimesAccGenomeRatio_NoPhasing_AllTypes * expectedAndObservedCounts_PlusPhasingInfo$fraction_phased

# mult overall count by phasing frac
expectedAndObservedCounts_PlusPhasingInfo$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes <- expectedAndObservedCounts_PlusPhasingInfo$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_NoPhasing_AllTypes * expectedAndObservedCounts_PlusPhasingInfo$fraction_phased

# get fractions phased to each parent in expectation:
expectedAndObservedCounts_PlusPhasingInfo$Expected_Maternal_Fraction_Phasing_AllTypes <- expectedAndObservedCounts_PlusPhasingInfo$Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes/expectedAndObservedCounts_PlusPhasingInfo$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes

expectedAndObservedCounts_PlusPhasingInfo$Expected_Paternal_Fraction_Phasing_AllTypes <- expectedAndObservedCounts_PlusPhasingInfo$Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes/expectedAndObservedCounts_PlusPhasingInfo$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes

### Poisson probabilities of TOTAL obs and expected phased counts to mat and pat #########
expectedAndObservedCounts_PlusPhasingInfo$Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal <- ppois(q=(expectedAndObservedCounts_PlusPhasingInfo$Observed_count_maternally_phased-1),lambda=expectedAndObservedCounts_PlusPhasingInfo$Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes,lower.tail=F) # takes genome acc and phasing success rates into account!


expectedAndObservedCounts_PlusPhasingInfo$Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal <- ppois(q=(expectedAndObservedCounts_PlusPhasingInfo$Observed_count_paternally_phased-1),lambda=expectedAndObservedCounts_PlusPhasingInfo$Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes,lower.tail=F) # takes genome acc and phasing success rates into account!

poissonPlot_TotalPhasedToEachParent <- 
  ggplot(expectedAndObservedCounts_PlusPhasingInfo,aes(x=child_id,y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal)))+
  geom_point(aes(color="maternally phased"),size=4,shape=1)+
  geom_point(aes(x=child_id,y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal),color="paternally phased"),size=4,shape=1)+
  theme_bw()+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  ggtitle("Testing whether total counts phased to each parent is incompatible with Jonsson predictions")

poissonPlot_TotalPhasedToEachParent

ggsave(paste0(outdir,"OurData.PoissonProbs.TotalPhasedToEachParent.allTypesTogether.pdf"),poissonPlot_TotalPhasedToEachParent,height=5, width=7)

write.table(expectedAndObservedCounts_PlusPhasingInfo,paste0(outdir,"OurData.expectedAndObservedCounts_PlusPhasingInfo.plusPoissonProbs.txt"),quote=F,sep = "\t",row.names = F)

#### melt and select only the fractions for plotting ####
expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY <- melt(expectedAndObservedCounts_PlusPhasingInfo,id.vars = c("child_id","pathogenic_carrier_parent"))

expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY <- expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$variable %in% c("Observed_fraction_paternally_phased","Observed_fraction_maternally_phased","Expected_Paternal_Fraction_Phasing_AllTypes","Expected_Maternal_Fraction_Phasing_AllTypes"),]

# reorder vars :
expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$variable <- factor(expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY$variable,levels=c("Expected_Maternal_Fraction_Phasing_AllTypes","Observed_fraction_maternally_phased","Expected_Paternal_Fraction_Phasing_AllTypes","Observed_fraction_paternally_phased"))
### PLOT EXPECTED MATERNAL AND PATERNAL PHASING FRACTIONS WITH OBSERVED #######
# 
phasingFractionPlot_ourData <- ggplot(na.omit(expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY),aes(x=child_id,y=value,fill=variable))+
  geom_col(position = "dodge")+
  facet_wrap(~paste0("carrier-parent: ",pathogenic_carrier_parent),scales="free_x")+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  theme(text = element_text(size=14),legend.title=element_blank())+
  ggtitle("Expected proportions phased to each parent from Jonsson regression * phasing success rate")+
  ylab("fraction")+
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(ncol = 2)) # make legend two columns
phasingFractionPlot_ourData
ggsave(paste0(outdir,"OurData.PhasedToEachParent.Fractions.pdf"),phasingFractionPlot_ourData,height=7,width=12)

write.table(expectedAndObservedCounts_PlusPhasingInfo_melt_FRACTIONSONLY,paste0(outdir,"OurData.PhasingFractions.ObservedExpected.txt"),quote=F,row.names=F,sep="\t")


####### PLOT PHASING COUNTS WITH EXPECTED AND OBSERVED ########
expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY <- melt(expectedAndObservedCounts_PlusPhasingInfo,id.vars = c("child_id","pathogenic_carrier_parent"))

expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY <- expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY[expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$variable %in% c("Observed_count_paternally_phased","Observed_count_maternally_phased","Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes","Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes"),]

# reorder vars :
expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$variable <- factor(expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY$variable,levels=c("Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes","Observed_count_maternally_phased","Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_AllTypes","Observed_count_paternally_phased"))

phasingCountPlot_ourData <- ggplot(na.omit(expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY),aes(x=child_id,y=value,fill=variable))+
  geom_col(position = "dodge")+
  facet_wrap(~paste0("carrier-parent: ",pathogenic_carrier_parent),scales="free_x")+
  scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  ylab("Count (after phasing)")+
  theme(text = element_text(size=14),legend.title=element_blank())+
  ggtitle("Expected counts phased to each parent from Jonsson regression * phasing success rate\ndifferences in denominator corrected for")+
  theme(legend.position = "bottom")+
  xlab("")+
  guides(fill = guide_legend(ncol = 2)) # make legend two columns
phasingCountPlot_ourData
ggsave(paste0(outdir,"OurData.PhasedToEachParent.Counts.pdf"),phasingCountPlot_ourData,height=7,width=12)

write.table(expectedAndObservedCounts_PlusPhasingInfo_melt_COUNTSONLY,paste0(outdir,"OurData.PhasingCounts.ObservedExpected.txt"),quote=F,sep = "\t",row.names = F)

#### write out fractions: 
##### compare expected phased spectra (merge with mutation counts) #########
# note :sherwood doesn't report phased spectra so no comparison with them here
# read in spectra:
ourPhasedSpectra <- read.table(paste0(dataIndir,"/ALLCHILDREN.todate.SpectrumCounts.PhasedToEachParent.long.txt"),header=T,sep="\t")

ourPhasedSpectra_wide <- pivot_wider(ourPhasedSpectra,id_cols = c("child_id","mut_type"),names_from = "parent_origin",values_from = "count")
# rename cols:
ourPhasedSpectra_wide <- ourPhasedSpectra_wide %>%
  rename(Observed_Maternal_Count=maternal,Observed_Paternal_Count=paternal)

# want to compare to the spectra from Jonsson, but need to apply phasing %s:
head(expectedMutationCounts_df_CpGsAddedIn) # ALL THESE COUNTS ARE PRE PHASING SUCCESS RATE!! # Expected_Maternal_Count and Expected_Paternal_Count PER MUT TYPE are from Jonsson but haven't had phasing correction yet. 
# want to mult by phasing success rate
# add in phasing success rates:
expectedAndObservedSpectra_Phased <- merge(select(expectedMutationCounts_df_CpGsAddedIn[expectedMutationCounts_df_CpGsAddedIn$child_id %in% unique(ourPhasedSpectra_wide$child_id),],c("child_id","Mut_Class_DontSeparateCpG","Expected_Maternal_Count_TimesAccGenomeRatio","Expected_Paternal_Count_TimesAccGenomeRatio")),ourPhasedSpectra_wide,by.x=c("child_id","Mut_Class_DontSeparateCpG"),by.y=c("child_id","mut_type"),all=T) # need all = T otherwise you drop some mutation types that are missing BUT IT DOES SUBSET INDIVIDUALS THAT ARE IN PHASED SPECTRA

# add in group info:
expectedAndObservedSpectra_Phased$group <- ""
expectedAndObservedSpectra_Phased[expectedAndObservedSpectra_Phased$child_id %in% biallelicMother,]$group <- "biallelic mother"
expectedAndObservedSpectra_Phased[expectedAndObservedSpectra_Phased$child_id %in% biallelicFather,]$group <- "biallelic father"
expectedAndObservedSpectra_Phased[expectedAndObservedSpectra_Phased$child_id %in% monoallelicMother,]$group <- "monoallelic mother"
# note that P1-P4 are not phased so aren't added to groups here 


# fill in NAs with 0s: (note you previously excluded inds that weren't phased yet so they aren't getting 0s filled in, it's just filling in missing mutation types)
expectedAndObservedSpectra_Phased <- expectedAndObservedSpectra_Phased %>%
  mutate(Observed_Maternal_Count = ifelse(is.na(Observed_Maternal_Count), 0, Observed_Maternal_Count),Observed_Paternal_Count = ifelse(is.na(Observed_Paternal_Count), 0, Observed_Paternal_Count)) # replace with 0 

# rename some columns to make clear what is pre/post phasing :
expectedAndObservedSpectra_Phased <- expectedAndObservedSpectra_Phased %>%
  rename(Expected_Maternal_Count_TimesAccGenomeRatio_PrePhase = Expected_Maternal_Count_TimesAccGenomeRatio,Expected_Paternal_Count_TimesAccGenomeRatio_PrePhase=Expected_Paternal_Count_TimesAccGenomeRatio)
# now merge in phasing success rate info:
expectedAndObservedSpectra_Phased_PlusPhasingInfo <- merge(expectedAndObservedSpectra_Phased,select(empiricalMutations_PHASED[!empiricalMutations_PHASED$child_ID %in% c("P1","P2","P3","P4"),],c("child_ID","fraction_phased","pathogenic_carrier_parent")),by.x="child_id",by.y="child_ID",all=T) # this gives % successes of phasing 

# mult by phasing fraction:
expectedAndObservedSpectra_Phased_PlusPhasingInfo$Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction <- expectedAndObservedSpectra_Phased_PlusPhasingInfo$Expected_Maternal_Count_TimesAccGenomeRatio_PrePhase* expectedAndObservedSpectra_Phased_PlusPhasingInfo$fraction_phased

expectedAndObservedSpectra_Phased_PlusPhasingInfo$Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction <- expectedAndObservedSpectra_Phased_PlusPhasingInfo$Expected_Paternal_Count_TimesAccGenomeRatio_PrePhase* expectedAndObservedSpectra_Phased_PlusPhasingInfo$fraction_phased



# get proportional spectra (for ma/pa inherited mutations separately):
expectedAndObservedSpectra_Phased_PlusPhasingInfo <- expectedAndObservedSpectra_Phased_PlusPhasingInfo %>%
  group_by(child_id) %>%
  mutate(Observed_Maternal_Fraction=Observed_Maternal_Count/sum(Observed_Maternal_Count),Observed_Paternal_Fraction=Observed_Paternal_Count/sum(Observed_Paternal_Count),Expected_Maternal_Fraction=Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction/sum(Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction),Expected_Paternal_Fraction=Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction/sum(Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction)) %>%
  ungroup()
# results in NA for C22 bc doesn't have any maternally phased sites

# also get obs/expected 
### divide obs counts by expected (Corrected for genome and phasing frac)
expectedAndObservedSpectra_Phased_PlusPhasingInfo$Observed_Over_Expected_Maternal <- expectedAndObservedSpectra_Phased_PlusPhasingInfo$Observed_Maternal_Count/expectedAndObservedSpectra_Phased_PlusPhasingInfo$Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction

expectedAndObservedSpectra_Phased_PlusPhasingInfo$Observed_Over_Expected_Paternal <- expectedAndObservedSpectra_Phased_PlusPhasingInfo$Observed_Paternal_Count/expectedAndObservedSpectra_Phased_PlusPhasingInfo$Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction

# skipping P1-P4 bc not phased
# melt: 
expectedAndObservedSpectra_Phased_PlusPhasingInfo_melt <- melt(expectedAndObservedSpectra_Phased_PlusPhasingInfo,id.vars = c("child_id","group","pathogenic_carrier_parent","Mut_Class_DontSeparateCpG"))

# subset COUNT vars we want to plot:
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SPECTRUMCOUNTS_melt <- expectedAndObservedSpectra_Phased_PlusPhasingInfo_melt[expectedAndObservedSpectra_Phased_PlusPhasingInfo_melt$variable %in% c("Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction","Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction","Observed_Maternal_Count","Observed_Paternal_Count"),]

# reorder variables:
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SPECTRUMCOUNTS_melt$variable <- factor(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SPECTRUMCOUNTS_melt$variable, levels=c("Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction","Observed_Maternal_Count","Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction","Observed_Paternal_Count"))

# exclude P1-P4 bc weren't phased (should already be excluded from dataframe)
spectrumPlot1 <- ggplot(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SPECTRUMCOUNTS_melt[!expectedAndObservedSpectra_Phased_PlusPhasingInfo_SPECTRUMCOUNTS_melt$child_id %in% c("P1","P2","P3","P4"),],aes(x=Mut_Class_DontSeparateCpG,y=value,fill=variable))+
  geom_col(position="dodge")+
  facet_wrap(~paste0(child_id,"\n",group))+
  theme_bw()+
  scale_fill_brewer(palette = "Paired")+
  ylab("mutation count")+
  theme(legend.title = element_blank())+
  ggtitle("expected counts from Jonsson adjusted for\ndiffs in acc. genome size and phasing success")
spectrumPlot1
ggsave(paste0(outdir,"OurData.PhasedToEachParent.SPECTRA.pdf"),spectrumPlot1,height=7,width=16)


#### subset to fraction variables we want to plot:
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SPECTRUMFRACTIONS_melt <- expectedAndObservedSpectra_Phased_PlusPhasingInfo_melt[expectedAndObservedSpectra_Phased_PlusPhasingInfo_melt$variable %in% c("Expected_Maternal_Fraction","Expected_Paternal_Fraction","Observed_Maternal_Fraction","Observed_Paternal_Fraction"),]

# reorder variables:
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SPECTRUMFRACTIONS_melt$variable <- factor(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SPECTRUMFRACTIONS_melt$variable, levels=c("Expected_Maternal_Fraction","Observed_Maternal_Fraction","Expected_Paternal_Fraction","Observed_Paternal_Fraction"))


spectrumPlot2 <- ggplot(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SPECTRUMFRACTIONS_melt[!expectedAndObservedSpectra_Phased_PlusPhasingInfo_SPECTRUMFRACTIONS_melt$child_id %in% c("P1","P2","P3","P4"),],aes(x=Mut_Class_DontSeparateCpG,y=value,fill=variable))+
  geom_col(position="dodge")+
  facet_wrap(~paste0(child_id,"\n",group))+
  theme_bw()+
  scale_fill_brewer(palette = "Paired")+
  ylab("mutation fraction")+
  theme(legend.title = element_blank())

spectrumPlot2
ggsave(paste0(outdir,"OurData.PhasedToEachParent.SPECTRA.FRACTIONS.pdf"),spectrumPlot2,height=7,width=16)
# missing all from C22 (none maternally phased) so you get warning: 
#Warning message:
#  Removed 6 rows containing missing values (geom_col). 


########### Plot obs/exp #################
# obs_over_expected_phased_spetrum_plot_OurData <- ggplot(expectedAndObservedSpectra_Phased_PlusPhasingInfo_melt[expectedAndObservedSpectra_Phased_PlusPhasingInfo_melt$variable %in% c("Observed_Over_Expected_Maternal","Observed_Over_Expected_Paternal"),],aes(x=gsub(pattern = "Observed_Over_Expected_",replacement = "",variable),y=value,fill=variable))+
#   geom_col(position="dodge")+
#   facet_grid(~gsub(pattern = "\\.",replacement = ">",Mut_Class_DontSeparateCpG)~paste0(child_id,"\n",group),switch="y")+
#   theme_bw()+
#   theme(legend.title = element_blank(),text=element_text(size=14),axis.text.x=element_text(angle=45,hjust=1),legend.position = "none",strip.placement = "outside")+
#   geom_hline(yintercept = 1,linetype="dashed")+
#   ggtitle("observed count / expected count\n(corr. for acc. genome and phasing succ rate)")+
#   ylab("obs count / expected count")+
#   xlab("")
# obs_over_expected_phased_spetrum_plot_OurData
# ggsave(paste0(outdir,"OurData.PhasedToEachParent.SPECTRA.obs_over_expected.pdf"),obs_over_expected_phased_spetrum_plot_OurData,height=12,width=16)


obs_over_expected_phased_spetrum_plot_OurData_HEATMAP <- ggplot(expectedAndObservedSpectra_Phased_PlusPhasingInfo_melt[expectedAndObservedSpectra_Phased_PlusPhasingInfo_melt$variable %in% c("Observed_Over_Expected_Maternal","Observed_Over_Expected_Paternal") & !expectedAndObservedSpectra_Phased_PlusPhasingInfo_melt$child_id %in% c("P1","P2","P3","P4"),],aes(y=child_id,x=gsub(pattern = "\\.",replacement = ">",Mut_Class_DontSeparateCpG),fill=value))+
  geom_tile()+
  facet_wrap(~variable)+
  geom_text(aes(label=round(value,2)))+
  theme_bw()+
  theme(legend.title = element_blank(),text=element_text(size=14),axis.text.x=element_text(angle=45,hjust=1),strip.placement = "outside")+
  ggtitle("observed count / expected count\n(corr. for acc. genome and phasing succ rate)")+
  ylab("")+
  xlab("")+
  scale_fill_gradient2(midpoint = 1,high="red",low="blue",mid="white")
obs_over_expected_phased_spetrum_plot_OurData_HEATMAP
ggsave(paste0(outdir,"OurData.PhasedToEachParent.SPECTRA.obs_over_expected.HEATMAP.pdf"),obs_over_expected_phased_spetrum_plot_OurData_HEATMAP,height=12,width=15)


############# Poisson probabilities: phased spectra per individual ###########


expectedAndObservedSpectra_Phased_PlusPhasingInfo$Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal <- ppois(q=(expectedAndObservedSpectra_Phased_PlusPhasingInfo$Observed_Maternal_Count-1),lambda=expectedAndObservedSpectra_Phased_PlusPhasingInfo$Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction,lower.tail=F) # takes genome acc and phasing success rates into account!

expectedAndObservedSpectra_Phased_PlusPhasingInfo$Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal <- ppois(q=(expectedAndObservedSpectra_Phased_PlusPhasingInfo$Observed_Paternal_Count-1),lambda=expectedAndObservedSpectra_Phased_PlusPhasingInfo$Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction,lower.tail=F)

poissonPlot_PhasedSpectra_PerIndividual <- ggplot(expectedAndObservedSpectra_Phased_PlusPhasingInfo[!expectedAndObservedSpectra_Phased_PlusPhasingInfo$child_id %in% c("P1","P2","P3","P4"),],aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal),color="maternally phased"))+
  geom_point(size=4,shape=1)+
  geom_point(aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal),color="paternally phased"),size=4,shape=1)+
  facet_wrap(~paste0(child_id,"\n",group))+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Poisson prob of observing rand. var > our (observed count-1)\naka >= our observed count\ngiven expected Jonsson count\ngenome acc correction AND phasing success rate corrections applied")+
  xlab("")+
  ylab("-log10(Poisson prob of observing greater or equal count given Jonsson)")+
  theme(text=element_text(size=14))
poissonPlot_PhasedSpectra_PerIndividual
ggsave(paste0(outdir,"OurData.PoissonProbs.Phased.Spectra.PerIndividual.GreaterThanEqualto.pdf"),poissonPlot_PhasedSpectra_PerIndividual,height=8,width=9)

write.table(expectedAndObservedSpectra_Phased_PlusPhasingInfo,paste0(outdir,"PHASED.expectedAndObservedSpectra_Phased_PerIndividual.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),quote=F,sep="\t",row.names = F)

#### sum up these phased spectra per group #######


# sum up per group, and exclude P1-P4 that aren't part of any group(?)
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups <- expectedAndObservedSpectra_Phased_PlusPhasingInfo %>%
  filter(!child_id %in% c("P1","P2","P3","P4")) %>% 
  group_by(group,Mut_Class_DontSeparateCpG) %>%
  summarise(Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERGROUP=sum(Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction),Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERGROUP=sum(Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction),Observed_Maternal_Count_SUMMEDPERGROUP=sum(Observed_Maternal_Count),Observed_Paternal_Count_SUMMEDPERGROUP=sum(Observed_Paternal_Count)) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(Observed_Over_Expected_Count_Maternal_SUMMEDPERGROUP=Observed_Maternal_Count_SUMMEDPERGROUP/Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERGROUP,Observed_Over_Expected_Count_Paternal_SUMMEDPERGROUP=Observed_Paternal_Count_SUMMEDPERGROUP/Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERGROUP) %>% # add obs/expected
  ungroup()

expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups_melt <- melt(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups,id.vars = c("group","Mut_Class_DontSeparateCpG"))

######## heatmap plot of phased spectra (obs/expected) ##########
GROUP_obs_over_expected_phased_spetrum_plot_OurData_HEATMAP <- ggplot(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups_melt[expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups_melt$variable %in% c("Observed_Over_Expected_Count_Maternal_SUMMEDPERGROUP","Observed_Over_Expected_Count_Paternal_SUMMEDPERGROUP"),],aes(y=group,x=gsub(pattern = "\\.",replacement = ">",Mut_Class_DontSeparateCpG),fill=value))+
  geom_tile()+
  facet_wrap(~variable)+
  geom_text(aes(label=round(value,2)))+
  theme_bw()+
  theme(legend.title = element_blank(),text=element_text(size=14),axis.text.x=element_text(angle=45,hjust=1),strip.placement = "outside")+
  ggtitle("observed count / expected count\n(corr. for acc. genome and phasing succ rate)\nphased spectra summed per group")+
  ylab("")+
  xlab("")+
  scale_fill_gradient2(midpoint = 1,high="red",low="blue",mid="white")
GROUP_obs_over_expected_phased_spetrum_plot_OurData_HEATMAP
ggsave(paste0(outdir,"OurData.PhasedToEachParent.SPECTRA.obs_over_expected.HEATMAP.PERGROUP.pdf"),GROUP_obs_over_expected_phased_spetrum_plot_OurData_HEATMAP,height=6,width=15)

######## Poisson probabilities: phased spectra per group ###########
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups$Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal <- ppois(q=(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups$Observed_Maternal_Count_SUMMEDPERGROUP-1),lambda=expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups$Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERGROUP,lower.tail=F) # takes genome acc and phasing success rates into account!

expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups$Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal <- ppois(q=(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups$Observed_Paternal_Count_SUMMEDPERGROUP-1),lambda=expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups$Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERGROUP,lower.tail=F) # takes genome acc and phasing success rates into account!


poissonPlot_PhasedSpectra_PerGroup <- ggplot(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups,aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal),color="maternally phased"))+
  geom_point(size=4,shape=1)+
  geom_point(aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal),color="paternally phased"),size=4,shape=1)+
  facet_wrap(~group)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Poisson prob of observing rand. var > our (observed count-1)\naka X >=obs count\ngiven expected Jonsson count\ngenome acc correction AND phasing success rate corrections applied\nSUMMED OVER GROUPS")+
  xlab("")+
  ylab("-log10(Poisson prob of observing greater or equal count given Jonsson)")+
  theme(text=element_text(size=14))
poissonPlot_PhasedSpectra_PerGroup
ggsave(paste0(outdir,"OurData.PoissonProbs.Phased.Spectra.PerGroup.GreaterThanEqualTo.pdf"),poissonPlot_PhasedSpectra_PerGroup,height=5,width=9)

write.table(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverGroups,paste0(outdir,"PHASED.expectedAndObservedSpectra_PerGroup.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),quote=F,row.names=F,sep="\t")

########## group phased spectra by family ##############
expectedAndObservedSpectra_Phased_PlusPhasingInfo$family <- ""
expectedAndObservedSpectra_Phased_PlusPhasingInfo[expectedAndObservedSpectra_Phased_PlusPhasingInfo$child_id %in% c("C11","C12"),]$family <- "Family 1\n(biallelic mother)"
expectedAndObservedSpectra_Phased_PlusPhasingInfo[expectedAndObservedSpectra_Phased_PlusPhasingInfo$child_id %in% c("C21","C22","C23"),]$family <- "Family 2\n(biallelic mother)"
expectedAndObservedSpectra_Phased_PlusPhasingInfo[expectedAndObservedSpectra_Phased_PlusPhasingInfo$child_id %in% c("C31","C32"),]$family <- "Family 3\n(biallelic father)"
expectedAndObservedSpectra_Phased_PlusPhasingInfo[expectedAndObservedSpectra_Phased_PlusPhasingInfo$child_id %in% c("C41","C42"),]$family <- "Family 4\n(monoallelic mother)"

# sum up per family, and exclude P1-P4 that aren't part of any group(?)
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies <- expectedAndObservedSpectra_Phased_PlusPhasingInfo %>%
  filter(!child_id %in% c("P1","P2","P3","P4")) %>% 
  group_by(family,Mut_Class_DontSeparateCpG) %>%
  summarise(Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERFAMILY=sum(Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction),Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERFAMILY=sum(Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction),Observed_Maternal_Count_SUMMEDPERFAMILY=sum(Observed_Maternal_Count),Observed_Paternal_Count_SUMMEDPERFAMILY=sum(Observed_Paternal_Count)) %>%
  ungroup() %>%
  group_by(family) %>%
  mutate(Observed_Over_Expected_Count_Maternal_SUMMEDPERFAMILY=Observed_Maternal_Count_SUMMEDPERFAMILY/Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERFAMILY,Observed_Over_Expected_Count_Paternal_SUMMEDPERFAMILY=Observed_Paternal_Count_SUMMEDPERFAMILY/Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERFAMILY) %>% # add obs/expected
  ungroup()


expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal <- ppois(q=(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$Observed_Maternal_Count_SUMMEDPERFAMILY-1),lambda=expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$Expected_Maternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERFAMILY,lower.tail=F) # takes genome acc and phasing success rates into account!

expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal <- ppois(q=(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$Observed_Paternal_Count_SUMMEDPERFAMILY-1),lambda=expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$Expected_Paternal_Count_TimesAccGenomeRatio_TimesPhasingFraction_SUMMEDPERFAMILY,lower.tail=F) # takes genome acc and phasing success rates into account!

# significance cutoff:
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$mat_significant <- "no"
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$pat_significant <- "no"
expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies[expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal<0.05,]$mat_significant <- "yes"

expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies[expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies$Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal<0.05,]$pat_significant <- "yes"

poissonPlot_PhasedSpectra_PerFamily <- ggplot(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies,aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Maternal),color="maternally phased",shape=mat_significant))+
  geom_point(size=4)+
  geom_point(aes(x=gsub("\\.",">",Mut_Class_DontSeparateCpG),y=-log10(Poisson_Probability_OfObservingGreaterThanOrEqualToCount_Paternal),color="paternally phased",shape=pat_significant),size=4)+
  facet_wrap(~family)+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Poisson prob of observing rand. var > our (observed count-1)\naka X >=obs count\ngiven expected Jonsson count\ngenome acc correction AND phasing success rate corrections applied\nSUMMED OVER FAMILIES")+
  xlab("")+
  ylab("-log10(Poisson prob of observing greater or equal count given Jonsson)")+
  theme(text=element_text(size=14))+
  scale_shape_manual(values=c(1,8))
poissonPlot_PhasedSpectra_PerFamily
ggsave(paste0(outdir,"OurData.PoissonProbs.Phased.Spectra.PerFamily.GreaterThanEqualTo.pdf"),poissonPlot_PhasedSpectra_PerFamily,height=8,width=9)

write.table(expectedAndObservedSpectra_Phased_PlusPhasingInfo_SummedOverFamilies,paste0(outdir,"PHASED.expectedAndObservedSpectra_PerFamily.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),quote=T,row.names=F,sep="\t") # quote = T to deal with family labels

############# want to plug in Sherwood data as well #############
# from Sherwood et al 2023 table 1. In Excel I already multiplied their fractions by counts to get mat/pat counts.

sherwoodData <- read.table("/Users/annabelbeichman/Documents/UW/Human_MUTYH/information/SherwoodMutyhPaper2023/TableOfDNMsandPhasing.txt",header=T,sep="\t")
head(sherwoodData)

SHERWOOD_IDs=unique(sherwoodData$Child.ID)
SHERWOOD_expectedMutationCountsList=list()
for(sherwood_child in SHERWOOD_IDs){
  mat_age=sherwoodData[sherwoodData$Child.ID==sherwood_child,]$Age.mother
  pat_age=sherwoodData[sherwoodData$Child.ID==sherwood_child,]$Age.father
  
  childResult <- expectedDNMSpectrumBasedOnParentalAgesFunc(jonssonTableS9,mat_age=mat_age,pat_age=pat_age)
  
  SHERWOOD_expectedMutationCountsList[[sherwood_child]] <- childResult
}

SHERWOOD_expectedMutationCounts_df <- bind_rows(SHERWOOD_expectedMutationCountsList,.id="child_id")

# get totals per child:
SHERWOOD_expectedMutationCounts_df_TOTALS = SHERWOOD_expectedMutationCounts_df %>%
  group_by(child_id) %>%
  summarise(Expected_Maternal_Count_TOTAL=sum(Expected_Maternal_Count),Expected_Paternal_Count_TOTAL=sum(Expected_Paternal_Count)) %>%
  ungroup() %>%
  mutate(EXPECTED_TOTAL_MAT_PLUS_PAT=Expected_Maternal_Count_TOTAL+Expected_Paternal_Count_TOTAL)

 # merge with observed
sherwoodData_withExpected <- merge(sherwoodData,SHERWOOD_expectedMutationCounts_df_TOTALS,by.x="Child.ID",by.y="child_id",all=T)

######### multiply jonsson expected totals by proportion phased SNVs in sherwood (not the mat/pat specific %s phased, but the proportion of overall snps phased. so if there were 70 snps originally, and 28% phased --> and Jonsson expects 40 dad and 30 mom, 40 and 30 would both by mult by .28....)
sherwoodData_withExpected$EXPECTED_TOTAL_MAT_PLUS_PAT_timesProportionPhasedInSherwood = sherwoodData_withExpected$Proportion.phased.SNVs * sherwoodData_withExpected$EXPECTED_TOTAL_MAT_PLUS_PAT

sherwoodData_withExpected$Expected_Maternal_Count_TOTAL_timesProportionPhasedInSherwood = sherwoodData_withExpected$Expected_Maternal_Count_TOTAL * sherwoodData_withExpected$Proportion.phased.SNVs

sherwoodData_withExpected$Expected_Paternal_Count_TOTAL_timesProportionPhasedInSherwood = sherwoodData_withExpected$Expected_Paternal_Count_TOTAL * sherwoodData_withExpected$Proportion.phased.SNVs

####### Poisson probs: sherwood obs and expected to each parent ######

# prob of observing >= maternal count 
# 20231115: adding rounding because Sherwood phased observations aren't integers due to being calculated based on phasing %s in their Table 1 (total counts are ok, and summed up spectra are ok)
sherwoodData_withExpected$Poisson_Probability_OfObservingGreaterThanEqualToCount_Maternal <- ppois(q=(round(sherwoodData_withExpected$Observed_Maternal_Count_phased)-1),lambda=sherwoodData_withExpected$Expected_Maternal_Count_TOTAL_timesProportionPhasedInSherwood,lower.tail=F)

# prob of observing >= paternal count 
# 20231115: adding rounding because Sherwood phased observations aren't integers due to being calculated based on phasing %s in their Table 1 (total counts are ok, and summed up spectra are ok)

sherwoodData_withExpected$Poisson_Probability_OfObservingGreaterThanEqualToCount_Paternal <- ppois(q=(round(sherwoodData_withExpected$Observed_Paternal_Count_phased)-1),lambda=sherwoodData_withExpected$Expected_Paternal_Count_TOTAL_timesProportionPhasedInSherwood,lower.tail=F)

Sherwood_TotalPhasedTOEachParent_PoissonPlot <- ggplot(sherwoodData_withExpected,aes(x=Child.ID,y=-log10(Poisson_Probability_OfObservingGreaterThanEqualToCount_Maternal)))+
  geom_point(aes(color="maternally phased"),size=4,shape=1)+
  geom_point(aes(x=Child.ID,y=-log10(Poisson_Probability_OfObservingGreaterThanEqualToCount_Paternal),color="paternally phased"),size=4,shape=1)+
  theme_bw()+
  facet_wrap(~label,scales="free")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  ggtitle("Testing whether Sherwood counts phased to each\nparent are sig. >= from Jonsson predictions")+
  xlab("")
Sherwood_TotalPhasedTOEachParent_PoissonPlot
ggsave(paste0(outdir,"Sherwood.PoissonProbs.TotalPhasedToEachParent.AllTypes.pdf"),Sherwood_TotalPhasedTOEachParent_PoissonPlot,height = 6, width=10)

####### get sherwood ratios of observed/expected counts phased to each parent #######

sherwoodData_withExpected$RATIO_Observed_Maternal_Count_phased_OVER_Expected_Maternal_Count_TOTAL_timesProportionPhasedInSherwood <- sherwoodData_withExpected$Observed_Maternal_Count_phased / sherwoodData_withExpected$Expected_Maternal_Count_TOTAL_timesProportionPhasedInSherwood

sherwoodData_withExpected$RATIO_Observed_Paternal_Count_phased_OVER_Expected_Paternal_Count_TOTAL_timesProportionPhasedInSherwood <- sherwoodData_withExpected$Observed_Paternal_Count_phased / sherwoodData_withExpected$Expected_Paternal_Count_TOTAL_timesProportionPhasedInSherwood

View(sherwoodData_withExpected)

################# melt the data for plotting ############################ 
sherwoodData_withExpected_melt <- melt(select(sherwoodData_withExpected,c("Child.ID", "Parent.with.germline.DNA.repair.variant","Observed_Maternal_Count_phased","Observed_Paternal_Count_phased"   ,"label","sherwood_spectrum_label","Expected_Maternal_Count_TOTAL_timesProportionPhasedInSherwood","Expected_Paternal_Count_TOTAL_timesProportionPhasedInSherwood"))) # excluding "Expected_Maternal_Count_TOTAL" , "Expected_Paternal_Count_TOTAL" because they are too high due to not accounting for phasing dropout (I think? need to talk with KH about how Jonsson deals with this)

head(sherwoodData_withExpected_melt)
unique(sherwoodData_withExpected_melt$variable)

sherwoodData_withExpected_melt$variable <- factor(sherwoodData_withExpected_melt$variable,levels=c("Expected_Maternal_Count_TOTAL_timesProportionPhasedInSherwood","Observed_Maternal_Count_phased","Expected_Paternal_Count_TOTAL_timesProportionPhasedInSherwood","Observed_Paternal_Count_phased","Observed_Over_Expected_Count_Maternal","Observed_Over_Expected_Count_Paternal"))

sherwoodPlot1 <- ggplot(sherwoodData_withExpected_melt,aes(x=Child.ID,y=value,fill=variable))+
  geom_col(position='dodge')+
  scale_fill_brewer(palette="Paired")+
  ggtitle("multiplied Jonsson expected counts by % of snps per individual that were phased in Sherwood")+
  facet_wrap(~label,scales="free" )+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+
  xlab("child ID")+ylab("mutation count")+
  theme(legend.title = element_blank())
sherwoodPlot1
ggsave(paste0(outdir,"Sherwood.ExpectedUnderJonsson.vsObserved.CorrectedForPhasingPercent.pdf"),sherwoodPlot1,height=8,width=14)

# plot absolute totals:
sherwoodPlot2 <- ggplot(sherwoodData_withExpected,aes(x=EXPECTED_TOTAL_MAT_PLUS_PAT,y=De.novo.SNV.burden,color=label))+
  geom_point(size=2)+
  geom_abline()+
  xlim(0,250)+ylim(0,250)+
  theme_bw()+
  ylab("Observed total SNV burden\n(before phasing)")+
  xlab("Expected total SNV burden\n(Jonsson S9)")+
  scale_color_brewer(palette="Dark2")+
  theme(legend.title = element_blank())

sherwoodPlot2
ggsave(paste0(outdir,"Sherwood.TotalCounts.ExpectedUnderJonsson.vsObserved.pdf"),sherwoodPlot2,height=4,width=6)


write.table(sherwoodData_withExpected_melt,paste0(outdir,"sherwoodData_withExpected_melt.ForPlottingPhasedColumnPlots.txt"),row.names=F,quote=T,sep="\t")

#### plot heatmap of ratio of phased 
############ compare expected sherwood spectra to Table S2 in sherwood ##############
#Get expected counts per type from Jonsson (have)
#Multiply by phase % from Sherwood
#Add mutyh individuals together and compare to their weird table in SI that has total counts for all individuals combined 

# this has per-mutation type expectations from Jonsson:
head(SHERWOOD_expectedMutationCounts_df)

# add together mat+pat in expectations:
SHERWOOD_expectedMutationCounts_df$EXPECTED_mat_plus_pat_beforePhaseCorrection <- SHERWOOD_expectedMutationCounts_df$Expected_Maternal_Count + SHERWOOD_expectedMutationCounts_df$Expected_Paternal_Count

# need to add labels from expected data:
SHERWOOD_expectedMutationCounts_df_plusLabel <- merge(SHERWOOD_expectedMutationCounts_df, select(sherwoodData,c("Child.ID","sherwood_spectrum_label","Proportion.phased.SNVs")),by.x="child_id",by.y="Child.ID")

##### DO *NOT* multiply by phasing success rate for the spectrum -- it's the non-phased spectrum counts in Sherwood (not separated into mat/pat) ######

#####  *EXCLUDE* weird MUTYH sample from spectrum because that is what Sherwood does: ######
SHERWOOD_expectedMutationCounts_df_plusLabel_excludeC11.1 <- SHERWOOD_expectedMutationCounts_df_plusLabel[!SHERWOOD_expectedMutationCounts_df_plusLabel$child_id=="MUTYH_C:II.1",] # excludes weird mutyh individaul to match sherwood

# sum up spectrum counts per category (POLE,POLD1, etc) ##########

SHERWOOD_expectedMutationCounts_SUMMED_PER_LABEL_TOMATCHSHERWOODSPECTRUM <- SHERWOOD_expectedMutationCounts_df_plusLabel_excludeC11.1 %>%
  group_by(sherwood_spectrum_label,Mut_Class) %>%
  summarise(Expected_Count_AllIndividualsSummed=sum(EXPECTED_mat_plus_pat_beforePhaseCorrection),n=n(),Expected_Maternal_Count_AllIndividualsSummed=sum(Expected_Maternal_Count),Expected_Paternal_Count_AllIndividualsSummed=sum(Expected_Paternal_Count)) %>%
  ungroup()
## adding in summed up mat/pat count expectations

##### add CpGs + C>Ts together for expected sherwood counts: ###
SHERWOOD_expectedMutationCounts_SUMMED_PER_LABEL_TOMATCHSHERWOODSPECTRUM <- SHERWOOD_expectedMutationCounts_SUMMED_PER_LABEL_TOMATCHSHERWOODSPECTRUM %>%
mutate(Mut_Class_DontSeparateCpG = Mut_Class)

# add in CpGs to C>T category:
SHERWOOD_expectedMutationCounts_SUMMED_PER_LABEL_TOMATCHSHERWOODSPECTRUM[SHERWOOD_expectedMutationCounts_SUMMED_PER_LABEL_TOMATCHSHERWOODSPECTRUM$Mut_Class=="CpG.TpG",]$Mut_Class_DontSeparateCpG <- "C.T"

SHERWOOD_expectedMutationCounts_SUMMED_PER_LABEL_TOMATCHSHERWOODSPECTRUM_CpGsAddedIn <- SHERWOOD_expectedMutationCounts_SUMMED_PER_LABEL_TOMATCHSHERWOODSPECTRUM %>%
  group_by(sherwood_spectrum_label,Mut_Class_DontSeparateCpG,n) %>%
  summarise(Expected_Count_AllIndividualsSummed=sum(Expected_Count_AllIndividualsSummed),Expected_Maternal_Count_AllIndividualsSummed=sum(Expected_Maternal_Count_AllIndividualsSummed),Expected_Paternal_Count_AllIndividualsSummed=sum(Expected_Paternal_Count_AllIndividualsSummed)) %>% # this should leave all counts unchanged except adding CpG+C.T
  ungroup()


############ read in sherwood observed (summed) spectrum ###########


sherwoodSummedSpectrumObserved <- read.table("/Users/annabelbeichman/Documents/UW/Human_MUTYH/information/SherwoodMutyhPaper2023/Sherwood.SummedUpSpectrum.TableS2.forR.txt",sep="\t",header=T)
sherwoodSummedSpectrumObserved
sherwoodSummedSpectrumObserved_melt <- melt(sherwoodSummedSpectrumObserved,id.vars = c("sherwood_spectrum_label","n"),variable.name = "Mut_Class",value.name = "Observed_Count")

sherwoodSummedSpectrumObserved_melt

##### merge sherwood spectrum and jonsson expected spectrum #######
# cpgs added to C>T count
sherwoodSummedSpectrumObserved_PlusExpectations <- merge(sherwoodSummedSpectrumObserved_melt,SHERWOOD_expectedMutationCounts_SUMMED_PER_LABEL_TOMATCHSHERWOODSPECTRUM_CpGsAddedIn,by.x=c("sherwood_spectrum_label","n","Mut_Class"),by.y=c("sherwood_spectrum_label","n","Mut_Class_DontSeparateCpG"))

dim(SHERWOOD_expectedMutationCounts_SUMMED_PER_LABEL_TOMATCHSHERWOODSPECTRUM_CpGsAddedIn)
dim(sherwoodSummedSpectrumObserved_melt)
dim(sherwoodSummedSpectrumObserved_PlusExpectations)

# View(sherwoodSummedSpectrumObserved_melt_PlusExpectations)
# add in obs/exp
sherwoodSummedSpectrumObserved_PlusExpectations$Observed_Over_Expected_Counts <- sherwoodSummedSpectrumObserved_PlusExpectations$Observed_Count/sherwoodSummedSpectrumObserved_PlusExpectations$Expected_Count_AllIndividualsSummed

# get obs and expected fractions:
sherwoodSummedSpectrumObserved_PlusExpectations <- sherwoodSummedSpectrumObserved_PlusExpectations %>%
  group_by(sherwood_spectrum_label,n) %>%
  mutate(Expected_Fraction = Expected_Count_AllIndividualsSummed/sum(Expected_Count_AllIndividualsSummed),Observed_Fraction=Observed_Count/sum(Observed_Count)) %>%
  ungroup()

# exclude mat/pat expectations from melting -- just wanted those for Kelley table:
sherwoodSummedSpectrumObserved_PlusExpectations_melt  <- melt(select(sherwoodSummedSpectrumObserved_PlusExpectations,c("sherwood_spectrum_label","n","Mut_Class","Expected_Count_AllIndividualsSummed","Observed_Count","Expected_Fraction",    "Observed_Fraction"    ,"Observed_Over_Expected_Counts")),id.vars = c("sherwood_spectrum_label","n","Mut_Class"),value.name = "value") # just pick the variables you want

sherwoodSummedSpectrumObserved_PlusExpectations_melt$variable <- factor(sherwoodSummedSpectrumObserved_PlusExpectations_melt$variable,levels=c("Expected_Count_AllIndividualsSummed","Observed_Count","Expected_Fraction",    "Observed_Fraction"    ,"Observed_Over_Expected_Counts"))

sherwoodSummedSpectrumObserved_PlusExpectations_melt$n_label <- paste0(sherwoodSummedSpectrumObserved_PlusExpectations_melt$sherwood_spectrum_label," (n = ",as.character(sherwoodSummedSpectrumObserved_PlusExpectations_melt$n),")")

sherwoodSpectrumPlot1 <- ggplot(sherwoodSummedSpectrumObserved_PlusExpectations_melt[sherwoodSummedSpectrumObserved_PlusExpectations_melt$variable %in% c("Expected_Count_AllIndividualsSummed","Observed_Count"),],aes(x=Mut_Class,y=value,fill=variable))+
  geom_col(position="dodge")+
  facet_wrap(~n_label,scales="free")+
  theme_bw()+
  scale_fill_brewer(palette = "Paired")+
  xlab("")+ylab("mutation count")+theme(legend.title = element_blank())+
  ggtitle("Comparing expected (based on Jonsson S9 regression) and observed (Sherwood) spectrum counts\nnote counts are summed over all individuals in each category as in Sherwood Table S2\nand individual MUTYH_CII.1 was excluded as in Sherwood")
sherwoodSpectrumPlot1

ggsave(paste0(outdir,"Sherwood.SPECTRA.summedOverIndividuals.ExpectedUnderJonsson.vsObserved.pdf"),sherwoodSpectrumPlot1,height=5,width=12)


####### plot spectrum fractions as well: #########
sherwoodSpectrumPlot2 <- ggplot(sherwoodSummedSpectrumObserved_PlusExpectations_melt[sherwoodSummedSpectrumObserved_PlusExpectations_melt$variable %in% c("Expected_Fraction",    "Observed_Fraction") ,],aes(x=Mut_Class,y=value,fill=variable))+
  geom_col(position="dodge")+
  facet_wrap(~n_label,scales="free")+
  theme_bw()+
  scale_fill_brewer(palette = "Paired")+
  xlab("")+ylab("mutation fraction")+theme(legend.title = element_blank())+
  ggtitle("Comparing expected (based on Jonsson S9 regression) and observed (Sherwood) spectrum fractions\n individual MUTYH_CII.1 was excluded as in Sherwood")
sherwoodSpectrumPlot2

ggsave(paste0(outdir,"Sherwood.SPECTRA.FRACTIONS.summedOverIndividuals.ExpectedUnderJonsson.vsObserved.pdf"),sherwoodSpectrumPlot2,height=5,width=12)
  


### obs/expected heatmap ######
sherwoodSpectrumPlot3_HEATMAP <- ggplot(sherwoodSummedSpectrumObserved_PlusExpectations_melt[sherwoodSummedSpectrumObserved_PlusExpectations_melt$variable=="Observed_Over_Expected_Counts",],aes(y=sherwood_spectrum_label,x=gsub("\\.",">",Mut_Class),fill=value))+
  geom_tile()+
  theme_bw()+
  scale_fill_gradient2(midpoint = 1,mid="white",high="red",low="blue")+
  geom_text(aes(label=round(value,2)))+
  ylab("")+
  xlab("")+
  theme(text=element_text(size=14))+
  ggtitle("Sherwood obs / expected counts\n(summed up spectra across inds)")
  
sherwoodSpectrumPlot3_HEATMAP
ggsave(paste0(outdir,"Sherwood.SPECTRA.OBS_OVER_EXP.summedOverIndividuals.HEATMAP.pdf"),sherwoodSpectrumPlot3_HEATMAP,height=5,width=7)


sherwoodSpectrumPlot3_HEATMAP_excludePOL <- ggplot(sherwoodSummedSpectrumObserved_PlusExpectations_melt[sherwoodSummedSpectrumObserved_PlusExpectations_melt$variable=="Observed_Over_Expected_Counts" & !sherwoodSummedSpectrumObserved_PlusExpectations_melt$sherwood_spectrum_label %in% c("POLE","POLD1"),],aes(y=sherwood_spectrum_label,x=gsub("\\.",">",Mut_Class),fill=value))+
  geom_tile()+
  theme_bw()+
  scale_fill_gradient2(midpoint = 1,mid="white",high="red",low="blue")+
  geom_text(aes(label=round(value,2)))+
  ylab("")+
  xlab("")+
  theme(text=element_text(size=14))+
  ggtitle("Sherwood obs / expected counts\n(summed up spectra across inds)")

sherwoodSpectrumPlot3_HEATMAP_excludePOL
ggsave(paste0(outdir,"Sherwood.SPECTRA.OBS_OVER_EXP.summedOverIndividuals.HEATMAP.ExclPOL.pdf"),sherwoodSpectrumPlot3_HEATMAP_excludePOL,height=5,width=7)

write.table(sherwoodSummedSpectrumObserved_PlusExpectations_melt,paste0(outdir,"sherwoodSummedSpectrumObserved_PlusExpectations_melt.FORHEATMAP.txt"),row.names = F,sep="\t",quote=F)



### Poisson probabilities: Sherwood full spectrum per group (no per ind available) ################
head(sherwoodSummedSpectrumObserved_PlusExpectations)

# note for sherwood: we don't have an acc genome estimate, but since they didn't do surrogate fams it's probs pretty close to Jonsson as our full families are   
# note summed spectrum obs values are integers, no rounding needed here
sherwoodSummedSpectrumObserved_PlusExpectations$Poisson_Probability_OfObservingGreaterThanEqualToCount <- ppois(q=(sherwoodSummedSpectrumObserved_PlusExpectations$Observed_Count-1),lambda=sherwoodSummedSpectrumObserved_PlusExpectations$Expected_Count_AllIndividualsSummed,lower.tail=F)

SHERWOOD_poissonProbPlot_PerGroup <- ggplot(sherwoodSummedSpectrumObserved_PlusExpectations,aes(x=gsub("\\.",">",Mut_Class),y=-log10(Poisson_Probability_OfObservingGreaterThanEqualToCount)))+
  geom_point(size=4,shape=1)+
  facet_wrap(~sherwood_spectrum_label,scales="free_y")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+
  ggtitle("Poisson prob of observing rand. var > (Sherwood observed count-1)\naka >= obs count\ngiven expected Jonsson count\nSPECTRA SUMMED OVER GROUPS")+
  xlab("")+
  theme(text=element_text(size=14))
SHERWOOD_poissonProbPlot_PerGroup
ggsave(paste0(outdir,"Sherwood.PoissonProbs.Full.Spectra.PerGroup.GreaterThanEqualTo.pdf"),SHERWOOD_poissonProbPlot_PerGroup,height=6,width=8)

write.table(sherwoodSummedSpectrumObserved_PlusExpectations,paste0(outdir,"Sherwood.expectedAndObservedSpectra_SummedOverGroups.PlusPoissonProbabilities.GreaterThanEqualTo.txt"),quote=F,row.names=F,sep="\t")

# write out version for Kelley that just has obs and expected for C>A
write.table(sherwoodSummedSpectrumObserved_PlusExpectations[sherwoodSummedSpectrumObserved_PlusExpectations$Mut_Class=="C.A",c("sherwood_spectrum_label","n","Mut_Class","Observed_Count","Expected_Count_AllIndividualsSummed","Expected_Maternal_Count_AllIndividualsSummed","Expected_Paternal_Count_AllIndividualsSummed")],paste0(outdir,"Sherwood.expectedAndObservedSpectra_SummedOverGroups.C.AOnly.forKelley.txt"),row.names = F,sep="\t",quote=F)

### Compare our results with Sherwood #############

  


# add more general label to Sherwood
sherwoodData_withExpected$label2 <- sherwoodData_withExpected$label
sherwoodData_withExpected[sherwoodData_withExpected$label %in% c("POLD1_Mother"),]$label2 <- "POLD1"
sherwoodData_withExpected[sherwoodData_withExpected$label %in% c("POLE_Mother","POLE_Father"),]$label2 <- "POLE"

totalPerChild_EXPECTED_and_OBSERVED$label2 <- totalPerChild_EXPECTED_and_OBSERVED$label
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$label=="monoallelic mother",]$label2 <- "MUTYH_Mother monoallelic"
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$label=="biallelic father",]$label2 <- "MUTYH_Father biallelic"
totalPerChild_EXPECTED_and_OBSERVED[totalPerChild_EXPECTED_and_OBSERVED$label=="biallelic mother",]$label2 <- "MUTYH_Mother biallelic"


# colors:
# colors from Dark2 palette
colorsForScatterplot=list("control"="#666666","POLE"="#E6AB02","POLD1"="#D95F02","Parent generation"="#A6761D","MUTYH_Mother monoallelic"=       "#1B9E77",       "MUTYH_Both monoallelic"=    "#66A61E","MUTYH_Mother biallelic"="#E7298A","MUTYH_Father biallelic"="#7570B3")


# Get Poisson Probs for Sherwood:
# note these total counts are integers, no rounding needed
sherwoodData_withExpected$Poisson_Probability_OfObservingGreaterThanEqualToCount <- ppois(q=(sherwoodData_withExpected$De.novo.SNV.burden-1),lambda=sherwoodData_withExpected$EXPECTED_TOTAL_MAT_PLUS_PAT,lower.tail=F)
# don't use phasing fraction here ^ just want total comparisons

compareTotalCounts_OurDataSherwood <- ggplot(sherwoodData_withExpected,aes(x=EXPECTED_TOTAL_MAT_PLUS_PAT,y=De.novo.SNV.burden,color=label2,shape="Sherwood"))+
  geom_point(size=4)+
  geom_point(data=totalPerChild_EXPECTED_and_OBSERVED,aes(x=EXPECTED_totalPerChild_TimesAccGenomeRatio,y=Observed_totalPerChild,color=label2,shape="our dataset"),size=4)+
  geom_text_repel(data=totalPerChild_EXPECTED_and_OBSERVED,aes(x=EXPECTED_totalPerChild_TimesAccGenomeRatio,y=Observed_totalPerChild,label=child_id))+
  geom_abline()+
  xlim(20,85)+ylim(20,250)+
  theme_bw()+
  ylab("Observed total DNM burden\n(before phasing)")+
  xlab("Expected total DNM burden\n(Jonsson S9)\nexp counts for our dataset have acc genome correction")+
  scale_color_manual(values=colorsForScatterplot)+
  theme(legend.title = element_blank())+
  scale_shape_manual(values=c(17,16))+
  ggtitle("Comparing our observed & expected counts to Sherwood's")+
  theme(text=element_text(size=14))

compareTotalCounts_OurDataSherwood
ggsave(paste0(outdir,"compareObsExpectedCounts.OurData.PlusSherwood.pdf"),compareTotalCounts_OurDataSherwood,height=8,width=12)


SherwoodPoissonTotalCountPlot <- ggplot(sherwoodData_withExpected,aes(x= Child.ID,y=-log10(Poisson_Probability_OfObservingGreaterThanEqualToCount)))+
  geom_point(size=4,shape=1)+
  facet_wrap(~label2,scales="free")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("Poisson Prob of Sherwood total counts being >= Jonsson expectations")

SherwoodPoissonTotalCountPlot
ggsave(paste0(outdir,"Sherwood.PoissonProb.TotalCounts.pdf"),SherwoodPoissonTotalCountPlot,height=7,width=10)

# facet by whether pol mother or father is carrier:
SherwoodPoissonTotalCountPlot2 <- ggplot(sherwoodData_withExpected,aes(x= Child.ID,y=-log10(Poisson_Probability_OfObservingGreaterThanEqualToCount)))+
  geom_point(size=4,shape=1)+
  facet_wrap(~label,scales="free")+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+theme(axis.text.x = element_text(angle=45,hjust=1))+
  ggtitle("Poisson Prob of Sherwood total counts being >= Jonsson expectations")

SherwoodPoissonTotalCountPlot2
ggsave(paste0(outdir,"Sherwood.PoissonProb.TotalCounts.plot2.pdf"),SherwoodPoissonTotalCountPlot2,height=7,width=10)

compareTotalCounts_OurDataSherwood_MUTYHonly <- ggplot(sherwoodData_withExpected[sherwoodData_withExpected$label2 %in% c("MUTYH_Both mono-allelic" ,"MUTYH_Mother biallelic"),],aes(x=EXPECTED_TOTAL_MAT_PLUS_PAT,y=De.novo.SNV.burden,color=label2,shape="Sherwood"))+
  geom_point(size=4)+
  geom_point(data=totalPerChild_EXPECTED_and_OBSERVED,aes(x=EXPECTED_totalPerChild_TimesAccGenomeRatio,y=Observed_totalPerChild,color=label2,shape="our dataset"),size=4)+
  geom_text_repel(data=totalPerChild_EXPECTED_and_OBSERVED,aes(x=EXPECTED_totalPerChild_TimesAccGenomeRatio,y=Observed_totalPerChild,label=child_id))+
  geom_abline()+
  xlim(20,85)+ylim(20,200)+
  theme_bw()+
  ylab("Observed total SNV burden\n(before phasing)")+
  xlab("Expected total SNV burden\n(Jonsson S9)\nexp counts for our dataset have acc genome correction")+
  scale_color_manual(values=colorsForScatterplot[which(names(colorsForScatterplot) %in% c("MUTYH_Both monoallelic" ,"MUTYH_Mother biallelic","MUTYH_Mother biallelic","MUTYH_Father biallelic","MUTYH_Mother monoallelic", "Parent generation"))])+
  theme(legend.title = element_blank())+
  scale_shape_manual(values=c(17,16))+
  ggtitle("Comparing our observed & expected counts to Sherwood's\nMUTYH only")+
  theme(text=element_text(size=14))

compareTotalCounts_OurDataSherwood_MUTYHonly


ggsave(paste0(outdir,"compareObsExpectedCounts.OurData.PlusSherwood.MUTYHOnly.pdf"),compareTotalCounts_OurDataSherwood_MUTYHonly,height=8,width=12)


##### exclude polE and one extreme mutyh indivdual from sherwood so we can get a closer look

compareTotalCounts_OurDataSherwood_ExcludePOLEAndExtremeMutyh <- ggplot(sherwoodData_withExpected[!sherwoodData_withExpected$label2 %in% c("POLE","MUTYH_Both mono-allelic"),],aes(x=EXPECTED_TOTAL_MAT_PLUS_PAT,y=De.novo.SNV.burden,color=label2,shape="Sherwood"))+
  geom_point(size=4)+
  geom_point(data=totalPerChild_EXPECTED_and_OBSERVED,aes(x=EXPECTED_totalPerChild_TimesAccGenomeRatio,y=Observed_totalPerChild,color=label2,shape="our dataset"),size=4)+
  geom_text_repel(data=totalPerChild_EXPECTED_and_OBSERVED,aes(x=EXPECTED_totalPerChild_TimesAccGenomeRatio,y=Observed_totalPerChild,label=child_id))+
  geom_abline()+
  xlim(20,85)+ylim(20,85)+
  theme_bw()+
  ylab("Observed total SNV burden\n(before phasing)")+
  xlab("Expected total SNV burden\n(Jonsson S9)\nexp counts for our dataset have acc genome correction")+
  scale_color_manual(values=colorsForScatterplot[which(!names(colorsForScatterplot) %in% c("POLE","MUTYH_Both mono-allelic"))])+ # subsetting colors 
  theme(legend.title = element_blank())+
  scale_shape_manual(values=c(8,16))+
  ggtitle("Comparing our observed & expected counts to Sherwood's\nexclude extreme POLE and sherwood mutyh individual")+
  theme(text=element_text(size=14))

compareTotalCounts_OurDataSherwood_ExcludePOLEAndExtremeMutyh


ggsave(paste0(outdir,"compareObsExpectedCounts.OurData.PlusSherwood.ExcludePOLEandOneExtremeMUTYHInd.pdf"),compareTotalCounts_OurDataSherwood_ExcludePOLEAndExtremeMutyh,height=8,width=12)

write.table(sherwoodData_withExpected,paste0(outdir,"sherwoodData_withExpected.txt"),quote=F,row.names=F,sep="\t")



######## Minimum Effect Sizes for Significance : Estimating excess parental contributions ##########
# E_(C>A,total)=E_(C>A,carrier parent)+ E_(C>A,non-carrier parent) 
# The Jonsson parental age model gives us expectations of DNMs coming from each parent
# Observed data:
#O_(C>A,total)=O_(C>A,carrier parent)+ O_(C>A,non-carrier parent) 
#But in our observed data, we don’t know the values of  O_(C>A,carrier parent) or O_(C>A,non-carrier parent) because we can’t assign the mutations as coming from mother or father.
# But let’s assume that all excess C>A is coming from the carrier parent such that 
#O*_(C>A,non-carrier parent)  = E_(C>A,non-carrier parent) , such that an estimate of the mutations in our observed data that comes from that parent is:
#  O*_(C>A,carrier parent) = O_(C>A,total)- E_(C>A,non-carrier parent)

#Then the effect size on the carrier parent can be estimated as the ratio of the estimate #of the observed count, over the expected:
#  Effect size = (O*_(C>A,carrier parent)  )/E_(C>A,carrier parent) 


# am not doing this for Family 4 or parent generation since don't have expectation of what monoallelic variant would do (plus Fam 4 data issues)

ExpectedCounts_ForFindingCACountThreshold <- select(expectedMutationCounts_df_CpGsAddedIn,c("child_id","Mut_Class_DontSeparateCpG","Expected_Maternal_Count_TimesAccGenomeRatio" , "Expected_Paternal_Count_TimesAccGenomeRatio","Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio")) # start with Jonsson Expectations for our families (not their observed counts, just their parental ages are relevant for the expectations); keeping in per-parent counts now (times acc genome ratio)

# filter to just C.A
ExpectedCounts_ForFindingCACountThreshold <- ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$Mut_Class_DontSeparateCpG=="C.A",]

## get qpois() estimates of threshold at 0.05 and add +1 (lower.tail=F) for each individual
ExpectedCounts_ForFindingCACountThreshold$minCACountForSignficance_qpois_0.05Threshold_UpperTail_plus1 <- qpois(0.05,ExpectedCounts_ForFindingCACountThreshold$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio,lower.tail=F) + 1  # +1 because that's the actual count we'd have to OBSERVE to get the p < 0.05, since otherwise is saying is significant if you observe something greater than that 


head(ExpectedCounts_ForFindingCACountThreshold)

## assign carrier parent labels:
ExpectedCounts_ForFindingCACountThreshold$CarrierParent <- ""
ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$child_id %in% c("C11","C12","C21" ,"C22" ,"C23"),]$CarrierParent <- "mother" # # excluding ,"C41" ,"C42" (monoallelic mother)
ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$child_id %in% c( "C31" , "C32"),]$CarrierParent <- "father"
#ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$child_id %in% c("P1" , "P2",  "P3" , "P4"),]$CarrierParent <- "both" # note both are monoallelic

# then want to subtract off expected contribution from non-carrier parent 
ExpectedCounts_ForFindingCACountThreshold$CarrierParentContributionNeededForSignificance <- NA


# when mother is carrier parent: 
ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="mother",]$CarrierParentContributionNeededForSignificance <- ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="mother",]$minCACountForSignficance_qpois_0.05Threshold_UpperTail_plus1 - ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="mother",]$Expected_Paternal_Count_TimesAccGenomeRatio # subtract off paternal contribution here because just want maternal contribution to include all the excess C>A

# when father is carrier parent:
# rename CarrierParentContribution to something more informative
ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="father",]$CarrierParentContributionNeededForSignificance <- ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="father",]$minCACountForSignficance_qpois_0.05Threshold_UpperTail_plus1 - ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="father",]$Expected_Maternal_Count_TimesAccGenomeRatio # subtract off paternal contribution here because just want maternal contribution to include all the excess C>A


ExpectedCounts_ForFindingCACountThreshold$minEffectSizeForSignificance_CarrierParent <- NA

# get ratio of mother carrier parent contribution to expected maternal contribution
ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="mother",]$minEffectSizeForSignificance_CarrierParent <- ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="mother",]$CarrierParentContributionNeededForSignificance/ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="mother",]$Expected_Maternal_Count_TimesAccGenomeRatio
  
# get ratio of father carrier parent contribution to expected paternal contribution 
ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="father",]$minEffectSizeForSignificance_CarrierParent <- ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="father",]$CarrierParentContributionNeededForSignificance/ExpectedCounts_ForFindingCACountThreshold[ExpectedCounts_ForFindingCACountThreshold$CarrierParent=="father",]$Expected_Paternal_Count_TimesAccGenomeRatio

# Calculate what this would do to the OVERALL C>A count ratio (obs/exp C>A)
ExpectedCounts_ForFindingCACountThreshold$overallCACount_Ratio_ThresholdForSignificance <- ExpectedCounts_ForFindingCACountThreshold$minCACountForSignficance_qpois_0.05Threshold_UpperTail_plus1 /ExpectedCounts_ForFindingCACountThreshold$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio

# also want to estimate the 'observed' carrier parent contribution based on the actual total number of C>A mutations observed (whether significant or not), assuming that the non-carrier contribution matches Jonsson:

## merge in observed counts:
## add in observed counts (Observed_Count)and overall C>A ratio(Observed_Over_Expected_Count_Total)
ExpectedCounts_ForFindingCACountThreshold_plusObserved <- merge(ExpectedCounts_ForFindingCACountThreshold,expectedAndObservedSpectra[expectedAndObservedSpectra$Mut_Class_DontSeparateCpG=="C.A" & !expectedAndObservedSpectra$child_id %in% c("P1" , "P2",  "P3" , "P4","C41","C42"),c("family","child_id","Mut_Class_DontSeparateCpG","Observed_Count","Observed_Over_Expected_Count_Total")],by=c("child_id","Mut_Class_DontSeparateCpG"))

##### subtract non-carrier expectation from observed count to get estimate of carrier parent contribution
# counts:
ExpectedCounts_ForFindingCACountThreshold_plusObserved$EstimateOfObservedCarrierParentContribution_assumingNonCarrierEqualsJonsson <- NA

# ratios:
ExpectedCounts_ForFindingCACountThreshold_plusObserved$EstimateOfObservedCarrierParentRatio_assumingNonCarrierEqualsJonsson <- NA

# counts:
ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="mother",]$EstimateOfObservedCarrierParentContribution_assumingNonCarrierEqualsJonsson <- 
  ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="mother",]$Observed_Count - ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="mother",]$Expected_Paternal_Count_TimesAccGenomeRatio # subtract paternal expectation from observed mutation count to get estimate of maternal contribution when mothers are carriers

# ratio:
ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="mother",]$EstimateOfObservedCarrierParentRatio_assumingNonCarrierEqualsJonsson <- 
  ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="mother",]$EstimateOfObservedCarrierParentContribution_assumingNonCarrierEqualsJonsson / ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="mother",]$Expected_Maternal_Count_TimesAccGenomeRatio # divided estimate of carrier parent contribution by their expectation

# count:
ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="father",]$EstimateOfObservedCarrierParentContribution_assumingNonCarrierEqualsJonsson <- 
  ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="father",]$Observed_Count - ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="father",]$Expected_Maternal_Count_TimesAccGenomeRatio # subtract maternal expectation from observed mutation count to get estimate of paternal contribution when fathers are carriers

# ratio
ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="father",]$EstimateOfObservedCarrierParentRatio_assumingNonCarrierEqualsJonsson <- 
  ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="father",]$EstimateOfObservedCarrierParentContribution_assumingNonCarrierEqualsJonsson / ExpectedCounts_ForFindingCACountThreshold_plusObserved[ExpectedCounts_ForFindingCACountThreshold_plusObserved$CarrierParent=="father",]$Expected_Paternal_Count_TimesAccGenomeRatio # divided estimate of carrier parent contribution by their expectation


#### effect size plots (per individual) ######

# overall C>A count:
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




write.table(ExpectedCounts_ForFindingCACountThreshold_plusObserved[!ExpectedCounts_ForFindingCACountThreshold_plusObserved$child_id %in% c("C41","C42","P1","P2","P3","P4"),],paste0(outdir,"MinimumEffectSize.ExpectedCounts_ForFindingCACountThreshold.Fams1-3Only.txt"),quote=T,sep="\t",row.names=F) # need quote=T 


  
########## minimum effect sizes: per family (+Sherwood biallelic mutyh family) ########## 
# figure this out per-family in terms of effect size on carrier parents
# 20231108: add Sherwood mutyh biallelic family to this analysis

# get per-individual counts then sum up per family (note these aren't mult by phasing fraction)
ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis <- select(expectedMutationCounts_df_CpGsAddedIn,c("child_id","Mut_Class_DontSeparateCpG","Expected_Maternal_Count_TimesAccGenomeRatio" , "Expected_Paternal_Count_TimesAccGenomeRatio","Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio"))

# restrict to just C>A and just families 1-3:
ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis <- ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis[ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis$Mut_Class_DontSeparateCpG=="C.A" & !ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis$child_id %in% c("C41","C42","P1","P2","P3","P4"),] 

# add family labels (1-3):
ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis$family <- ""

ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis[ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis$child_id %in% c("C11","C12"),]$family <- "Family 1\n(biallelic mother)"

ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis[ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis$child_id %in% c("C21","C22","C23"),]$family <- "Family 2\n(biallelic mother)"

ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis[ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis$child_id %in% c("C31","C32"),]$family <- "Family 3\n(biallelic father)"



head(ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis)
# sum up expectations per family:
ExpectedCounts_ForFindingCACountThreshold_perFamily <- ExpectedCounts_ForFindingCACountThreshold_ForFamilyAnalysis %>%
  group_by(Mut_Class_DontSeparateCpG,family) %>%
  summarise(Expected_Maternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY=sum(Expected_Maternal_Count_TimesAccGenomeRatio),Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY=sum(Expected_Paternal_Count_TimesAccGenomeRatio),Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY=sum(Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio)) %>%
  ungroup()

##### add in Sherwood biallelic family expectations for min effect size analysis: ######
sherwoodSummedSpectrumObserved_PlusExpectations_ForMinEffectSizeAnalysis <- sherwoodSummedSpectrumObserved_PlusExpectations[sherwoodSummedSpectrumObserved_PlusExpectations$sherwood_spectrum_label=="MUTYH_biallelic" & sherwoodSummedSpectrumObserved_PlusExpectations$Mut_Class=="C.A",c("sherwood_spectrum_label","Mut_Class","Expected_Count_AllIndividualsSummed","Expected_Maternal_Count_AllIndividualsSummed","Expected_Paternal_Count_AllIndividualsSummed","Observed_Count")]

sherwoodSummedSpectrumObserved_PlusExpectations_ForMinEffectSizeAnalysis$family <- "Sherwood et al.\n(biallelic mother)"

# make column names consistent with ExpectedCounts_ForFindingCACountThreshold_perFamily so can bind rows:
sherwoodSummedSpectrumObserved_PlusExpectations_ForMinEffectSizeAnalysis <- sherwoodSummedSpectrumObserved_PlusExpectations_ForMinEffectSizeAnalysis %>%
  rename(Mut_Class_DontSeparateCpG = Mut_Class,Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY=Expected_Count_AllIndividualsSummed,Expected_Maternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY=Expected_Maternal_Count_AllIndividualsSummed,Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY=Expected_Paternal_Count_AllIndividualsSummed,Observed_Count_SUMMEDPERFAMILY=Observed_Count)
# note Sherwood wasn't corrected for acc genome size, just using variable name for convenience for binding rows 


# want to add to ExpectedCounts_ForFindingCACountThreshold_perFamily

ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood <- bind_rows(ExpectedCounts_ForFindingCACountThreshold_perFamily,select(sherwoodSummedSpectrumObserved_PlusExpectations_ForMinEffectSizeAnalysis,names(ExpectedCounts_ForFindingCACountThreshold_perFamily))) # pulling the same columns that are in ExpectedCounts_ForFindingCACountThreshold_perFamily from Sherwood using select() here

# get qpois thresholds for overall C>A count to be significant: (+1)

ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$minCACountForSignficance_qpois_0.05Threshold_UpperTail_plus1_PERFAMILY <- qpois(0.05,ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY,lower.tail=F) + 1 # +1 because that's the actual count we'd have to OBSERVE to get the 0.05, since otherwise is saying is significant if you observe something greater than that 

## assign carrier parent labels:
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent <- ""

ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$family %in% c("Family 1\n(biallelic mother)","Family 2\n(biallelic mother)","Sherwood et al.\n(biallelic mother)"),]$CarrierParent <- "mother" 

ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$family=="Family 3\n(biallelic father)",]$CarrierParent <- "father"

# then want to subtract off expected contribution from non-carrier parent: 
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParentContributionNeededForSignificance_PERFAMILY <- NA


# when mother is carrier parent: subtract the non-carrier father (Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY) contributions (per family) from the per-family qpois sig threshold 
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="mother",]$CarrierParentContributionNeededForSignificance_PERFAMILY <- ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="mother",]$minCACountForSignficance_qpois_0.05Threshold_UpperTail_plus1_PERFAMILY - ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="mother",]$Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY # subtract off paternal contribution from min total C>A needed for significance because just want maternal contribution to include all the excess C>A

# when father is carrier parent:subtract the non-carrier mother (Expected_Maternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY) contributions (per family) from the per-family qpois sig threshold 
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="father",]$CarrierParentContributionNeededForSignificance_PERFAMILY <- ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="father",]$minCACountForSignficance_qpois_0.05Threshold_UpperTail_plus1_PERFAMILY - ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="father",]$Expected_Maternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY # subtract off MATERNAL contribution here because just want paternal contribution to include all the excess C>A


# get ratio effect size in the carrier parent (~"observed" (assuming all excess from carrier)/expectd carrier parent contribution)

ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$minEffectSizeForSignificance_CarrierParent_PERFAMILY <- NA

# divide the carrier parent contribution needed for significance (per fam) by the expected counts from that parent (Jonsson; summed per family)
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="mother",]$minEffectSizeForSignificance_CarrierParent_PERFAMILY <- ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="mother",]$CarrierParentContributionNeededForSignificance_PERFAMILY/ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="mother",]$Expected_Maternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY

# get ratio of father carrier parent contribution to expected paternal contribution 
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="father",]$minEffectSizeForSignificance_CarrierParent_PERFAMILY <- ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="father",]$CarrierParentContributionNeededForSignificance_PERFAMILY/ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$CarrierParent=="father",]$Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY

# Calculate what this would do to the OVERALL C>A count 
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$overallCACount_Ratio_ThresholdForSignificance <- ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$minCACountForSignficance_qpois_0.05Threshold_UpperTail_plus1_PERFAMILY /ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood$Expected_Maternal_Count_PLUS_Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY

# also want to estimate the 'observed' carrier parent contribution based on the actual total number of C>A mutations observed (whether significant or not), assuming that the non-carrier contribution matches Jonsson:
## merge in observed counts from our data:
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved <- merge(ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood,expectedAndObservedSpectra_SummedOverFamilies[expectedAndObservedSpectra_SummedOverFamilies$Mut_Class_DontSeparateCpG=="C.A" & !expectedAndObservedSpectra_SummedOverFamilies$family %in% c("Family 4\n(monoallelic mother)","Parent generation\n(monoallelic mother+father)") ,c("family","Mut_Class_DontSeparateCpG","Observed_Count_SUMMEDPERFAMILY")],by=c("family","Mut_Class_DontSeparateCpG"),all.x = T) # excluding fam4 and probands

# add in Sherwood Observed value: (32 C>A )
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$family=="Sherwood et al.\n(biallelic mother)",]$Observed_Count_SUMMEDPERFAMILY <- sherwoodSummedSpectrumObserved_PlusExpectations_ForMinEffectSizeAnalysis$Observed_Count_SUMMEDPERFAMILY

##### subtract non-carrier expectation from observed count to get estimate of carrier parent contribution
# counts:
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$EstimateOfObservedCarrierParentContribution_assumingNonCarrierEqualsJonsson_PERFAMILY <- NA

# ratios:
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$EstimateOfObservedCarrierParentRatio_assumingNonCarrierEqualsJonsson_PERFAMILY <- NA


ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="mother",]$EstimateOfObservedCarrierParentContribution_assumingNonCarrierEqualsJonsson_PERFAMILY <- 
  ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="mother",]$Observed_Count_SUMMEDPERFAMILY - ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="mother",]$Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY # subtract paternal expectation (summed per family) from observed mutation count (summed per family) to get estimate of maternal contribution when mothers are carriers

# calculate ratio of carrier parent contribution / expectation: 
ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="mother",]$EstimateOfObservedCarrierParentRatio_assumingNonCarrierEqualsJonsson_PERFAMILY <- 
  ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="mother",]$EstimateOfObservedCarrierParentContribution_assumingNonCarrierEqualsJonsson_PERFAMILY / ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="mother",]$Expected_Maternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY # divided estimate of carrier parent contribution by their expectation to get a ratio 


ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="father",]$EstimateOfObservedCarrierParentContribution_assumingNonCarrierEqualsJonsson_PERFAMILY <- 
  ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="father",]$Observed_Count_SUMMEDPERFAMILY - ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="father",]$Expected_Maternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY # subtract maternal expectation from observed mutation count to get estimate of paternal contribution when fathers are carriers (both summed per family)


ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="father",]$EstimateOfObservedCarrierParentRatio_assumingNonCarrierEqualsJonsson_PERFAMILY <- 
  ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="father",]$EstimateOfObservedCarrierParentContribution_assumingNonCarrierEqualsJonsson_PERFAMILY / ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved[ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved$CarrierParent=="father",]$Expected_Paternal_Count_TimesAccGenomeRatio_SUMMEDPERFAMILY # divided estimate of carrier parent contribution by their expectation


#### effect size plots (per family) ######

# overall C>A count:
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
ggsave(paste0(outdir,"minimumEffectSizesForSignificance.PerFamily.EffectSizeOnOverallCACount.PlusSherwood.pdf"),minEffectSizePlot_onOverallCACount_Fams1_3_PERFAMILY,height=3.5,width=4.5)

# ratios that this impact on overall C>A implies in carrier parent: 
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
ggsave(paste0(outdir,"minimumEffectSizesForSignificance.PerFamily.AllMutTypes.Ratios.PlusSherwood.pdf"),minEffectSizePlot_Fams_1_3_Ratios_PERFAMILY,height=3.5,width=4.5)




write.table(ExpectedCounts_ForFindingCACountThreshold_perFamily_PlusSherwood_PlusObserved,paste0(outdir,"ExpectedCounts_ForFindingCACountThreshold.Fams1-3Only.PERFAMILY.plusSherwood.txt"),quote=T,sep="\t",row.names=F) # quote=T to deal with family labels


########### save R environment: ###########
save.image(file = paste0(outdir,todaysdate,".environment.RData"))

