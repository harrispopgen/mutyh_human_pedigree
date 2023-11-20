require(dplyr)
require(reshape2)
require(tidyr)
# code to count up maternally and paternally phased mutations from Candice's unfazed tables

indir="/Users/annabelbeichman/Documents/UW/Human_MUTYH/results/Unfazed/230912_minimal_phased_mutations_IGV_filtered/" # these have bugfix, minimal VAF, and have had IGV red calls excluded

totalCountsPhasingSuccess=read.table(paste0(indir,"total_counts_with_parents_IGV_filtered.txt"),header=T,sep="\t")
# this also has accessible bases 

totalCountsPhasingSuccess <- totalCountsPhasingSuccess %>%
  rename(fraction_phased = percent_phased) # fixing naming since this is a fraction not a percent
# read in children into one big table
phasedMutations_allKids <- data.frame()

children=unique(totalCountsPhasingSuccess[!totalCountsPhasingSuccess$child_ID %in% c("P1","P2","P3","P4"),]$child_ID) # need to exclude Probands b/c not phased
children
for(child in children){
  muts = read.table(paste0(indir,"/phased_mutations_per_individual/",child,"_mutations_IGV_filtered.txt"),header=T,sep="\t",comment.char = "",strip.white = T,blank.lines.skip = T)
  muts$child_id <- child
  phasedMutations_allKids <- bind_rows(phasedMutations_allKids,muts)
  
}


#### want paternal/maternal totals:
phasedMutations_allKids_SummedUp <- phasedMutations_allKids %>% 
  group_by(child_id,parent_origin) %>%
  summarise(count=n()) %>%
  ungroup()

phasedMutations_allKids_SummedUp_wide <- pivot_wider(phasedMutations_allKids_SummedUp,id_cols = "child_id",names_from = "parent_origin",values_from = "count")



phasedMutations_allKids_SummedUp_wide <- phasedMutations_allKids_SummedUp_wide %>%
  rename(count_paternally_phased=paternal,count_maternally_phased=maternal)

# fill in NAs with 0s
phasedMutations_allKids_SummedUp_wide[is.na(phasedMutations_allKids_SummedUp_wide)] <- 0

# get fraction phased to each parent: 
phasedMutations_allKids_SummedUp_wide <- phasedMutations_allKids_SummedUp_wide %>%
  mutate(fraction_maternally_phased=count_maternally_phased/(count_paternally_phased+count_maternally_phased),fraction_paternally_phased=count_paternally_phased/(count_maternally_phased+count_paternally_phased))

# merge with total info:
fullDF_allInfo <- merge(totalCountsPhasingSuccess[!totalCountsPhasingSuccess$child_ID %in% c("P1","P2","P3","P4"),],phasedMutations_allKids_SummedUp_wide,by.x="child_ID",by.y="child_id",all=T)

fullDF_allInfo

write.table(fullDF_allInfo,paste0(indir,"ALLCHILDREN.TotalCounts.PhasedToEachParent.WithPhasingSucessRate.txt"),quote=F,sep="\t",row.names = F)

############# write out spectrum info ############ 

phasedMutations_allKids_SummedUp_SPECTRUM <- phasedMutations_allKids %>% 
  group_by(child_id,parent_origin,mut_type) %>%
  summarise(count=n()) %>%
  ungroup()

write.table(phasedMutations_allKids_SummedUp_SPECTRUM,paste0(indir,"ALLCHILDREN.todate.SpectrumCounts.PhasedToEachParent.long.txt"),quote=F,sep="\t",row.names = F)

