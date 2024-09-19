########### Sigfit ##################
require(tidyverse)
require(reshape2)
require(dplyr)
require(ggplot2)
require(tidyr)
require(sigfit)
require(spgs)  # for reverse complements
require(gtools) # for mixedsort
require(RColorBrewer)
require(scales)
require(sigfit)

wd="/Users/annabelbeichman/Documents/UW/Human_MUTYH/results/signatureFitting/"
data("cosmic_signatures_v3.2") # this pulls from sigfit

### 3mer mutation spectrum from David
# need to exclude 'red' calls that were flagged by IGV and total up per person/per type
mutationTypesPerSite <- read.table(paste0(wd,"240918_subsetted_threemerSpectrum.csv"),sep=",",header=T) ## THIS IS NOW UPDATED TO BE DS corrected 20240822 and had bug fix on 20240918

head(mutationTypesPerSite)

############ sum up mutation types per sample across all sites, excluding red ######
mutationSpectrum <- mutationTypesPerSite %>%
  filter(Type!="red") %>%
  group_by(Sample,ms) %>%
  summarize(total_count=n()) %>%
  ungroup()

# dim(mutationTypesPerSite[mutationTypesPerSite$Type=="red",]) #261 are red


dim(mutationTypesPerSite) # 931 --> 606 after ds corretion --> 642 after P-bug fix 
sum(mutationSpectrum$total_count) #670 --> 606 after ds correction and candice removed red sites so 0 red sites are present now; cool --> 642 after bug fix

head(mutationTypesPerSite)

#### make a sigfit like mutation label. currently is ACT>G , but make it ACT>AGT ####
mutationSpectrum$mutation_label <- paste0(substr(mutationSpectrum$ms,1,4),substr(mutationSpectrum$ms,1,1),substr(mutationSpectrum$ms,5,5),substr(mutationSpectrum$ms,3,3))
# pasting together first 4 chars (XYZ>A) and then first char (X), then 5th char (A), then third char (Z)


######## format for SIGFIT ######
# make into a wide table, with rownames = sample names, colnames= mutation labels,  and reorder columns to match sigfit cosmic signatures
mutationSpectrum_wider <-  pivot_wider(mutationSpectrum,id_cols="Sample",names_from="mutation_label",values_from ="total_count")

head(mutationSpectrum_wider)

# need to fill in NA --> 0
mutationSpectrum_wider[is.na(mutationSpectrum_wider)] <- 0

head(mutationSpectrum_wider)


# add species labels as rowname (important! )
mutationSpectrum_wider <- mutationSpectrum_wider %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="Sample") # adds row names

head(mutationSpectrum_wider)


# some types are totally missing (0 in all samples)
# need to fill those in 
missingmutationtypes <- colnames(cosmic_signatures_v3.2)[!colnames(cosmic_signatures_v3.2) %in% colnames(mutationSpectrum_wider)]
missingmutationtypes
# "CCG>CGG" "ATG>AAG" "ATC>AGC" "GTA>GGA" "GTG>GGG" "GTT>GGT" "TTC>TGC"
# "CCG>CGG" "ATG>AAG" "ATC>AGC" "ATT>AGT" "GTA>GGA" "GTG>GGG" "GTT>GGT" "TTC>TGC" # after ds correction
# "CCG>CGG" "ATG>AAG" "ATC>AGC" "ATT>AGT" "GTA>GGA" "GTG>GGG" "GTT>GGT" "TTC>TGC" after bug fix
mutationSpectrum_wider[,missingmutationtypes] <- 0 # this adds those missing mutation types as empty columns at end of df

# reorder columns to match the order of the cosmic signatures
mutationSpectrum_wider <- mutationSpectrum_wider[,colnames(cosmic_signatures_v3.2)]

# 
write.table(mutationSpectrum_wider,paste0(wd,"sigfitFormat3merspectrum.txt"),sep="\t",row.names=T,quote=F)

########## make a mutation spectrum per family ###########
family1 <- c("C11","C12")
family2 <- c("C21","C22","C23")
family3 <- c("C31","C32")
family4 <- c("C41","C42")
probands <- c("P1","P3","P4")

mutationTypesPerSite$family <- ""
mutationTypesPerSite[mutationTypesPerSite$Sample %in% family1,]$family  <- "family I (bi-mom)"
mutationTypesPerSite[mutationTypesPerSite$Sample %in% family2,]$family  <- "family II (bi-mom)"
mutationTypesPerSite[mutationTypesPerSite$Sample %in% family3,]$family  <- "family III (bi-dad)"
mutationTypesPerSite[mutationTypesPerSite$Sample %in% family4,]$family  <- "family IV (mono-mom)"

mutationTypesPerSite[mutationTypesPerSite$Sample %in% probands,]$family  <- "probands (bi-mom-and-dad)"


mutationSpectrum_perFam <- mutationTypesPerSite %>%
  filter(Type!="red") %>%
  group_by(family,ms) %>%
  summarize(total_count=n()) %>%
  ungroup()

#### make a sigfit like mutation label. currently is ACT>G , but make it ACT>AGT ####
mutationSpectrum_perFam$mutation_label <- paste0(substr(mutationSpectrum_perFam$ms,1,4),substr(mutationSpectrum_perFam$ms,1,1),substr(mutationSpectrum_perFam$ms,5,5),substr(mutationSpectrum_perFam$ms,3,3))
# pasting together first 4 chars (XYZ>A) and then first char (X), then 5th char (A), then third char (Z)


######## format per-fam mutation spectrum for SIGFIT ######
# make into a wide table, with rownames = sample names, colnames= mutation labels,  and reorder columns to match sigfit cosmic signatures
mutationSpectrum_perFam_wider <-  pivot_wider(mutationSpectrum_perFam,id_cols="family",names_from="mutation_label",values_from ="total_count")

head(mutationSpectrum_perFam_wider)

# need to fill in NA --> 0
mutationSpectrum_perFam_wider[is.na(mutationSpectrum_perFam_wider)] <- 0

head(mutationSpectrum_perFam_wider)


# add species labels as rowname (important! )
mutationSpectrum_perFam_wider <- mutationSpectrum_perFam_wider %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="family") # adds row names

head(mutationSpectrum_perFam_wider)

data("cosmic_signatures_v3.2") # this pulls from sigfit

# some types are totally missing (0 in all samples)
# need to fill those in 
missingmutationtypes_perFam <- colnames(cosmic_signatures_v3.2)[!colnames(cosmic_signatures_v3.2) %in% colnames(mutationSpectrum_perFam_wider)]
# "CCG>CGG" "ATG>AAG" "ATC>AGC" "GTA>GGA" "GTG>GGG" "GTT>GGT" "TTC>TGC"
missingmutationtypes_perFam
mutationSpectrum_perFam_wider[,missingmutationtypes_perFam] <- 0 # this adds those missing mutation types as empty columns at end of df

# reorder columns to match the order of the cosmic signatures
mutationSpectrum_perFam_wider <- mutationSpectrum_perFam_wider[,colnames(cosmic_signatures_v3.2)]

# 
write.table(mutationSpectrum_perFam_wider,paste0(wd,"sigfitFormat3merspectrum.PERFAMILY.txt"),sep="\t",row.names=T,quote=F)

########## make a spectrum per group #############
# groups to sum over: 
biallelicMother=c("C11","C12","C23","C21","C22") # two fams combined
biallelicFather=c("C31","C32")
monoallelicMother=c("C41","C42")
monoallelic_BothParents=c("P1","P2","P3","P4")

mutationTypesPerSite$group <- ""
mutationTypesPerSite[mutationTypesPerSite$Sample %in% biallelicMother,]$group  <- "bi-allelic mother"
mutationTypesPerSite[mutationTypesPerSite$Sample %in% biallelicFather,]$group  <- "bi-allelic father"
mutationTypesPerSite[mutationTypesPerSite$Sample %in% monoallelicMother,]$group  <- "mono-allelic mother"
mutationTypesPerSite[mutationTypesPerSite$Sample %in% monoallelic_BothParents,]$group  <- "mono-allelic mother+father"


mutationSpectrum_perGroup <- mutationTypesPerSite %>%
  filter(Type!="red") %>%
  group_by(group,ms) %>%
  summarize(total_count=n()) %>%
  ungroup()

#### make a sigfit like mutation label. currently is ACT>G , but make it ACT>AGT ####
mutationSpectrum_perGroup$mutation_label <- paste0(substr(mutationSpectrum_perGroup$ms,1,4),substr(mutationSpectrum_perGroup$ms,1,1),substr(mutationSpectrum_perGroup$ms,5,5),substr(mutationSpectrum_perGroup$ms,3,3))
# pasting together first 4 chars (XYZ>A) and then first char (X), then 5th char (A), then third char (Z)


######## format per-group mutation spectrum for SIGFIT ######
# make into a wide table, with rownames = sample names, colnames= mutation labels,  and reorder columns to match sigfit cosmic signatures
mutationSpectrum_perGroup_wider <-  pivot_wider(mutationSpectrum_perGroup,id_cols="group",names_from="mutation_label",values_from ="total_count")

head(mutationSpectrum_perGroup_wider)

# need to fill in NA --> 0
mutationSpectrum_perGroup_wider[is.na(mutationSpectrum_perGroup_wider)] <- 0

head(mutationSpectrum_perGroup_wider)


# add species labels as rowname (important! )
mutationSpectrum_perGroup_wider <- mutationSpectrum_perGroup_wider %>%
  remove_rownames %>% # have to remove any preexisting rownames 
  column_to_rownames(var="group") # adds row names

head(mutationSpectrum_perGroup_wider)

data("cosmic_signatures_v3.2") # this pulls from sigfit

# some types are totally missing (0 in all samples)
# need to fill those in 
missingmutationtypes_perGroup <- colnames(cosmic_signatures_v3.2)[!colnames(cosmic_signatures_v3.2) %in% colnames(mutationSpectrum_perGroup_wider)]
# "CCG>CGG" "ATG>AAG" "ATC>AGC" "GTA>GGA" "GTG>GGG" "GTT>GGT" "TTC>TGC"
missingmutationtypes_perGroup
mutationSpectrum_perGroup_wider[,missingmutationtypes_perGroup] <- 0 # this adds those missing mutation types as empty columns at end of df

# reorder columns to match the order of the cosmic signatures
mutationSpectrum_perGroup_wider <- mutationSpectrum_perGroup_wider[,colnames(cosmic_signatures_v3.2)]

# 
write.table(mutationSpectrum_perGroup_wider,paste0(wd,"sigfitFormat3merspectrum.PERGROUP.txt"),sep="\t",row.names=T,quote=F)

######### also format for SPE ##########
# see this example file for how to format (must be exactly like this!!!)
spe_examplefile <- read.table("/Users/annabelbeichman/Documents/UW/Human_MUTYH/results/signatureFitting/sigprofilerextractor/exampleData.OrderYourDataThisWay.txt",header = T,sep = "\t")
head(spe_examplefile)


##### SPE: per individual ######
mutationSpectrum$mutation_label_SPEFormat <- paste0(substr(mutationSpectrum$ms,1,1),"[",substr(mutationSpectrum$ms,2,2),">",substr(mutationSpectrum$ms,5,5),"]",substr(mutationSpectrum$ms,3,3))

mutationSpectrum_SPE <- pivot_wider(mutationSpectrum,names_from=Sample,id_cols=mutation_label_SPEFormat,values_from = total_count )
head(mutationSpectrum_SPE)

# rename cols to match 'Mutation Types' col name
colnames(mutationSpectrum_SPE) <- c("Mutation Types",names(select(mutationSpectrum_SPE,-"mutation_label_SPEFormat")))


# need to add empty types:
missingTypes_spe <- data.frame(spe_examplefile$Mutation.Types[!spe_examplefile$Mutation.Types %in% mutationSpectrum_SPE$`Mutation Types`])
#  "A[T>A]G" "A[T>G]C" "C[C>G]G" "G[T>G]A" "G[T>G]G" "G[T>G]T" "T[T>G]C"
colnames(missingTypes_spe) <- "Mutation Types"
head(missingTypes_spe)
mutationSpectrum_SPE <- bind_rows(missingTypes_spe,mutationSpectrum_SPE)
head(mutationSpectrum_SPE)
# add 0s for nas
mutationSpectrum_SPE[is.na(mutationSpectrum_SPE)] <-0


# THIS Puts rows in same order:
mutationSpectrum_SPE <- mutationSpectrum_SPE[match(spe_examplefile$Mutation.Types, mutationSpectrum_SPE$`Mutation Types`),]

mutationSpectrum_SPE$`Mutation Types`==spe_examplefile$Mutation.Types # good. 

if(sum(mutationSpectrum_SPE$`Mutation Types`!=spe_examplefile$Mutation.Types)!=0){
  print("rows are in wrong order!")
}

write.table(mutationSpectrum_SPE,paste0(wd,"AllSpectra.FormattedForSigProfilerExtractor.PERINDIVIDUAL.txt"),row.names = F,quote=F,sep="\t")
# need to 





###### SPE: per family ######
mutationSpectrum_perFam$mutation_label_SPEFormat <- paste0(substr(mutationSpectrum_perFam$ms,1,1),"[",substr(mutationSpectrum_perFam$ms,2,2),">",substr(mutationSpectrum_perFam$ms,5,5),"]",substr(mutationSpectrum_perFam$ms,3,3))

mutationSpectrum_perFam_SPE <- pivot_wider(mutationSpectrum_perFam,names_from=family,id_cols=mutation_label_SPEFormat,values_from = total_count )
head(mutationSpectrum_perFam_SPE)

# rename cols to match 'Mutation Types' col name
colnames(mutationSpectrum_perFam_SPE) <- c("Mutation Types",names(select(mutationSpectrum_perFam_SPE,-"mutation_label_SPEFormat")))


# need to add empty types:
missingTypes_spe <- data.frame(spe_examplefile$Mutation.Types[!spe_examplefile$Mutation.Types %in% mutationSpectrum_perFam_SPE$`Mutation Types`])
#  "A[T>A]G" "A[T>G]C" "C[C>G]G" "G[T>G]A" "G[T>G]G" "G[T>G]T" "T[T>G]C"
colnames(missingTypes_spe) <- "Mutation Types"
head(missingTypes_spe)
mutationSpectrum_perFam_SPE <- bind_rows(missingTypes_spe,mutationSpectrum_perFam_SPE)
head(mutationSpectrum_perFam_SPE)
# add 0s for nas
mutationSpectrum_perFam_SPE[is.na(mutationSpectrum_perFam_SPE)] <-0


# THIS Puts rows in same order:
mutationSpectrum_perFam_SPE <- mutationSpectrum_perFam_SPE[match(spe_examplefile$Mutation.Types, mutationSpectrum_perFam_SPE$`Mutation Types`),]

mutationSpectrum_perFam_SPE$`Mutation Types`==spe_examplefile$Mutation.Types # good. 

if(sum(mutationSpectrum_perFam_SPE$`Mutation Types`!=spe_examplefile$Mutation.Types)!=0){
  print("rows are in wrong order!")
}

write.table(mutationSpectrum_perFam_SPE,paste0(wd,"AllSpectra.FormattedForSigProfilerExtractor.PERFAMILY.txt"),row.names = F,quote=F,sep="\t")



#### SPE per group #####
mutationSpectrum_perGroup$mutation_label_SPEFormat <- paste0(substr(mutationSpectrum_perGroup$ms,1,1),"[",substr(mutationSpectrum_perGroup$ms,2,2),">",substr(mutationSpectrum_perGroup$ms,5,5),"]",substr(mutationSpectrum_perGroup$ms,3,3))

mutationSpectrum_perGroup_SPE <- pivot_wider(mutationSpectrum_perGroup,names_from=group,id_cols=mutation_label_SPEFormat,values_from = total_count )
head(mutationSpectrum_perGroup_SPE)

# rename cols to match 'Mutation Types' col name
colnames(mutationSpectrum_perGroup_SPE) <- c("Mutation Types",names(select(mutationSpectrum_perGroup_SPE,-"mutation_label_SPEFormat")))


# need to add empty types:
missingTypes_spe <- data.frame(spe_examplefile$Mutation.Types[!spe_examplefile$Mutation.Types %in% mutationSpectrum_perGroup_SPE$`Mutation Types`])
#  "A[T>A]G" "A[T>G]C" "C[C>G]G" "G[T>G]A" "G[T>G]G" "G[T>G]T" "T[T>G]C"
colnames(missingTypes_spe) <- "Mutation Types"
head(missingTypes_spe)
mutationSpectrum_perGroup_SPE <- bind_rows(missingTypes_spe,mutationSpectrum_perGroup_SPE)
head(mutationSpectrum_perGroup_SPE)
# add 0s for nas
mutationSpectrum_perGroup_SPE[is.na(mutationSpectrum_perGroup_SPE)] <-0


# THIS Puts rows in same order:
mutationSpectrum_perGroup_SPE <- mutationSpectrum_perGroup_SPE[match(spe_examplefile$Mutation.Types, mutationSpectrum_perGroup_SPE$`Mutation Types`),]

mutationSpectrum_perGroup_SPE$`Mutation Types`==spe_examplefile$Mutation.Types # good. 

if(sum(mutationSpectrum_perGroup_SPE$`Mutation Types`!=spe_examplefile$Mutation.Types)!=0){
  print("rows are in wrong order!")
}

write.table(mutationSpectrum_perGroup_SPE,paste0(wd,"AllSpectra.FormattedForSigProfilerExtractor.PERGROUP.txt"),row.names = F,quote=F,sep="\t")

