#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KMA Sockeye Mixtures 2014-2016 ####
# Kyle Shedd Wed Mar 30 08:53:54 2016
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#### Introduction ####
# The goal of this script is to analyze sockeye salmon mixtures from the KMA
# commercial harvest from 2014-2016 using a coastwide sockeye baseline
# containing 473 populations in 15 reporting groups characterized by 48 SNPs
# grouped into 46 loci (2 sets of linked SNPs are combined as haplotypes).
# All mixtures are to be analyzed with the program BAYES.

#### Specific Objectives ####
# This script will:
# 1) Import mixture data
# 2) Define spatio-temporal strata
# 3) Perform a data QC on mixtures
# 4) Prepare BAYES input files
# 5) Summarize BAYES results
# 6) Generate plots and tables of results

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
username <- "krshedd"
password <- "********"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Repeat Late August ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Initially analyzed Late August harvest for 2014 from Uyak, Uganik, and Karluk
# on the WASSIP baseline with different reporting groups (EASSIP).

## Get collection SILLYs
LateAugustMixtures2014 <- dget(file = "Objects/EASSIP-LateAugust2014/LateAugustMixtures2014.txt")

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye2011_96SNPs", username = username, password = password)

## Save original LocusControl
loci96 <- LocusControl$locusnames
mito.loci <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl.txt")
dput(x = loci96, file = "Objects/loci96.txt")
dput(x = mito.loci, file = "Objects/mito.loci.txt")

#~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = LateAugustMixtures2014, username = username, password = password)
rm(username, password)
objects(pattern = "\\.gcl")

## Save unaltered .gcl's as back-up:
dir.create("Raw genotypes/OriginalCollections")
dir.create("Raw genotypes/OriginalCollections/KMALateAugust2014")
invisible(sapply(LateAugustMixtures2014, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/KMALateAugust2014/" , silly, ".txt", sep = ''))} )); beep(8)

## Original sample sizes by SILLY
collection.size.original <- sapply(LateAugustMixtures2014, function(silly) get(paste(silly, ".gcl", sep = ""))$n)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Confirming samples sizes
sapply(LateAugustMixtures2014, function(silly) {table(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE)} )


## SKARLC
str(SKARLC14.gcl$attributes$CAPTURE_DATE)
sort(unique(SKARLC14.gcl$attributes$CAPTURE_DATE))
table(SKARLC14.gcl$attributes$CAPTURE_DATE)

SKARLC_Aug26IDs <- AttributesToIDs.GCL(silly = "SKARLC14", attribute = "CAPTURE_DATE", matching = sort(unique(SKARLC14.gcl$attributes$CAPTURE_DATE))[11])

SKARLC_Aug26IDs <- list(as.numeric(na.omit(SKARLC_Aug26IDs)))
names(SKARLC_Aug26IDs) <- "SKARLC14"

PoolCollections.GCL("SKARLC14", loci = loci96, IDs = SKARLC_Aug26IDs, newname = "SKARLC14_Aug26")
SKARLC_Aug26.gcl$n  # 285


## SUGANC
str(SUGANC14.gcl$attributes$CAPTURE_DATE)
sort(unique(SUGANC14.gcl$attributes$CAPTURE_DATE))
table(SUGANC14.gcl$attributes$CAPTURE_DATE)

SUGANC_Aug27IDs <- AttributesToIDs.GCL(silly = "SUGANC14", attribute = "CAPTURE_DATE", matching = sort(unique(SUGANC14.gcl$attributes$CAPTURE_DATE))[6])

SUGANC_Aug27IDs <- list(as.numeric(na.omit(SUGANC_Aug27IDs)))
names(SUGANC_Aug27IDs) <- "SUGANC14"

PoolCollections.GCL("SUGANC14", loci = loci96, IDs = SUGANC_Aug27IDs, newname = "SUGANC14_Aug27")
SUGANC_Aug27.gcl$n  # 285


## SUYAKC
str(SUYAKC14.gcl$attributes$CAPTURE_DATE)
sort(unique(SUYAKC14.gcl$attributes$CAPTURE_DATE))
table(SUYAKC14.gcl$attributes$CAPTURE_DATE)

SUYAKC_Aug26IDs <- AttributesToIDs.GCL(silly = "SUYAKC14", attribute = "CAPTURE_DATE", matching = sort(unique(SUYAKC14.gcl$attributes$CAPTURE_DATE))[22])

SUYAKC_Aug26IDs <- list(as.numeric(na.omit(SUYAKC_Aug26IDs)))
names(SUYAKC_Aug26IDs) <- "SUYAKC14"

PoolCollections.GCL("SUYAKC14", loci = loci96, IDs = SUYAKC_Aug26IDs, newname = "SUYAKC14_Aug26")
SUYAKC_Aug26.gcl$n  # 285

LateAugustMixtures2014Strata <- c("SKARLC14_Aug26", "SUGANC14_Aug27", "SUYAKC14_Aug26")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

LateAugustMixtures

LateAugustMixtures_SampleSizes <- matrix(data = NA, nrow = length(LateAugustMixtures), ncol = 5, 
                                         dimnames = list(LateAugustMixtures, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_LateAugustMixtures_SampleSizebyLocus <- SampSizeByLocus.GCL(LateAugustMixtures, loci96)
min(Original_LateAugustMixtures_SampleSizebyLocus) ## 267/285
apply(Original_LateAugustMixtures_SampleSizebyLocus,1,min)/apply(Original_LateAugustMixtures_SampleSizebyLocus,1,max) # Good, 0.947


#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_LateAugustMixtures_ColSize <- sapply(paste(LateAugustMixtures, ".gcl", sep = ''), function(x) get(x)$n)
LateAugustMixtures_SampleSizes[, "Genotyped"] <- Original_LateAugustMixtures_ColSize


### Alternate
## Indentify alternate species individuals
LateAugustMixtures_Alternate <- FindAlternateSpecies.GCL(sillyvec = LateAugustMixtures, species = "sockeye")

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport = LateAugustMixtures_Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)

## Get number of individuals per silly after removing alternate species individuals
ColSize_LateAugustMixtures_PostAlternate <- sapply(paste(LateAugustMixtures, ".gcl", sep = ''), function(x) get(x)$n)
LateAugustMixtures_SampleSizes[, "Alternate"] <- Original_LateAugustMixtures_ColSize-ColSize_LateAugustMixtures_PostAlternate


### Missing
## Remove individuals with >20% missing data
LateAugustMixtures_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = LateAugustMixtures, loci = loci96, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_LateAugustMixtures_PostMissLoci <- sapply(paste(LateAugustMixtures, ".gcl", sep = ''), function(x) get(x)$n)
LateAugustMixtures_SampleSizes[, "Missing"] <- ColSize_LateAugustMixtures_PostAlternate-ColSize_LateAugustMixtures_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
LateAugustMixtures_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = LateAugustMixtures, loci = loci96, quantile = NULL, minproportion = 0.95)
LateAugustMixtures_DuplicateCheckReportSummary <- sapply(LateAugustMixtures, function(x) LateAugustMixtures_DuplicateCheck95MinProportion[[x]]$report)

## Remove duplicate individuals
LateAugustMixtures_RemovedDups <- RemoveDups.GCL(LateAugustMixtures_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_LateAugustMixtures_PostDuplicate <- sapply(paste(LateAugustMixtures, ".gcl", sep = ''), function(x) get(x)$n)
LateAugustMixtures_SampleSizes[, "Duplicate"] <- ColSize_LateAugustMixtures_PostMissLoci-ColSize_LateAugustMixtures_PostDuplicate


### Final
LateAugustMixtures_SampleSizes[, "Final"] <- ColSize_LateAugustMixtures_PostDuplicate
LateAugustMixtures_SampleSizes

write.xlsx(LateAugustMixtures_SampleSizes, file = "Output/LateAugustMixtures_SampleSizes.xlsx")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Combine Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~