#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### KMA Sockeye Mixtures 2014-2016 ####
# Kyle Shedd Wed Mar 30 08:53:54 2016
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Introduction ####
# The goal of this script is to analyze sockeye salmon mixtures from the KMA
# commercial harvest from 2014-2016 using a coastwide sockeye baseline
# containing 473 populations in 15 reporting groups characterized by 48 SNPs
# grouped into 46 loci (2 sets of linked SNPs are combined as haplotypes).
# All mixtures are to be analyzed with the program BAYES.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
LateAugustMixtures2014 <- dget(file = "Objects/EASSIP-LateAugust2014/LateAugustMixtures2014.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("EASSIP-LateAugust2014", "OriginalLocusControl.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered mixtures
invisible(sapply(LateAugustMixtures2014, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections/KMALateAugust2014/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Confirming samples sizes by date
sapply(LateAugustMixtures2014, function(silly) {table(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE)} )


## SKARLC
str(SKARLC14.gcl$attributes$CAPTURE_DATE)
sort(unique(SKARLC14.gcl$attributes$CAPTURE_DATE))
max(sort(unique(SKARLC14.gcl$attributes$CAPTURE_DATE)))
table(SKARLC14.gcl$attributes$CAPTURE_DATE)

SKARLC14_Aug26IDs <- AttributesToIDs.GCL(silly = "SKARLC14", attribute = "CAPTURE_DATE", matching = max(sort(unique(SKARLC14.gcl$attributes$CAPTURE_DATE))))

SKARLC14_Aug26IDs <- list(as.numeric(na.omit(SKARLC14_Aug26IDs)))
names(SKARLC14_Aug26IDs) <- "SKARLC14"

PoolCollections.GCL(collections = "SKARLC14", loci = loci96, IDs = SKARLC14_Aug26IDs, newname = "SKARLC14_Aug26")
SKARLC14_Aug26.gcl$n  # 285


## SUGANC
str(SUGANC14.gcl$attributes$CAPTURE_DATE)
sort(unique(SUGANC14.gcl$attributes$CAPTURE_DATE))
max(sort(unique(SUGANC14.gcl$attributes$CAPTURE_DATE)))
table(SUGANC14.gcl$attributes$CAPTURE_DATE)

SUGANC14_Aug27IDs <- AttributesToIDs.GCL(silly = "SUGANC14", attribute = "CAPTURE_DATE", matching = max(sort(unique(SUGANC14.gcl$attributes$CAPTURE_DATE))))

SUGANC14_Aug27IDs <- list(as.numeric(na.omit(SUGANC14_Aug27IDs)))
names(SUGANC14_Aug27IDs) <- "SUGANC14"

PoolCollections.GCL(collections = "SUGANC14", loci = loci96, IDs = SUGANC14_Aug27IDs, newname = "SUGANC14_Aug27")
SUGANC14_Aug27.gcl$n  # 285


## SUYAKC
str(SUYAKC14.gcl$attributes$CAPTURE_DATE)
sort(unique(SUYAKC14.gcl$attributes$CAPTURE_DATE))
max(sort(unique(SUYAKC14.gcl$attributes$CAPTURE_DATE)))
table(SUYAKC14.gcl$attributes$CAPTURE_DATE)

SUYAKC14_Aug26IDs <- AttributesToIDs.GCL(silly = "SUYAKC14", attribute = "CAPTURE_DATE", matching = max(sort(unique(SUYAKC14.gcl$attributes$CAPTURE_DATE))))

SUYAKC14_Aug26IDs <- list(as.numeric(na.omit(SUYAKC14_Aug26IDs)))
names(SUYAKC14_Aug26IDs) <- "SUYAKC14"

PoolCollections.GCL(collections = "SUYAKC14", loci = loci96, IDs = SUYAKC14_Aug26IDs, newname = "SUYAKC14_Aug26")
SUYAKC14_Aug26.gcl$n  # 285

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mixture sillyvec
LateAugustMixtures2014Strata <- c("SKARLC14_Aug26", "SUGANC14_Aug27", "SUYAKC14_Aug26")
dput(x = LateAugustMixtures2014Strata, file = "Objects/LateAugustMixtures2014Strata.txt")

# Confirm sample sizes
sapply(LateAugustMixtures2014Strata, function(silly) get(paste(silly, ".gcl", sep = ""))$n)

# dput mixture sillys
dir.create("Raw genotypes/OriginalCollections_Strata")
invisible(sapply(LateAugustMixtures2014Strata, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_Strata/" , silly, ".txt", sep = ''))} )); beep(8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

LateAugustMixtures2014Strata

LateAugustMixtures2014Strata_SampleSizes <- matrix(data = NA, nrow = length(LateAugustMixtures2014Strata), ncol = 5, 
                                                   dimnames = list(LateAugustMixtures2014Strata, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_LateAugustMixtures2014Strata_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = LateAugustMixtures2014Strata, loci = loci96)
min(Original_LateAugustMixtures2014Strata_SampleSizebyLocus)  ## 267/285
apply(Original_LateAugustMixtures2014Strata_SampleSizebyLocus, 1, min) / apply(Original_LateAugustMixtures2014Strata_SampleSizebyLocus, 1, max) # Good, 0.947


#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_LateAugustMixtures2014Strata_ColSize <- sapply(paste(LateAugustMixtures2014Strata, ".gcl", sep = ''), function(x) get(x)$n)
LateAugustMixtures2014Strata_SampleSizes[, "Genotyped"] <- Original_LateAugustMixtures2014Strata_ColSize


### Alternate
## Indentify alternate species individuals
LateAugustMixtures2014Strata_Alternate <- FindAlternateSpecies.GCL(sillyvec = LateAugustMixtures2014Strata, species = "sockeye")

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport = LateAugustMixtures2014Strata_Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)

## Get number of individuals per silly after removing alternate species individuals
ColSize_LateAugustMixtures2014Strata_PostAlternate <- sapply(paste(LateAugustMixtures2014Strata, ".gcl", sep = ''), function(x) get(x)$n)
LateAugustMixtures2014Strata_SampleSizes[, "Alternate"] <- Original_LateAugustMixtures2014Strata_ColSize-ColSize_LateAugustMixtures2014Strata_PostAlternate


### Missing
## Remove individuals with >20% missing data
LateAugustMixtures2014Strata_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = LateAugustMixtures2014Strata, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_LateAugustMixtures2014Strata_PostMissLoci <- sapply(paste(LateAugustMixtures2014Strata, ".gcl", sep = ''), function(x) get(x)$n)
LateAugustMixtures2014Strata_SampleSizes[, "Missing"] <- ColSize_LateAugustMixtures2014Strata_PostAlternate-ColSize_LateAugustMixtures2014Strata_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
LateAugustMixtures2014Strata_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = LateAugustMixtures2014Strata, loci = loci96, quantile = NULL, minproportion = 0.95)
LateAugustMixtures2014Strata_DuplicateCheckReportSummary <- sapply(LateAugustMixtures2014Strata, function(x) LateAugustMixtures2014Strata_DuplicateCheck95MinProportion[[x]]$report)

## Remove duplicate individuals
LateAugustMixtures2014Strata_RemovedDups <- RemoveDups.GCL(LateAugustMixtures2014Strata_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_LateAugustMixtures2014Strata_PostDuplicate <- sapply(paste(LateAugustMixtures2014Strata, ".gcl", sep = ''), function(x) get(x)$n)
LateAugustMixtures2014Strata_SampleSizes[, "Duplicate"] <- ColSize_LateAugustMixtures2014Strata_PostMissLoci-ColSize_LateAugustMixtures2014Strata_PostDuplicate


### Final
LateAugustMixtures2014Strata_SampleSizes[, "Final"] <- ColSize_LateAugustMixtures2014Strata_PostDuplicate
LateAugustMixtures2014Strata_SampleSizes

write.xlsx(LateAugustMixtures2014Strata_SampleSizes, file = "Output/LateAugustMixtures2014Strata_SampleSizes.xlsx")
dput(x = LateAugustMixtures2014Strata_SampleSizes, file = "Objects/LateAugustMixtures2014Strata_SampleSizes.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Combine Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### loci46
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get the final locus set from baseline file to avoid any errors
loci46 <- dget(file = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects/loci46.txt")
dput(x = loci46, file = "Objects/loci46.txt")

## Combine loci
combined.loci.46 <- sapply(grep(pattern = "\\.", x = loci46, value = TRUE), function(locus) {unlist(strsplit(x = locus, split = "\\."))}, simplify = FALSE)
sapply(combined.loci.46, function(loci2combine) {CombineLoci.GCL(sillyvec = LateAugustMixtures2014Strata, markerset = loci2combine, update = TRUE)} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### loci89
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get the final locus set from baseline file to avoid any errors
loci89 <- dget(file = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects/loci89.txt")
dput(x = loci89, file = "Objects/loci89.txt")

## Combine loci
combined.loci.89 <- sapply(grep(pattern = "\\.", x = loci89, value = TRUE), function(locus) {unlist(strsplit(x = locus, split = "\\."))}, simplify = FALSE)
# Which haven't already been combined?
combined.loci89.not46 <- sapply(grep(pattern = "\\.", x = loci89, value = TRUE)[!grep(pattern = "\\.", x = loci89, value = TRUE) %in% grep(pattern = "\\.", x = loci46, value = TRUE)], function(locus) {unlist(strsplit(x = locus, split = "\\."))}, simplify = FALSE)
sapply(combined.loci89.not46, function(loci2combine) {CombineLoci.GCL(sillyvec = LateAugustMixtures2014Strata, markerset = loci2combine, update = TRUE)} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save PostQC/Combined loci .gcl's as back-up:
dir.create("Raw genotypes/OriginalCollections_Strata_PostQC_CombinedLoci")
invisible(sapply(LateAugustMixtures2014Strata, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_Strata_PostQC_CombinedLoci/" , silly, ".txt", sep = ''))} )); beep(8)

## Dput LocusControl
dput(x = LocusControl, file = "Objects/LocusControl99.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = LateAugustMixtures2014Strata, loci = loci96[-mito.loci], path = "Genepop/LateAugustMixtures2014Strata_93nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/LateAugustMixtures2014Strata_93nuclearloci.txt.P")

#~~~~~~~~~~~~~~~~~~
# Plot Fis values (looking for gross excess of hets [-Fis]) and look for markers out of HWE

# by(data = HWE$DataByPop, INDICES = HWE$DataByPop$Pop, FUN = function(x) {
#   plot(sort(x[, "WC Fis"]), type = "h", lwd = 5, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40", main = "Silly")
#   abline(h = 0, lwd = 5)
# } )

sapply(as.character(unique(HWE$DataByPop$Pop)), function(pop) {
  x <- subset(x = HWE$DataByPop, subset = Pop == pop)
  plot(sort(x[, "WC Fis"]), type = "h", lwd = 5, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40", main = pop)
  abline(h = 0, lwd = 5)
  x[x$PValue[!is.na(x$PValue)] < 0.05, ]
}, USE.NAMES = TRUE, simplify = FALSE
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/LocusControl99.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("EASSIP-LateAugust2014", "OriginalLocusControl.txt", "LocusControl98.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered mixtures
invisible(sapply(LateAugustMixtures2014Strata, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_Strata_PostQC_CombinedLoci/", silly, ".txt", sep = "")), pos = 1)} )); beep(4)
objects(pattern = "\\.gcl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get baseline objects needed for MSA
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects")

KMA473Pops15FlatPrior <- dget(file = "KMA473Pops15FlatPrior.txt")
KMA473PopsInits <- dget(file = "KMA473PopsInits.txt")
KMA473PopsGroupVec15 <- dget(file = "KMA473PopsGroupVec15.txt")
KMA473Pops <- dget(file = "KMA473Pops.txt")
KMA15GroupsPC <- dget(file = "PCGroups15.txt")
KMA47346Baseline <- dget(file = "KMA47346Baseline.txt")
KMA47389Baseline <- dget(file = "KMA47389Baseline.txt")
WASSIPSockeyeSeeds <- dget("V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt")
KMA473CommonNames <- dget(file = "CommonNames473.txt")

setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MSA files for BAYES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 89 Loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping Mixture files
KMA47389MixtureFormat <- CreateMixture.GCL(sillys = LateAugustMixtures2014Strata[1], loci = loci89, IDs = NULL, mixname = LateAugustMixtures2014Strata[1],
                                           dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(KMA47389MixtureFormat, file = "Objects/KMA47389MixtureFormat.txt")

sapply(LateAugustMixtures2014Strata, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci89, IDs = NULL, mixname = Mix, dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)})

## Dumping Control files
sapply(LateAugustMixtures2014Strata, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci89, mixname = Mix, basename = "KMA473Pops89Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits, dir = "BAYES/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47389MixtureFormat, basefortran = KMA47389Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(LateAugustMixtures2014Strata, function(Mix) {dir.create(paste("BAYES/Output/", Mix, sep = ""))})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 46 Loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping Mixture files
KMA47346MixtureFormat <- CreateMixture.GCL(sillys = LateAugustMixtures2014Strata[1], loci = loci46, IDs = NULL, mixname = LateAugustMixtures2014Strata[1],
                                           dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)
dput(KMA47346MixtureFormat, file = "Objects/KMA47346MixtureFormat.txt")

sapply(LateAugustMixtures2014Strata, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/Mixture", type = "BAYES", PT = FALSE)})

## Dumping Control files
sapply(LateAugustMixtures2014Strata, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits, dir = "BAYES/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(LateAugustMixtures2014Strata, function(Mix) {dir.create(paste("BAYES/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~