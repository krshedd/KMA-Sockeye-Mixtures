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
KMA15Colors <- dget(file = "Colors15.txt")
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")

## Copy these baseline objects and put them in the Mixtures/Objects directory

setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects")
file.copy(from = c("KMA473Pops15FlatPrior.txt", "KMA473PopsInits.txt", "KMA473PopsGroupVec15.txt", "KMA473Pops.txt", "PCGroups15.txt", "CommonNames473.txt", "Colors15.txt",
                   "V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt"), 
          to = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects")
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")


KMA15GroupsPC2Rows <- PCGroups15Rows2 <- c("West of\nChignik", "Black\nLake", "Chignik\nLake", "U. Station\nAkalura", "Frazer\n", "Ayakulik\n", "Karluk\n", "Uganik\n", "Northwest\nKodiak", "Afognak\n", "Eastside\nKodiak", "Saltery\n", "Cook\nInlet", "PWS\n", "South of\nCape Suckling")
dput(x = KMA15GroupsPC2Rows, file = "Objects/KMA15GroupsPC2Rows.txt")
# Note, I changed the names for CommonNames473.txt -> KMA473CommonNames.txt & PCGroups15 -> KMA15GroupsPC.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MSA files for BAYES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 89 Loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping Mixture files
dir.create("BAYES/Late August 89loci")
KMA47389MixtureFormat <- CreateMixture.GCL(sillys = LateAugustMixtures2014Strata[1], loci = loci89, IDs = NULL, mixname = LateAugustMixtures2014Strata[1],
                                           dir = "BAYES/Late August 89loci/Mixture", type = "BAYES", PT = FALSE)
dput(KMA47389MixtureFormat, file = "Objects/KMA47389MixtureFormat.txt")

sapply(LateAugustMixtures2014Strata, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci89, IDs = NULL, mixname = Mix, dir = "BAYES/Late August 89loci/Mixture", type = "BAYES", PT = FALSE)})

## Dumping Control files
sapply(LateAugustMixtures2014Strata, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci89, mixname = Mix, basename = "KMA473Pops89Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits, dir = "BAYES/Late August 89loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47389MixtureFormat, basefortran = KMA47389Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(LateAugustMixtures2014Strata, function(Mix) {dir.create(paste("BAYES/Late August 89loci/Output/", Mix, sep = ""))})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 46 Loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping Mixture files
KMA47346MixtureFormat <- CreateMixture.GCL(sillys = LateAugustMixtures2014Strata[1], loci = loci46, IDs = NULL, mixname = LateAugustMixtures2014Strata[1],
                                           dir = "BAYES/Late August 46loci/Mixture", type = "BAYES", PT = FALSE)
dput(KMA47346MixtureFormat, file = "Objects/KMA47346MixtureFormat.txt")

sapply(LateAugustMixtures2014Strata, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/Late August 46loci/Mixture", type = "BAYES", PT = FALSE)})

## Dumping Control files
sapply(LateAugustMixtures2014Strata, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits, dir = "BAYES/Late August 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(LateAugustMixtures2014Strata, function(Mix) {dir.create(paste("BAYES/Late August 46loci/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 89 Loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LateAugust2014Strata_89loci_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                                      maindir = "BAYES/Late August 89loci/Output", 
                                                                      mixvec = LateAugustMixtures2014Strata, prior = "",  
                                                                      ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(LateAugust2014Strata_89loci_Estimates, file = "Estimates objects/LateAugust2014Strata_89loci_Estimates.txt")
dput(LateAugust2014Strata_89loci_Estimates$Stats, file = "Estimates objects/LateAugust2014Strata_89loci_EstimatesStats.txt")

LateAugust2014Strata_89loci_Estimates <- dget(file = "Estimates objects/LateAugust2014Strata_89loci_Estimates.txt")
LateAugust2014Strata_89loci_EstimatesStats <- dget(file = "Estimates objects/LateAugust2014Strata_89loci_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(LateAugust2014Strata_89loci_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(LateAugust2014Strata_89loci_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})

# Quick look at raw posterior output
str(LateAugust2014Strata_89loci_Estimates$Output)
LateAugust2014StrataHeader <- c("Karluk Section August 26, 2014",
                                "Uganik Section August 27, 2014",
                                "Uyak Section August 26, 2014")
PlotPosterior <- function(output, header, groups, colors = NULL, set.mfrow, thin = 10){
  if(is.null(colors)) {colors <- rep("black", length(groups))}
  
  par(mfrow = set.mfrow, mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 1.1))
  
  thin = 10  # thins the posterior to create plots faster
  invisible(sapply(seq(length(output)), function(j) {
    invisible(sapply(seq(groups), function(i) {
      RG <- output[[j]][, i]
      RG <- RG[seq(1, length(RG), thin)]
      plot(RG, type = "l", ylim = c(0,1), xlab = "", ylab = "")
      abline(v = seq(0, length(RG), length(RG)/5))
      text(x = length(RG)/2, y = 0.98, labels = groups[i], col = colors[i], cex = 1.2, font = 2)} ))
    mtext(text = LateAugust2014StrataHeader[j], side = 3, outer = TRUE, cex = 1.5)
    mtext(text = "Iteration (5 chains each)", side = 1, outer = TRUE, cex = 1.5, line = 1.5)
    mtext(text = "Posterior", side = 2, outer = TRUE, cex = 1.5, line = 1.5)
  }))
  
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
}

PlotPosterior(output = LateAugust2014Strata_89loci_Estimates$Output, groups = KMA15GroupsPC, colors = KMA15Colors, header = LateAugust2014StrataHeader, set.mfrow = c(5, 3), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 46 Loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LateAugust2014Strata_46loci_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                                      maindir = "BAYES/Late August 46loci/Output", 
                                                                      mixvec = LateAugustMixtures2014Strata, prior = "",  
                                                                      ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(LateAugust2014Strata_46loci_Estimates, file = "Estimates objects/LateAugust2014Strata_46loci_Estimates.txt")
dput(LateAugust2014Strata_46loci_Estimates$Stats, file = "Estimates objects/LateAugust2014Strata_46loci_EstimatesStats.txt")

LateAugust2014Strata_46loci_Estimates <- dget(file = "Estimates objects/LateAugust2014Strata_46loci_Estimates.txt")
LateAugust2014Strata_46loci_EstimatesStats <- dget(file = "Estimates objects/LateAugust2014Strata_46loci_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(LateAugust2014Strata_46loci_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(LateAugust2014Strata_46loci_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})

# Quick look at raw posterior output
str(LateAugust2014Strata_46loci_Estimates$Output)
LateAugust2014StrataHeader <- c("Karluk Section August 26, 2014",
                                "Uganik Section August 27, 2014",
                                "Uyak Section August 26, 2014")

PlotPosterior(output = LateAugust2014Strata_46loci_Estimates$Output, groups = KMA15GroupsPC, colors = KMA15Colors, header = LateAugust2014StrataHeader, set.mfrow = c(5, 3), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
par(mar = c(5.6, 4.6, 3.1, 2.1))
sapply(LateAugustMixtures2014Strata, function(Mix) {
  plot(LateAugust2014Strata_46loci_Estimates$Stats[[Mix]][, "median"], cex = 3, pch = 16, col = KMA15Colors, ylab = "Proportion of Mixture", 
       ylim = c(0, 1), xlab = "", axes = FALSE, main = Mix, cex.main = 2, cex.lab = 1.5)
  arrows(x0 = seq(KMA15GroupsPC), y0 = LateAugust2014Strata_46loci_Estimates$Stats[[Mix]][, "5%"], x1 = seq(KMA15GroupsPC), 
         y1 = LateAugust2014Strata_46loci_Estimates$Stats[[Mix]][, "95%"], angle = 90, code = 3, length = 0.2, lwd = 2)
  points(LateAugust2014Strata_46loci_Estimates$Stats[[Mix]][, "median"], cex = 3, pch = 16, col = KMA15Colors)
  axis(side = 2)
  axis(side = 1, labels = NA, at = seq(KMA15GroupsPC))
  text(x = (seq(KMA15GroupsPC)) - 0.35, y = 0, labels = KMA15GroupsPC2Rows, srt = 60, pos = 1, offset = 3.5, xpd = TRUE)
})
par(mar = c(5.1, 4.1, 4.1, 2.1))

## Make violin plots of posteriors with RGs sorted

LateAugustMixtures2014StrataDetail <- c("Outer Karluk (255-20)\n8/26/14 (n = 284)", "Uganik (253-11, 13)\n8/27/14 (n = 282)", "Uyak (254-10, 20)\n8/26/14 (n = 280)")
names(LateAugustMixtures2014StrataDetail) <- LateAugustMixtures2014Strata

ViolinPlot <- function(estimates, groups, colors, header, wex = 1) {
  while(!require(vioplot, quietly = TRUE)){install.packages("vioplot")}
  
  mixtures <- names(estimates$Stats)
  
  par(mar = c(5.1, 4.6, 3.6, 1.1))
  sapply(mixtures, function(Mix) {
    plot(estimates$Stats[[Mix]][, "median"], cex = 3, pch = 16, col = colors, ylab = "Proportion of Mixture", ylim = c(0, 1), xlab = "", axes = FALSE, main = header[[Mix]], cex.main = 2, cex.lab = 1.5)
    sapply(seq(groups), function(i) {vioplot2(estimates$Output[[Mix]][, i], at = i, horizontal = FALSE, col = colors[i], border = TRUE, drawRect = FALSE, rectCol = colors[i], add = TRUE, wex = wex, lwd = 2)})
    arrows(x0 = seq(groups), y0 = estimates$Stats[[Mix]][, "5%"], x1 = seq(groups), y1 = estimates$Stats[[Mix]][, "95%"], angle = 90, code = 3, length = 0.2, lwd = 2)
    points(estimates$Stats[[Mix]][, "median"], cex = 2, pch = 21, col = "white", bg = colors, lwd = 3)
    axis(side = 2, lwd = 3, cex.axis = 1.5)
    text(x = (seq(groups)) - 0.35, y = 0, labels = groups, srt = 60, pos = 1, offset = 2.5, xpd = TRUE)
    axis(side = 1, labels = NA, at = seq(groups), pos = 0, lwd = 2, tick = FALSE)
    abline(h = 0, lwd = 3)
  } )
  
}

ViolinPlot(estimates = LateAugust2014Strata_89loci_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = LateAugustMixtures2014StrataDetail)
ViolinPlot(estimates = LateAugust2014Strata_46loci_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = LateAugustMixtures2014StrataDetail)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare markersets, side by side barplots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
require(devEMF)
require(gplots)
require(plotrix)


LateAugust2014Strata_89loci_EstimatesStats <- dget(file = "Estimates objects/LateAugust2014Strata_89loci_EstimatesStats.txt")
LateAugust2014Strata_46loci_EstimatesStats <- dget(file = "Estimates objects/LateAugust2014Strata_46loci_EstimatesStats.txt")

KMA15GroupsPC
KMA15GroupsPC2Rows
LateAugustMixtures2014StrataDetail

cex.lab <- 1.5
cex.yaxis <- 1.2
cex.xaxis <- 0.5
cex.title <- 2

mixtures <- LateAugustMixtures2014Strata
for(mixture in mixtures) {
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  Barplot1 <- barplot2(height = t(matrix(c(LateAugust2014Strata_89loci_EstimatesStats[[mixture]][, "median"],
                                           LateAugust2014Strata_46loci_EstimatesStats[[mixture]][, "median"]), ncol = 2)* 100), 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(matrix(c(LateAugust2014Strata_89loci_EstimatesStats[[mixture]][, "5%"],
                                         LateAugust2014Strata_46loci_EstimatesStats[[mixture]][, "5%"]), ncol = 2)* 100), 
                       ci.u = t(matrix(c(LateAugust2014Strata_89loci_EstimatesStats[[mixture]][, "95%"],
                                         LateAugust2014Strata_46loci_EstimatesStats[[mixture]][, "95%"]), ncol = 2)* 100), 
                       ylim = c(0, 100), col = colorpanel(2, low = "blue", high = "white"), cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(x = 1, y = 95, legend = c("loci89", "loci46"), fill = colorpanel(2, low = "blue", high = "white"), border = "black", bty = "n", cex = 1, title = LateAugustMixtures2014StrataDetail[1])
  abline(h = 0)
  mtext(text = "Percentage of Catch", side = 2, line = 2.5, cex = cex.lab)
  mtext(text = "Reporting Group", side = 1, line = 3, cex = cex.lab)
  mtext(text = KMA15GroupsPC2Rows, side = 1, line = 1, at = colMeans(Barplot1), adj = 0.5, cex = cex.xaxis)
  mtext(text = LateAugustMixtures2014StrataDetail[mixture], side = 3, cex = cex.title)
}








plot.ci.extract <- function(x, mixture, group, stat.col) {
  list.index <- grep(pattern = mixture, x = names(x))  
  sapply(list.index, function(i) {x[[i]][group, stat.col]} )
}



emf(file = "V:/Presentations/Regional/4_Westward/Sockeye/COMFISH/KarlukScenarioFigure.emf", width = 7, height = 7, family = "sans", bg = "white")

# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             1, 3,
                             1, 4,
                             5, 6), nrow = 4, ncol = 2, byrow = TRUE)

PCGroups15
PCGroups15Rows2 <- c("West of\nChignik", "Black\nLake", "Chignik\nLake", "U. Station\nAkalura", "Frazer\n", "Ayakulik\n", "Karluk\n", "Uganik\n", "Northwest\nKodiak", "Afognak\n", "Eastside\nKodiak", "Saltery\n", "Cook\nInlet", "PWS\n", "South of\nCape Suckling")
cex.lab <- 1.5
cex.yaxis <- 1.2
cex.xaxis <- 0.5
ci.lwd <- 2.5
TrueKarlukProportions <- MixtureProofTestProportions[FisheryProofTestScenarioNames[3], ] * 100



layout(mat = layoutmat, widths = c(0.1, 1, 1), heights = c(0.9, 0.9, 1, 0.1))
par(mar = rep(0, 4))

# Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)

# Hypothetical Late Karluk Replicate 1
par(mar = c(1, 1, 1, 1))
Barplot1 <- barplot2(height = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats$late.Karluk1[, "median"] * 100, 
                     beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                     ci.l = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats$late.Karluk1[, "5%"] * 100, 
                     ci.u = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats$late.Karluk1[, "95%"] * 100, 
                     ylim = c(0, 100), col = colorpanel(1, low = "blue", high = "white"), cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
segments(x0 = Barplot1-0.5, y0 = TrueKarlukProportions, x1 = Barplot1+0.5, y1 = TrueKarlukProportions, col = "red", lwd = ci.lwd)
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = "Hypothetical August", x = "topleft", fill = colorpanel(1, low = "blue", high = "white"), border = "black", bty = "n", cex = 1, title="Replicate 1")
abline(h = 0)

# Hypothetical Late Karluk Replicate 2
par(mar = c(1, 1, 1, 1))
Barplot2 <- barplot2(height = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats$late.Karluk2[, "median"] * 100, 
                     beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                     ci.l = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats$late.Karluk2[, "5%"] * 100, 
                     ci.u = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats$late.Karluk2[, "95%"] * 100, 
                     ylim = c(0, 100), col = colorpanel(1, low = "blue", high = "white"), cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
segments(x0 = Barplot1-0.5, y0 = TrueKarlukProportions, x1 = Barplot1+0.5, y1 = TrueKarlukProportions, col = "red", lwd = ci.lwd)
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = "Hypothetical August", x = "topleft", fill = colorpanel(1, low = "blue", high = "white"), border = "black", bty = "n", cex = 1, title="Replicate 2")
abline(h = 0)

# Hypothetical Late Karluk Replicate 3
par(mar = c(2, 1, 1, 1))
Barplot3 <- barplot2(height = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats$late.Karluk3[, "median"] * 100, 
                     beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                     ci.l = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats$late.Karluk3[, "5%"] * 100, 
                     ci.u = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats$late.Karluk3[, "95%"] * 100, 
                     ylim = c(0, 100), col = colorpanel(1, low = "blue", high = "white"), cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
segments(x0 = Barplot1-0.5, y0 = TrueKarlukProportions, x1 = Barplot1+0.5, y1 = TrueKarlukProportions, col = "red", lwd = ci.lwd)
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = "Hypothetical August", x = "topleft", fill = colorpanel(1, low = "blue", high = "white"), border = "black", bty = "n", cex = 1, title="Replicate 3")
abline(h = 0)
mtext(text = PCGroups15Rows2, side = 1, line = 1, at = Barplot3[, 1], adj = 0.5, cex = cex.xaxis)

## Blank Corner
par(mar = rep(0, 4))
plot.new()

## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.5, labels = "Reporting Group", cex = cex.lab)


dev.off()







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare markersets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Define layout structure
FiveByThreeMatrix <- matrix(data=c(1,2,3,4,
                                   1,5,6,7,
                                   1,8,9,10,
                                   1,11,12,13,
                                   1,14,15,16,
                                   17,18,18,18),nrow=6,ncol=4,byrow=TRUE)
FiveByThreeLayout <- layout(mat=FiveByThreeMatrix,widths=c(0.2,rep(1,3)),heights=c(1,1,1,1,1,0.2))
layout.show(n=18)

#~~~~~~~~~~~~~~~~~~
# Create 5%, median, 95% extractor function
plot.ci.extract <- function(x, mixture, group, stat.col) {
  list.index <- grep(pattern = mixture, x = names(x))  
  sapply(list.index, function(i) {x[[i]][group, stat.col]} )
}

#~~~~~~~~~~~~~~~~~~
# Create a figure wrapup funtion
require(gplots)
require(devEMF)

percent = TRUE
mixtures = LateAugustMixtures2014Strata

col = c("black", "grey70")
cex.pch = 1.5
ci.lwd = 2
sfrac = 0.02


for(mixture in mixtures) {
  
  x.loci89 <- dget(file = "Estimates objects/LateAugust2014Strata_89loci_EstimatesStats.txt")
  x.loci46 <- dget(file = "Estimates objects/LateAugust2014Strata_46loci_EstimatesStats.txt")
  
  if(percent){
    percent.multiplier <- 100
    x.loci89 <- sapply(x.loci89, function(rpt) {rpt[, 1:5] * percent.multiplier}, simplify = FALSE)
    x.loci46 <- sapply(x.loci46, function(rpt) {rpt[, 1:5] * percent.multiplier}, simplify = FALSE)
    file = paste("Figures/", mixture, ".Percent.emf", sep = '')
    y.label = "Percent"
  } else {
    percent.multiplier <- 1
    file = paste("Figures/", mixture, ".emf", sep = '')
    y.label = "Proportion"
  }
  
  MixtureProofTestProportions <- t(dget(file = "Objects/MixtureProofTestProportions.txt")) * percent.multiplier
  
  emf(file = file, width = 6.5, height = 8, family = "Times")
  #png(file = paste(filedir, "/", mixture, ".png", sep = ""), width = 6.5, height = 7.5, units = "in", res = 1200, pointsize = 4, family = "Times")
  layout(mat = FiveByThreeMatrix, widths = c(0.3, rep(1, 3)), heights = c(1, 1, 1, 1, 1, 0.4))
  
  # Plot 1 - y-axis title
  par(mar = c(0, 0, 0, 0))
  par(srt = 90)
  plot.new()
  text(labels = y.label, x = 0.2, y = 0.5, adj = c(0.5, 0.5), font = 1, cex = 2)
  
  # Plot 2-16 - Groups
  sapply(seq(groups), function(j) {
    
    y <- c(t(matrix(data = c(plot.ci.extract(x = x.loci89, mixture = mixture, group = groups[j], stat.col = "median"),
                             plot.ci.extract(x = x.loci46, mixture = mixture, group = groups[j], stat.col = "median")),
                    nrow = 5)))
    ui <- c(t(matrix(data = c(plot.ci.extract(x = x.loci89, mixture = mixture, group = groups[j], stat.col = "95%"),
                              plot.ci.extract(x = x.loci46, mixture = mixture, group = groups[j], stat.col = "95%")),
                     nrow = 5)))
    li <- c(t(matrix(data = c(plot.ci.extract(x = x.loci89, mixture = mixture, group = groups[j], stat.col = "5%"),
                              plot.ci.extract(x = x.loci46, mixture = mixture, group = groups[j], stat.col = "5%")),
                     nrow = 5)))
    
    par(mar = c(0, 0, 0, 0))
    par(srt = 0)
    if(j == 1) {
      plotCI(x = c(1:2, 4:5, 7:8, 10:11, 13:14), y = y, ui = ui, li = li, xlim = c(0.5, 14.5), ylim = c(0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", cex = cex.pch, col = col, lwd = ci.lwd, sfrac = sfrac)
      legend("topleft", legend = c("89 loci", "46 loci"), fill = col, bty = 'n', cex = 1.2)
    }
    
    if(j %in% c(4, 7, 10)) {
      plotCI(x = c(1:2, 4:5, 7:8, 10:11, 13:14), y = y, ui = ui, li = li, xlim = c(0.5, 14.5), ylim = c(0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", cex = cex.pch, col = col, lwd = ci.lwd, sfrac = sfrac)
    }
    
    if(j == 13) {
      plotCI(x = c(1:2, 4:5, 7:8, 10:11, 13:14), y = y, ui = ui, li = li, xlim = c(0.5, 14.5), ylim = c(0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", cex = cex.pch, col = col, lwd = ci.lwd, sfrac = sfrac)
      axis(side = 1, at = seq(from = 1.5, by = 3, length.out = 5), labels = 1:5)
    }
    
    if(j %in% c(14, 15)) {
      plotCI(x = c(1:2, 4:5, 7:8, 10:11, 13:14), y = y, ui = ui, li = li, xlim = c(0.5, 14.5), ylim = c(0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", yaxt = "n", cex = cex.pch, col = col, lwd = ci.lwd, sfrac = sfrac)
      axis(side = 1, at = seq(from = 1.5, by = 3, length.out = 5), labels = 1:5)
    }
    
    if(j %in% c(2, 3, 5, 6, 8, 9, 11, 12)) {
      plotCI(x = c(1:2, 4:5, 7:8, 10:11, 13:14), y = y, ui = ui, li = li, xlim = c(0.5, 14.5), ylim = c(0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", yaxt = "n", cex = cex.pch, col = col, lwd = ci.lwd, sfrac = sfrac)
    }
    
    abline(h = MixtureProofTestProportions[j, mixture], col = "red", lwd = 2)
    abline(h = 0, lwd = 1, lty = 2)
    text(PCGroups15[j], x = 14, y = 0.6 * percent.multiplier, adj = 1, cex = 1.2)
  } )
  
  # Plot 17 - Blank Corner
  plot.new()
  
  ## Plot 18 - x-axis title
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(labels = "Replicate", x = 0.5, y = 0.25, adj = c(0.5, 0.5), font = 1, cex = 2)
  
  dev.off()
}






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Doing first table old school to create the sheet

# THD and KRS agree that these genetic samples represent Harvest at the end of the Late stratum (8/25-8/29)
Karluk2014LateAugustHarvests <- sum(as.numeric(readClipboard()))
dput(x = Karluk2014LateAugustHarvests, file = "Objects/Karluk2014LateAugustHarvests.txt")

Uganik2014LateAugustHarvests <- sum(as.numeric(readClipboard()))
dput(x = Uganik2014LateAugustHarvests, file = "Objects/Uganik2014LateAugustHarvests.txt")

Uyak2014LateAugustHarvests <- sum(as.numeric(readClipboard()))
dput(x = Uyak2014LateAugustHarvests, file = "Objects/Uyak2014LateAugustHarvests.txt")


LateAugustMixtures2014Strata
KMAsockeyeNames <- EASSIP14Groups[GroupOrder]
LateAugustMixtures2014Strata_SampleSizes


# Create Harvest Object
LateAugustHarvests <- c(Karluk2014LateAugustHarvests, Uganik2014LateAugustHarvests, Uyak2014LateAugustHarvests)
names(LateAugustHarvests) <- LateAugustMixtures2014Strata
LateAugust_Estimates_Harvest <- setNames(object = lapply(LateAugustMixtures2014Strata, function (Mix) {LateAugust_Estimates$Stats[[Mix]] * LateAugustHarvests[[Mix]]}), nm = LateAugustMixtures2014Strata)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Defining caption variables
require(xlsx)

SheetNames <- c("Karluk", "Uganik", "Uyak")
names(SheetNames) <- LateAugustMixtures2014Strata

# Karluk
Stratum <- "stratum 4 (August 25-29; Harvest=55,637; n=284)"
Area <- "Karluk Section, "
Year <- "2014"

# Uyak
Stratum <- "stratum 4 (August 25-29; Harvest=70,084; n=280)"
Area <- "Uyak Section, "
Year <- "2014"

# Uganik
Stratum <- "stratum 4 (August 25-29; Harvest=126,310; n=282)"
Area <- "Uganik Section, "
Year <- "2014"




Mix <- LateAugustMixtures2014Strata[2]

KarlukCaption1 <- paste("Table X.-Estimates of stock composition (%) and stock-specific harvest for temporal ",Stratum," of the ",Area,Year,".  Estimates include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.",sep='')
Disclaimer <- cbind("Note: Stock composition estimates may not sum to 100% and stock-specific harvest estimates may not sum to the total harvest due to rounding error.",'','','','','','','','','','','','')

KarlukTable1=cbind(KarlukCaption1,'','','','','','','','','','','','')
KarlukTable1=rbind(KarlukTable1,cbind("","Stock Composition",'','','','','','',"Stock-specific Harvest",'','','',''))
KarlukTable1=rbind(KarlukTable1,cbind('','',"90% CI",'','','','','','',"90% CI",'','',''))
KarlukTable1=rbind(KarlukTable1,cbind("Reporting Group","Median","5%","95%","P=0","Mean","SD",'',"Median","5%","95%","Mean","SD"))
## West of Chignik
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[1],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[1],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[1],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[1],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[1],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[1],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[1],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[1],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[1],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[1],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[1],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[1],2],format="f",digits=0,big.mark=",")))
## Black Lake
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[2],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[2],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[2],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[2],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[2],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[2],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[2],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[2],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[2],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[2],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[2],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[2],2],format="f",digits=0,big.mark=",")))
## Chignik Lake
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[3],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[3],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[3],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[3],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[3],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[3],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[3],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[3],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[3],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[3],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[3],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[3],2],format="f",digits=0,big.mark=",")))
## Upper Station / Akalura
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[4],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[4],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[4],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[4],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[4],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[4],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[4],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[4],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[4],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[4],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[4],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[4],2],format="f",digits=0,big.mark=",")))
## Frazer
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[5],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[5],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[5],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[5],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[5],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[5],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[5],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[5],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[5],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[5],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[5],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[5],2],format="f",digits=0,big.mark=",")))
## Ayakulik
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[6],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[6],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[6],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[6],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[6],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[6],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[6],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[6],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[6],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[6],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[6],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[6],2],format="f",digits=0,big.mark=",")))
## Karluk
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[7],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[7],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[7],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[7],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[7],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[7],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[7],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[7],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[7],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[7],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[7],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[7],2],format="f",digits=0,big.mark=",")))
## Uganik
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[8],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[8],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[8],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[8],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[8],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[8],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[8],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[8],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[8],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[8],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[8],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[8],2],format="f",digits=0,big.mark=",")))
## Northwest Minor
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[9],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[9],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[9],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[9],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[9],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[9],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[9],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[9],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[9],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[9],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[9],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[9],2],format="f",digits=0,big.mark=",")))
## Afognak
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[10],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[10],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[10],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[10],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[10],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[10],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[10],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[10],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[10],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[10],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[10],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[10],2],format="f",digits=0,big.mark=",")))
## Eastside Minor
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[11],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[11],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[11],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[11],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[11],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[11],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[11],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[11],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[11],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[11],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[11],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[11],2],format="f",digits=0,big.mark=",")))
## Saltery
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[12],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[12],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[12],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[12],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[12],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[12],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[12],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[12],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[12],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[12],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[12],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[12],2],format="f",digits=0,big.mark=",")))
## Cook Inlet
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[13],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[13],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[13],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[13],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[13],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[13],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[13],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[13],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[13],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[13],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[13],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[13],2],format="f",digits=0,big.mark=",")))
## PWS / Copper River
KarlukTable1=rbind(KarlukTable1,cbind(KMAsockeyeNames[14],
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[14],3]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[14],4]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[14],5]*100,format="f",digits=1),formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[14],6],format="f",digits=2),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[14],1]*100,format="f",digits=1),
                                      formatC(x=LateAugust_Estimates$Stats[[Mix]][GroupOrder[14],2]*100,format="f",digits=1),'',formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[14],3],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[14],4],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[14],5],format="f",digits=0,big.mark=","),
                                      formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[14],1],format="f",digits=0,big.mark=","),formatC(x=LateAugust_Estimates_Harvest[[Mix]][GroupOrder[14],2],format="f",digits=0,big.mark=",")))
## Harvest[[Mix]] sum row
KarlukTable1=rbind(KarlukTable1,cbind('','','','','','','','','','','Total',formatC(x=LateAugustHarvests[[Mix]],format="f",digits=0,big.mark=","),''),Disclaimer)
write.xlsx(x=as.data.frame(KarlukTable1),file="Estimates tables/LateAugust2014Tables.xlsx",col.names=F,row.names=F,append=T,sheetName=SheetNames[[Mix]])




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Look at harvest by Section for Karluk, Uyak, and Uganik ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KarlukHarvest <- as.numeric(readClipboard())
KarlukHarvest <- data.frame("Date" = seq(from = as.Date("06/01/2014", "%m/%d/%Y"), by = "days", length.out = length(KarlukHarvest)), "Harvest" = KarlukHarvest)
dput(x = KarlukHarvest, file = "Objects/KarlukHarvest.txt")

UyakHarvest <- as.numeric(readClipboard())
UyakHarvest <- data.frame("Date" = seq(from = as.Date("06/01/2014", "%m/%d/%Y"), by = "days", length.out = length(UyakHarvest)), "Harvest" = UyakHarvest)
dput(x = UyakHarvest, file = "Objects/UyakHarvest.txt")

UganikHarvest <- as.numeric(readClipboard())
UganikHarvest <- data.frame("Date" = seq(from = as.Date("06/01/2014", "%m/%d/%Y"), by = "days", length.out = length(UganikHarvest)), "Harvest" = UganikHarvest)
dput(x = UganikHarvest, file = "Objects/UganikHarvest.txt")

max(c(KarlukHarvest$Harvest, UyakHarvest$Harvest, UganikHarvest$Harvest), na.rm = TRUE)

plot(Harvest ~ Date, data = KarlukHarvest, type = "l", lwd = 3, col ="purple", bty = "n", ylim = c(0, 50000), cex.axis = 1.5, cex.lab = 1.7, main = "2014 Harvest", cex.main = 2.5)
points(Harvest ~ Date, data = UyakHarvest, type = "l", lwd = 3, col ="grey50")
points(Harvest ~ Date, data = UganikHarvest, type = "l", lwd = 3, col ="orange")
legend("topleft", legend = c("Karluk", "Uyak", "Uganik"), fill = c("purple", "grey50", "orange"), bty = "n", cex = 2)
abline(v = as.Date("08/26/2014", "%m/%d/%Y"), lwd = 5)
text(x = as.Date("08/26/2014", "%m/%d/%Y"), y = 4e4, labels = "Last samples taken\n8/26-27/2014", pos = 2, cex = 1.5)

polygon(x = c(as.Date("08/26/2014", "%m/%d/%Y"), as.Date("08/26/2014", "%m/%d/%Y"), as.Date("10/02/2014", "%m/%d/%Y"), as.Date("10/02/2014", "%m/%d/%Y")), y = c(-2e3, 5e4, 5e4, -2e3), col=rgb(0.5, 0.5, 0.5, alpha = 0.7))
text(x = as.Date("08/26/2014", "%m/%d/%Y"), y = 5.1e4, labels = "Unsampled", pos = 4, cex = 1.5, offset = 2)