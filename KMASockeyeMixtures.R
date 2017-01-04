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
#### Repeat Late August 2014 Analyses ####
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
#### Clean workspace; dget .gcl objects and Locus Control
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl96.txt")
LateAugustMixtures2014 <- dget(file = "Objects/EASSIP-LateAugust2014/LateAugustMixtures2014.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("EASSIP-LateAugust2014", "OriginalLocusControl.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered mixtures
invisible(sapply(LateAugustMixtures2014, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections/KMALateAugust2014/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata
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
#### Data QC/Massage
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
#### Combine Loci
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
#### Clean workspace; dget .gcl objects and Locus Control
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


## Get post-QC, stratified, combined loci, mixtures
invisible(sapply(LateAugustMixtures2014Strata, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_Strata_PostQC_CombinedLoci/", silly, ".txt", sep = "")), pos = 1)} )); beep(4)
objects(pattern = "\\.gcl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get baseline objects needed for MSA
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects")

KMA473Pops15FlatPrior <- dget(file = "KMA473Pops15FlatPrior.txt")
KMA473PopsInits <- dget(file = "KMA473PopsInits.txt")
KMA473PopsGroupVec14 <- dget(file = "KMA473PopsGroupVec14.txt")
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
file.copy(from = c("KMA473Pops15FlatPrior.txt", "KMA473PopsInits.txt", "KMA473PopsGroupVec14.txt", "KMA473Pops.txt", "PCGroups15.txt", "CommonNames473.txt", "Colors15.txt",
                   "V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt"), 
          to = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects")
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")


KMA15GroupsPC2Rows <- PCGroups15Rows2 <- c("West of\nChignik", "Black\nLake", "Chignik\nLake", "U. Station\nAkalura", "Frazer\n", "Ayakulik\n", "Karluk\n", "Uganik\n", "Northwest\nKodiak", "Afognak\n", "Eastside\nKodiak", "Saltery\n", "Cook\nInlet", "PWS\n", "South of\nCape Suckling")
dput(x = KMA15GroupsPC2Rows, file = "Objects/KMA15GroupsPC2Rows.txt")
# Note, I changed the names for CommonNames473.txt -> KMA473CommonNames.txt & PCGroups15 -> KMA15GroupsPC.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MSA files for BAYES
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
                        groupvec = KMA473PopsGroupVec14, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits, dir = "BAYES/Late August 89loci/Control",
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
                        groupvec = KMA473PopsGroupVec14, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits, dir = "BAYES/Late August 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(LateAugustMixtures2014Strata, function(Mix) {dir.create(paste("BAYES/Late August 46loci/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Output
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
#### Plot results
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
dput(x = LateAugustMixtures2014StrataDetail, file = "Objects/LateAugustMixtures2014StrataDetail.txt")

ViolinPlot <- function(estimates, groups, colors, header, wex = 1, thin = 10) {
  while(!require(vioplot, quietly = TRUE)){install.packages("vioplot")}
  
  mixtures <- names(estimates$Stats)
  
  par(mar = c(5.1, 4.6, 3.6, 1.1))
  sapply(mixtures, function(Mix) {
    plot(estimates$Stats[[Mix]][, "median"], cex = 3, pch = 16, col = colors, ylab = "Proportion of Mixture", ylim = c(0, 1), xlab = "", axes = FALSE, main = header[[Mix]], cex.main = 2, cex.lab = 1.5)
    sapply(seq(groups), function(i) {vioplot2(estimates$Output[[Mix]][seq(from = 1, to = nrow(estimates$Output[[Mix]]), by = thin), i], at = i, horizontal = FALSE, col = colors[i], border = TRUE, drawRect = FALSE, rectCol = colors[i], add = TRUE, wex = wex, lwd = 2)})
    arrows(x0 = seq(groups), y0 = estimates$Stats[[Mix]][, "5%"], x1 = seq(groups), y1 = estimates$Stats[[Mix]][, "95%"], angle = 90, code = 3, length = 0.2, lwd = 2)
    points(estimates$Stats[[Mix]][, "median"], cex = 2, pch = 21, col = "white", bg = colors, lwd = 3)
    axis(side = 2, lwd = 3, cex.axis = 1.5)
    text(x = (seq(groups)) - 0.35, y = 0, labels = groups, srt = 60, pos = 1, offset = 2.5, xpd = TRUE)
    axis(side = 1, labels = NA, at = seq(groups), pos = 0, lwd = 2, tick = FALSE)
    abline(h = 0, lwd = 3)
  } )
  
}

ViolinPlot(estimates = LateAugust2014Strata_89loci_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = LateAugustMixtures2014StrataDetail, thin = 10)
ViolinPlot(estimates = LateAugust2014Strata_46loci_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = LateAugustMixtures2014StrataDetail)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare markersets, side by side barplots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries
require(devEMF)
require(gplots)
require(plotrix)

# Bring in EASSIP estimates to compare as well
LateAugust2014_EASSIP_Estimates <- dget(file = "Estimates objects/LateAugust2014_EASSIP_Estimates.txt")
LateAugust2014_EASSIP_EstimatesStats <- LateAugust2014_EASSIP_Estimates$Stats

names(LateAugust2014_EASSIP_EstimatesStats) <- LateAugustMixtures2014Strata

for(mix in names(LateAugust2014_EASSIP_EstimatesStats)) {
  LateAugust2014_EASSIP_EstimatesStats[[mix]] <- LateAugust2014_EASSIP_EstimatesStats[[mix]][c(1, 2, 3, 12, 11, 4, 5, 7, 6, 8, 9, 10, 13, 14), ]
  LateAugust2014_EASSIP_EstimatesStats[[mix]] <- rbind(LateAugust2014_EASSIP_EstimatesStats[[mix]], rep(0, ncol(LateAugust2014_EASSIP_EstimatesStats[[mix]])))
  rownames(LateAugust2014_EASSIP_EstimatesStats[[mix]]) <- KMA15GroupsPC
}
str(LateAugust2014_EASSIP_EstimatesStats)

# Get estimates stats for loci89 and loci46 on KMA baseline
LateAugust2014Strata_89loci_EstimatesStats <- dget(file = "Estimates objects/LateAugust2014Strata_89loci_EstimatesStats.txt")
LateAugust2014Strata_46loci_EstimatesStats <- dget(file = "Estimates objects/LateAugust2014Strata_46loci_EstimatesStats.txt")


KMA15GroupsPC
KMA15GroupsPC2Rows
LateAugustMixtures2014StrataDetail

cex.lab <- 1.2
cex.yaxis <- 1
cex.xaxis <- 0.45
cex.title <- 1.5
ci.lwd <- 2
n.base <- 3

#~~~~~~~~~~~~~~~~~~
# Make .emf of barplots
for(mixture in LateAugustMixtures2014Strata) {
  emf(file = paste("Figures/", mixture, "_barplot.emf", sep = ''), width = 5.75, height = 5.5)
  par(mar = c(2.6, 3.6, 3.1, 0.6))
  Barplot1 <- barplot2(height = t(matrix(c(LateAugust2014_EASSIP_EstimatesStats[[mixture]][, "median"],
                                           LateAugust2014Strata_89loci_EstimatesStats[[mixture]][, "median"],
                                           LateAugust2014Strata_46loci_EstimatesStats[[mixture]][, "median"]), ncol = n.base)* 100), 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(matrix(c(LateAugust2014_EASSIP_EstimatesStats[[mixture]][, "5%"],
                                         LateAugust2014Strata_89loci_EstimatesStats[[mixture]][, "5%"],
                                         LateAugust2014Strata_46loci_EstimatesStats[[mixture]][, "5%"]), ncol = n.base)* 100), 
                       ci.u = t(matrix(c(LateAugust2014_EASSIP_EstimatesStats[[mixture]][, "95%"],
                                         LateAugust2014Strata_89loci_EstimatesStats[[mixture]][, "95%"],
                                         LateAugust2014Strata_46loci_EstimatesStats[[mixture]][, "95%"]), ncol = n.base)* 100), 
                       ylim = c(0, 100), col = colorpanel(n.base, low = "blue", high = "white"), cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(x = 0, y = 100, legend = c("EASSIP loci91", "KMA loci89", "KMA loci46"), fill = colorpanel(n.base, low = "blue", high = "white"), border = "black", bty = "n", cex = 1, title = NULL)
  abline(h = 0, xpd = FALSE)
  mtext(text = "Percentage of Catch", side = 2, line = 2.5, cex = cex.lab)
  mtext(text = "Reporting Group", side = 1, line = 1.5, cex = cex.lab)
#   text(x = colMeans(Barplot1), y = -1, labels = KMA15GroupsPC2Rows, cex = cex.xaxis, srt = 90, xpd = TRUE, adj = 1)
  mtext(text = KMA15GroupsPC2Rows, side = 1, line = 0.5, at = colMeans(Barplot1), adj = 0.5, cex = cex.xaxis)
  mtext(text = LateAugustMixtures2014StrataDetail[mixture], side = 3, cex = cex.title)
  dev.off()
}

#~~~~~~~~~~~~~~~~~~
# Make .emf of barplots
for(mixture in LateAugustMixtures2014Strata) {
  png(file = paste("Figures/", mixture, "_barplot.png", sep = ''), width = 5.75, height = 5.5, units = "in", res = 300)
  par(mar = c(2.6, 3.6, 3.1, 0.6))
  Barplot1 <- barplot2(height = t(matrix(c(LateAugust2014_EASSIP_EstimatesStats[[mixture]][, "median"],
                                           LateAugust2014Strata_89loci_EstimatesStats[[mixture]][, "median"],
                                           LateAugust2014Strata_46loci_EstimatesStats[[mixture]][, "median"]), ncol = n.base)* 100), 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(matrix(c(LateAugust2014_EASSIP_EstimatesStats[[mixture]][, "5%"],
                                         LateAugust2014Strata_89loci_EstimatesStats[[mixture]][, "5%"],
                                         LateAugust2014Strata_46loci_EstimatesStats[[mixture]][, "5%"]), ncol = n.base)* 100), 
                       ci.u = t(matrix(c(LateAugust2014_EASSIP_EstimatesStats[[mixture]][, "95%"],
                                         LateAugust2014Strata_89loci_EstimatesStats[[mixture]][, "95%"],
                                         LateAugust2014Strata_46loci_EstimatesStats[[mixture]][, "95%"]), ncol = n.base)* 100), 
                       ylim = c(0, 100), col = colorpanel(n.base, low = "blue", high = "white"), cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(x = 0, y = 100, legend = c("EASSIP loci91", "KMA loci89", "KMA loci46"), fill = colorpanel(n.base, low = "blue", high = "white"), border = "black", bty = "n", cex = 1, title = NULL)
  abline(h = 0, xpd = FALSE)
  mtext(text = "Percentage of Catch", side = 2, line = 2.5, cex = cex.lab)
  mtext(text = "Reporting Group", side = 1, line = 1.5, cex = cex.lab)
  #   text(x = colMeans(Barplot1), y = -1, labels = KMA15GroupsPC2Rows, cex = cex.xaxis, srt = 90, xpd = TRUE, adj = 1)
  mtext(text = KMA15GroupsPC2Rows, side = 1, line = 0.5, at = colMeans(Barplot1), adj = 0.5, cex = cex.xaxis)
  mtext(text = LateAugustMixtures2014StrataDetail[mixture], side = 3, cex = cex.title)
  dev.off()
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pull genotypes from LOKI 2014/2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get collection SILLYs
KMA2014 <- c("SALITC14", "SAYAKC14", "SKARLC14", "SUGANC14", "SUYAKC14")  # No Igvak samples, no fishing
dput(x = KMA2014, file = "Objects/KMA2014.txt")
KMA2015 <- c("SALITC15", "SAYAKC15", "SKARLC15", "SUGANC15", "SUYAKC15")  # Good Igvak samples, but very low fishing fishing
dput(x = KMA2015, file = "Objects/KMA2015.txt")
# KMA2016 <- c("SALITC16", "SAYAKC16", "SIGVAC16", "SKARLC16", "SUGANC16", "SUYAKC16")  # We'll see
# dput(x = KMA2016, file = "Objects/KMA2016.txt")

## Create Locus Control
CreateLocusControl.GCL(markersuite = "Sockeye_Kodiak_48SNPs", username = username, password = password)

## Save original LocusControl
loci48 <- LocusControl$locusnames
mito.loci48 <- which(LocusControl$ploidy == 1)

dput(x = LocusControl, file = "Objects/OriginalLocusControl48.txt")
dput(x = loci48, file = "Objects/loci48.txt")
dput(x = mito.loci48, file = "Objects/mito.loci48.txt")

#~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
# LOKI2R.GCL(sillyvec = c(KMA2014, KMA2015), username = username, password = password)
sapply(c(KMA2014, KMA2015), function(silly) {LOKI2R.GCL(sillyvec = silly, username = username, password = password)} )  # looping through due to heap space error when getting all sillys at once
rm(username, password)
objects(pattern = "\\.gcl")

## Save unaltered .gcl's as back-up:
invisible(sapply(c(KMA2014, KMA2015), function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/" , silly, ".txt", sep = ''))} )); beep(8)

## Original sample sizes by SILLY
collection.size.original <- sapply(c(KMA2014, KMA2015), function(silly) get(paste(silly, ".gcl", sep = ""))$n)
matrix(data = collection.size.original, ncol = 2, dimnames = list(c("Alitak", "Ayakulik", "Karluk", "Uganik", "Uyak"), 2014:2015))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl48.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("EASSIP-LateAugust2014", "OLD 15RG", "OriginalLocusControl96.txt", "OriginalLocusControl48.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered mixtures
invisible(sapply(c(KMA2014, KMA2015), function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## All fish have a capture date?
sapply(c(KMA2014, KMA2015), function(silly) {sum(is.na(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE))} )  # Zeros are good

## Confirming samples sizes by date
sapply(c(KMA2014, KMA2015), function(silly) {table(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE)} )


## Get dataframes of strata dates
KMA.Strata.Dates.2014Alitak <- read.xlsx(file = "MixtureStrataDates.xlsx", sheetName = "2014Alitak", header = TRUE)
KMA.Strata.Dates.2014Ayak <- read.xlsx(file = "MixtureStrataDates.xlsx", sheetName = "2014Ayak", header = TRUE)
KMA.Strata.Dates.2014KarlUganUyak <- read.xlsx(file = "MixtureStrataDates.xlsx", sheetName = "2014KarlUganUyak", header = TRUE)
KMA.Strata.Dates.2015 <- read.xlsx(file = "MixtureStrataDates.xlsx", sheetName = "2015", header = TRUE)


## Function to define strata by dates (date.df)

# # Inputs
# silly <- "SKARLC14"
# date.df <- KMA.Strata.Dates.2014KarlUganUyak
# loci <- loci48

PoolCollectionsByDateDF <- function(silly, date.df, loci) {
  sapply(silly, function(mix) {
    mix.dates <- unique(as.Date(get(paste(mix, ".gcl", sep = ''))$attributes$CAPTURE_DATE))
    by(data = date.df, INDICES = date.df$Strata, function(x) {
      IDs <- AttributesToIDs.GCL(silly = mix, attribute = "CAPTURE_DATE", matching = mix.dates[mix.dates >= x$Begin & mix.dates <= x$End])
      IDs <- list(na.omit(IDs))
      names(IDs) <- mix
      PoolCollections.GCL(collections = mix, loci = loci, IDs = IDs, newname = paste(mix, as.character(x$Strata), sep = "_"))
      list("First Last Fish" = range(as.numeric(unlist(IDs))), "n" = get(paste(mix, "_", as.character(x$Strata), ".gcl", sep = ''))$n)
    } )
  }, simplify = FALSE, USE.NAMES = TRUE)
}

# # Example
# PoolCollectionsByDateDF(silly = LateAugustMixtures2014[c(1, 3)], date.df = KMA.Strata.Dates.2014KarlUganUyak, loci = loci48)

PoolCollectionsByDateDF(silly = KMA2014[1], date.df = KMA.Strata.Dates.2014Alitak, loci = loci48)
PoolCollectionsByDateDF(silly = KMA2014[2], date.df = KMA.Strata.Dates.2014Ayak, loci = loci48)
PoolCollectionsByDateDF(silly = KMA2014[3:5], date.df = KMA.Strata.Dates.2014KarlUganUyak, loci = loci48)
PoolCollectionsByDateDF(silly = KMA2015, date.df = KMA.Strata.Dates.2015, loci = loci48)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mixture sillyvec
KMA2014Strata <- unlist(strsplit(x = grep(pattern = "14_", x = objects(pattern = "\\.gcl"), value = TRUE), split = "\\.gcl"))
KMA2015Strata <- unlist(strsplit(x = grep(pattern = "15_", x = objects(pattern = "\\.gcl"), value = TRUE), split = "\\.gcl"))
KMA2014_2015Strata <- c(KMA2014Strata, KMA2015Strata)

dput(x = KMA2014Strata, file = "Objects/KMA2014Strata.txt")
dput(x = KMA2015Strata, file = "Objects/KMA2015Strata.txt")
dput(x = KMA2014_2015Strata, file = "Objects/KMA2014_2015Strata.txt")

# Confirm sample sizes
sapply(KMA2014_2015Strata, function(silly) get(paste(silly, ".gcl", sep = ""))$n)

# View as tables by year
require(reshape)
samp.df.2014 <- data.frame(t(sapply(KMA2014Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2014, location ~ temporal.strata)

samp.df.2015 <- data.frame(t(sapply(KMA2015Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2015, location ~ temporal.strata)

# dput mixture sillys
invisible(sapply(KMA2014_2015Strata, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_Strata/" , silly, ".txt", sep = ''))} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

KMA2014_2015Strata

KMA2014_2015Strata_SampleSizes <- matrix(data = NA, nrow = length(KMA2014_2015Strata), ncol = 5, 
                                         dimnames = list(KMA2014_2015Strata, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_KMA2014_2015Strata_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = KMA2014_2015Strata, loci = loci48)
min(Original_KMA2014_2015Strata_SampleSizebyLocus)  ## 267/285
apply(Original_KMA2014_2015Strata_SampleSizebyLocus, 1, min) / apply(Original_KMA2014_2015Strata_SampleSizebyLocus, 1, max)  ## Good, 0.947

Original_KMA2014_2015Strata_PercentbyLocus <- apply(Original_KMA2014_2015Strata_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_KMA2014_2015Strata_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_KMA2014_2015Strata_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


#### Check individuals
## View Histogram of Failure Rate by Strata
invisible(sapply(KMA2014_2015Strata, function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  failure <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(ind == "0") / length(ind)} )
  hist(x = failure, main = mix, xlab = "Failure Rate", col = 8, xlim = c(0, 1), ylim = c(0, 20), breaks = seq(from = 0, to = 1, by = 0.02))
  abline(v = 0.2, lwd = 3)
}))

## Experimental method to look for a way to remove "missing" loci fish based on a strata specific failure rate (rather than 0.8 globally, as done normally)
# Experimented with loess, smooth.split, polynomial regressions, outlier (Franz), and Cooks Distance methods to model outliers.
# Thought of a clever idea to look at the "shoulder" of the cumulative success rate, but this concept is a bit weird.
# Planning to keep it simple!
invisible(sapply(KMA2014_2015Strata, function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  success <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(!ind == "0") / length(ind)} )
  plot(x = sort(success, decreasing = TRUE), main = mix, xlab = "Sorted Individuatls", col = 8, type = "p", pch = 16, ylim = c(0, 1), ylab = "Success Rate")
  
  cutoff.floor = 0.9
  cutoff <- min(sort(success, decreasing = TRUE)[sort(success, decreasing = TRUE) > (seq(success) / length(success))])
  cutoff <- pmin(cutoff, cutoff.floor)
  points(y = sort(success, decreasing = TRUE)[sort(success, decreasing = TRUE) < cutoff],
         x = seq(success)[sort(success, decreasing = TRUE) < cutoff], pch = 16, col = "black")
  
  abline(h = 0.8, col = "red", lwd = 3)
  abline(h = cutoff.floor, col = "red", lwd = 3, lty = 3)
  segments(x0 = 0, x1 = length(success), y0 = 0, y1 = 1, lwd = 3)
  
  # points(smooth.spline(sort(success, decreasing = TRUE), spar = 0.2), type = "l", col = 1, lwd = 3)
  # 
  # y.loess <- loess(sort(success, decreasing = TRUE) ~ seq(success), span = 0.75, data.frame(x = seq(success), y = sort(success, decreasing = TRUE)))
  # y.predict <- predict(y.loess, data.frame(x = seq(success)))
  # lines(x = seq(success), y = y.predict, col = "black", lwd = 3)
  # 
  # y.loess <- loess(sort(success, decreasing = TRUE) ~ seq(success), span = 0.5, data.frame(x = seq(success), y = sort(success, decreasing = TRUE)))
  # y.predict <- predict(y.loess, data.frame(x = seq(success)))
  # lines(x = seq(success), y = y.predict, col = "blue", lwd = 3)
  # 
  # y.loess <- loess(sort(success, decreasing = TRUE) ~ seq(success), span = 0.3, data.frame(x = seq(success), y = sort(success, decreasing = TRUE)))
  # y.predict <- predict(y.loess, data.frame(x = seq(success)))
  # lines(x = seq(success), y = y.predict, col = "green", lwd = 3)
  # 
  # fit1 <- lm(sort(success, decreasing = TRUE) ~ poly(seq(success) + 100, 1, raw = TRUE))
  # fit2 <- lm(sort(success, decreasing = TRUE) ~ poly(seq(success) + 100, 2, raw = TRUE))
  # fit3 <- lm(sort(success, decreasing = TRUE) ~ poly(seq(success) + 100, 3, raw = TRUE))
  # fit4 <- lm(sort(success, decreasing = TRUE) ~ poly(seq(success) + 100, 4, raw = TRUE))
  # 
  # points(y = predict(fit1, data.frame(x = seq(success) + 100)), x = seq(success), type = "l", col = 1, lwd = 3)
  # points(y = predict(fit2, data.frame(x = seq(success) + 100)), x = seq(success), type = "l", col = 2, lwd = 3)
  # points(y = predict(fit3, data.frame(x = seq(success) + 100)), x = seq(success), type = "l", col = 3, lwd = 3)
  # points(y = predict(fit4, data.frame(x = seq(success) + 100)), x = seq(success), type = "l", col = 4, lwd = 3)

  text(x = 0, y = 0.5, labels = paste("cutoff", sum(success < cutoff), sep = "_"), cex = 1.2, pos = 4)
  text(x = 0, y = 0.4, labels = paste("cutoff.floor", sum(success < cutoff.floor), sep = "_"), pos = 4)
  text(x = 0, y = 0.3, labels = paste("0.8", sum(success < 0.8), sep = "_"), pos = 4)
  
}))


### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_KMA2014_2015Strata_ColSize <- sapply(paste(KMA2014_2015Strata, ".gcl", sep = ''), function(x) get(x)$n)
KMA2014_2015Strata_SampleSizes[, "Genotyped"] <- Original_KMA2014_2015Strata_ColSize


### Alternate
## Indentify alternate species individuals
KMA2014_2015Strata_Alternate <- FindAlternateSpecies.GCL(sillyvec = KMA2014_2015Strata, species = "sockeye")

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport = KMA2014_2015Strata_Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)

## Get number of individuals per silly after removing alternate species individuals
ColSize_KMA2014_2015Strata_PostAlternate <- sapply(paste(KMA2014_2015Strata, ".gcl", sep = ''), function(x) get(x)$n)
KMA2014_2015Strata_SampleSizes[, "Alternate"] <- Original_KMA2014_2015Strata_ColSize-ColSize_KMA2014_2015Strata_PostAlternate


### Missing
## Remove individuals with >20% missing data
KMA2014_2015Strata_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = KMA2014_2015Strata, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_KMA2014_2015Strata_PostMissLoci <- sapply(paste(KMA2014_2015Strata, ".gcl", sep = ''), function(x) get(x)$n)
KMA2014_2015Strata_SampleSizes[, "Missing"] <- ColSize_KMA2014_2015Strata_PostAlternate-ColSize_KMA2014_2015Strata_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
KMA2014_2015Strata_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = KMA2014_2015Strata, loci = loci48, quantile = NULL, minproportion = 0.95)
KMA2014_2015Strata_DuplicateCheckReportSummary <- sapply(KMA2014_2015Strata, function(x) KMA2014_2015Strata_DuplicateCheck95MinProportion[[x]]$report)
KMA2014_2015Strata_DuplicateCheckReportSummary

## Remove duplicate individuals
KMA2014_2015Strata_RemovedDups <- RemoveDups.GCL(KMA2014_2015Strata_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_KMA2014_2015Strata_PostDuplicate <- sapply(paste(KMA2014_2015Strata, ".gcl", sep = ''), function(x) get(x)$n)
KMA2014_2015Strata_SampleSizes[, "Duplicate"] <- ColSize_KMA2014_2015Strata_PostMissLoci-ColSize_KMA2014_2015Strata_PostDuplicate


### Final
KMA2014_2015Strata_SampleSizes[, "Final"] <- ColSize_KMA2014_2015Strata_PostDuplicate
KMA2014_2015Strata_SampleSizes

write.xlsx(KMA2014_2015Strata_SampleSizes, file = "Output/KMA2014_2015Strata_SampleSizes.xlsx")
dput(x = KMA2014_2015Strata_SampleSizes, file = "Objects/KMA2014_2015Strata_SampleSizes.txt")



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
sapply(combined.loci.46, function(loci2combine) {CombineLoci.GCL(sillyvec = KMA2014_2015Strata, markerset = loci2combine, update = TRUE)} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save PostQC/Combined loci .gcl's as back-up:
invisible(sapply(KMA2014_2015Strata, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_Strata_PostQC_CombinedLoci/" , silly, ".txt", sep = ''))} )); beep(8)

## Dput LocusControl
dput(x = LocusControl, file = "Objects/LocusControl50.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = KMA2014_2015Strata, loci = loci48[-mito.loci48], path = "Genepop/KMA2014_2015Strata_46nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/KMA2014_2015Strata_46nuclearloci.txt.P")

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
LocusControl <- dget(file = "Objects/LocusControl50.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("EASSIP-LateAugust2014", "OriginalLocusControl96.txt", "OriginalLocusControl48.txt", "LocusControl98.txt", "LocusControl99.txt", "OLD 15RG")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get post-QC, stratified, combined loci, mixtures
invisible(sapply(KMA2014_2015Strata, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_Strata_PostQC_CombinedLoci/", silly, ".txt", sep = "")), pos = 1)} )); beep(4)
objects(pattern = "\\.gcl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get/Create New 14RG MSA Objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ## Get baseline objects needed for MSA
# setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects")
# 
# KMA473Pops15FlatPrior <- dget(file = "KMA473Pops15FlatPrior.txt")
# KMA473PopsInits <- dget(file = "KMA473PopsInits.txt")
# KMA473PopsGroupVec15 <- dget(file = "KMA473PopsGroupVec15.txt")
# KMA473Pops <- dget(file = "KMA473Pops.txt")
# KMA15GroupsPC <- dget(file = "PCGroups15.txt")
# KMA47346Baseline <- dget(file = "KMA47346Baseline.txt")
# KMA47389Baseline <- dget(file = "KMA47389Baseline.txt")
# WASSIPSockeyeSeeds <- dget(file = "V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt")
# KMA473CommonNames <- dget(file = "CommonNames473.txt")
# KMA15Colors <- dget(file = "Colors15.txt")
# setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# 
# ## Copy these baseline objects and put them in the Mixtures/Objects directory
# 
# setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects")
# file.copy(from = c("KMA473Pops15FlatPrior.txt", "KMA473PopsInits.txt", "KMA473PopsGroupVec14.txt", "KMA473Pops.txt", "PCGroups15.txt", "KMA47346Baseline.txt", "KMA47389Baseline.txt", "CommonNames473.txt", "Colors15.txt",
#                    "V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt"), 
#           to = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects")
# setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# 
# 
# KMA15GroupsPC2Rows <- PCGroups15Rows2 <- c("West of\nChignik", "Black\nLake", "Chignik\nLake", "U. Station\nAkalura", "Frazer\n", "Ayakulik\n", "Karluk\n", "Uganik\n", "Northwest\nKodiak", "Afognak\n", "Eastside\nKodiak", "Saltery\n", "Cook\nInlet", "PWS\n", "South of\nCape Suckling")
# dput(x = KMA15GroupsPC2Rows, file = "Objects/KMA15GroupsPC2Rows.txt")
# # Note, I changed the names for CommonNames473.txt -> KMA473CommonNames.txt & PCGroups15 -> KMA15GroupsPC.txt


## Create baseline objects needed for MSA
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects")

KMA473PopsGroupVec14 <- dget(file = "KMA473PopsGroupVec14.txt")
KMA473Pops14FlatPrior <- Prior.GCL(groupvec = KMA473PopsGroupVec14, groupweights = rep(1 / 14, 14), minval = 0.01)
  dput(x = KMA473Pops14FlatPrior, file = "KMA473Pops14FlatPrior.txt")
KMA473PopsInits <- dget(file = "KMA473PopsInits.txt")
KMA14GroupsPC <- dget(file = "KMA14GroupsPC.txt")
KMA14GroupsPC2Rows <- c("West of\nChignik", "Black\nLake", "Chignik\nLake", "U. Station\nAkalura", 
                        "Frazer\nAyakulik", "Karluk\n", "Uganik\n", "Northwest\nKodiak", 
                        "Afognak\n", "Eastside\nKodiak", "Saltery\n", "Cook\nInlet", 
                        "PWS\n", "South of\nCape Suckling")
  dput(x = KMA14GroupsPC2Rows, file = "KMA14GroupsPC2Rows.txt")
KMA47346Baseline <- dget(file = "KMA47346Baseline.txt")
WASSIPSockeyeSeeds <- dget("V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt")
KMA473CommonNames <- dget(file = "CommonNames473.txt")
KMA14Colors <- c(KMA15Colors[1:4], "dodgerblue", KMA15Colors[7:15])
  dput(x = KMA14Colors, file = "KMA14Colors.txt")
KMA14ColorsRGB <- t(col2rgb(KMA14Colors))
rownames(KMA14ColorsRGB) <- KMA14GroupsPC
  dput(x = KMA14ColorsRGB, file = "KMA14ColorsRGB.txt")

  

setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects")
dput(x = KMA473PopsGroupVec14, file = "KMA473PopsGroupVec14.txt")
dput(x = KMA473Pops14FlatPrior, file = "KMA473Pops14FlatPrior.txt")
dput(x = KMA14GroupsPC, file = "KMA14GroupsPC.txt")
dput(x = KMA14GroupsPC2Rows, file = "KMA14GroupsPC2Rows.txt")
dput(x = KMA14Colors, file = "KMA14Colors.txt")
dput(x = KMA14ColorsRGB, file = "KMA14ColorsRGB.txt")

setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")



dir.create(path = "BAYES/2014-2015 Mixtures 46loci 14RG")
sapply(c("Control", "Mixture", "Output"), function(folder) {dir.create(path = paste(getwd(), "BAYES/2014-2015 Mixtures 46loci 14RG", folder, sep = "/"))} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA files for BAYES 2014 Early Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Starting with a Regionally flat prior, then rolling prior afterwards
KMA473Pops14FlatPrior


KMA2014Strata_1_Early <- grep(pattern = "1_Early", x = KMA2014Strata, value = TRUE)
Round1Mixtures_2014 <- c(KMA2014Strata_1_Early, "SALITC14_2_Middle")
dput(x = Round1Mixtures_2014, file = "Objects/Round1Mixtures_2014.txt")

## Dumping Mixture files
KMA47346MixtureFormat <- CreateMixture.GCL(sillys = KMA2014Strata_1_Early[1], loci = loci46, IDs = NULL, mixname = KMA2014Strata_1_Early[1],
                                           dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)
dput(KMA47346MixtureFormat, file = "Objects/KMA47346MixtureFormat.txt")

sapply(Round1Mixtures_2014, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round1Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = KMA473Pops14FlatPrior, initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round1Mixtures_2014, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
                                                              maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
                                                              mixvec = Round1Mixtures_2014, prior = "",  
                                                              ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round1Mixtures_2014_Estimates, file = "Estimates objects/Round1Mixtures_2014_Estimates.txt")
dput(Round1Mixtures_2014_Estimates$Stats, file = "Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")

Round1Mixtures_2014_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2014_Estimates.txt")
Round1Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round1Mixtures_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round1Mixtures_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})
sapply(Round1Mixtures_2014[c(3, 4, 2, 1, 5)], function(Mix) {
  BarPlot <- barplot2(Round1Mixtures_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round1Mixtures_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", type = "h", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
  })

# Quick look at raw posterior output
str(Round1Mixtures_2014_Estimates$Output)
Round1Mixtures_2014_Header <- setNames(object = c("Ayakulik Section Early June 1-27, 2014",
                                                  "Karluk Section Early June 1-27, 2014",
                                                  "Uganik Section Early June 1-27, 2014",
                                                  "Uyak Section Early June 1-27, 2014",
                                                  "Alitak Section Middle June 28-July 25, 2014"), 
                                       nm = Round1Mixtures_2014)
dput(x = Round1Mixtures_2014_Header, file = "Objects/Round1Mixtures_2014_Header.txt")


PlotPosterior <- function(mixvec = NULL, output, header = NULL, groups, colors = NULL, set.mfrow, thin = 10, chains = 5){
  if(is.null(colors)) {colors <- rep("black", length(groups))}
  if(is.null(mixvec)) {mixvec <- names(output)}
  
  par(mfrow = set.mfrow, mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 1.1))
  
  invisible(sapply(mixvec, function(Mix) {
    invisible(sapply(seq(groups), function(i) {
      RG <- output[[Mix]][, i]
      RG <- RG[seq(1, length(RG), thin)]
      plot(RG, type = "l", ylim = c(0,1), xlab = "", ylab = "")
      abline(v = seq(0, length(RG), length(RG)/chains), xpd = FALSE)
      text(x = length(RG)/2, y = 0.96, labels = groups[i], col = colors[i], cex = 1.2, font = 2)} ))
    if(prod(set.mfrow) - length(groups) == 1) {plot.new()}
    if(prod(set.mfrow) - length(groups) == 2) {plot.new(); plot.new()}
    if(!is.null(header)) {
      mtext(text = header[Mix], side = 3, outer = TRUE, cex = 1.5)
    } else {
      mtext(text = Mix, side = 3, outer = TRUE, cex = 1.5)
    }
    mtext(text = paste("Iteration (", chains, " chain[s] each)", sep = ''), side = 1, outer = TRUE, cex = 1.5, line = 1.5)
    mtext(text = "Posterior", side = 2, outer = TRUE, cex = 1.5, line = 1.5)
  }))
  
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
}
dput(x = PlotPosterior, file = "Objects/PlotPosterior.txt")

PlotPosterior(mixvec = Round1Mixtures_2014[c(3, 4, 2, 1, 5)], output = Round1Mixtures_2014_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round1Mixtures_2014_Header, set.mfrow = c(5, 3), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot <- function(mixvec, estimatesstats, groups, groups2rows = NULL, colors, header) {
  
  if(is.null(groups2rows)) {groups2rows <- groups}
  if(length(mixvec) != length(header)) {stop("Header is not the same length as mixvec!")}
  for(i in seq(mixvec)) {
    header[i] <- paste(header[i], "\nn = ", get(paste(mixvec[i], ".gcl", sep = ''))$n, sep = '')
  }
  if("Stats" %in% names(estimatesstats)) {estimatesstats <- estimatesstats$Stats}
  
  par(mfrow = c(1, 1), mar = c(5.6, 4.6, 5.1, 2.1))
  sapply(mixvec, function(Mix) {
    plot(estimatesstats[[Mix]][, "median"], cex = 3, pch = 16, col = colors, ylab = "Proportion of Mixture", 
         ylim = c(0, 1), xlab = "", axes = FALSE, main = header[Mix], cex.main = 2, cex.lab = 1.5)
    arrows(x0 = seq(groups), y0 = estimatesstats[[Mix]][, "5%"], x1 = seq(groups), 
           y1 = estimatesstats[[Mix]][, "95%"], angle = 90, code = 3, length = 0.2, lwd = 2)
    points(estimatesstats[[Mix]][, "median"], cex = 3, pch = 16, col = colors)
    axis(side = 2)
    axis(side = 1, labels = NA, at = seq(groups))
    text(x = (seq(groups)) - 0.35, y = 0, labels = groups2rows, srt = 60, pos = 1, offset = 3.5, xpd = TRUE)
  })
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}
dput(x = QuickPlot, file = "Objects/QuickPlot.txt")

QuickPlot(mixvec = Round1Mixtures_2014, estimatesstats = Round1Mixtures_2014_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round1Mixtures_2014_Header)


## Barplots
QuickBarplot <- function(mixvec, estimatesstats, groups, groups2rows = NULL, header) {
  while(!require(gplots, quietly = TRUE)){install.packages("gplots")}
  
  # if(is.null(groups2rows)) {groups2rows <- groups}
  if(length(mixvec) != length(header)) {stop("Header is not the same length as mixvec!")}
  for(i in seq(mixvec)) {
    header[i] <- paste(header[i], "\nn = ", get(paste(mixvec[i], ".gcl", sep = ''))$n, sep = '')
  }
  if("Stats" %in% names(estimatesstats)) {estimatesstats <- estimatesstats$Stats}
  
  par(mfrow = c(1, 1), mar = c(7.1, 5.1, 4.1, 2.1), oma = rep(0, 4))
  sapply(mixvec, function(Mix) {
    Barplot <- barplot2(height = estimatesstats[[Mix]][, "median"] * 100,
                        beside = TRUE, plot.ci = TRUE, ci.lwd = 1,
                        ci.l = estimatesstats[[Mix]][, "5%"] * 100,
                        ci.u = estimatesstats[[Mix]][, "95%"] * 100,
                        ylim = c(0, 100), col = "blue", yaxt = 'n', xaxt = 'n',
                        main = header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
    axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
    abline(h = 0, xpd = FALSE)
    
    if(is.null(groups2rows)) {
      text(x = Barplot[, 1], y = -1, labels = groups, srt = 90, adj =  1, xpd = TRUE, cex = 0.7)
    } else {
      mtext(text = groups2rows, side = 1, line = 1, at = Barplot[, 1], adj = 0.5, cex = 0.7)
    }
  })
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}
dput(x = QuickBarplot, file = "Objects/QuickBarplot.txt")
QuickBarplot(mixvec = Round1Mixtures_2014, estimatesstats = Round1Mixtures_2014_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round1Mixtures_2014_Header)


## Make violin plots of posteriors with RGs sorted

ViolinPlot <- function(mixvec = NULL, estimates, groups, colors, header, wex = 1, thin = 10) {
  while(!require(vioplot, quietly = TRUE)){install.packages("vioplot")}
  
  if(is.null(mixvec)) {mixvec <- names(estimates$Stats)}
  
  par(mar = c(5.6, 4.6, 3.6, 1.1))
  sapply(mixvec, function(Mix) {
    plot(estimates$Stats[[Mix]][, "median"], cex = 3, pch = 16, col = colors, ylab = "Proportion of Mixture", ylim = c(0, 1), xlab = "", axes = FALSE, main = header[[Mix]], cex.main = 2, cex.lab = 1.5)
    sapply(seq(groups), function(i) {vioplot2(estimates$Output[[Mix]][seq(from = 1, to = nrow(estimates$Output[[Mix]]), by = thin), i], at = i, horizontal = FALSE, col = colors[i], border = TRUE, drawRect = FALSE, rectCol = colors[i], add = TRUE, wex = wex, lwd = 2)})
    arrows(x0 = seq(groups), y0 = estimates$Stats[[Mix]][, "5%"], x1 = seq(groups), y1 = estimates$Stats[[Mix]][, "95%"], angle = 90, code = 3, length = 0.2, lwd = 2)
    points(estimates$Stats[[Mix]][, "median"], cex = 2, pch = 21, col = "white", bg = colors, lwd = 3)
    axis(side = 2, lwd = 3, cex.axis = 1.5)
    text(x = (seq(groups)) - 0.35, y = 0, labels = groups, srt = 60, pos = 1, offset = 2.5, xpd = TRUE)
    axis(side = 1, labels = NA, at = seq(groups), pos = 0, lwd = 2, tick = FALSE)
    abline(h = 0, lwd = 3, xpd = FALSE)
  } )
  
}
dput(x = ViolinPlot, file = "Objects/ViolinPlot.txt")


ViolinPlot(estimates = Round1Mixtures_2014_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round1Mixtures_2014_Header)
rm(Round1Mixtures_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 1 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA473PopsGroupVec31 <- as.numeric(readClipboard())
dput(x = KMA473PopsGroupVec31, file = "Objects/KMA473PopsGroupVec31.txt")
KMA31GroupsPC <- readClipboard()
dput(x = KMA31GroupsPC, file = "Objects/KMA31GroupsPC.txt")

Round1Mixtures_2014_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/40K Iterations", 
                                                                   mixvec = Round1Mixtures_2014, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE)
str(Round1Mixtures_2014_31RG_Estimates)
QuickBarplot(mixvec = Round1Mixtures_2014[c(3, 4, 2, 1, 5)], estimatesstats = Round1Mixtures_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round1Mixtures_2014_Header[c(3, 4, 2, 1, 5)])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA files for BAYES 2014 Middle Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Indetify Strata to Run
KMA2014Strata_2_Middle <- grep(pattern = "2_Middle", x = KMA2014Strata, value = TRUE)
Round2Mixtures_2014 <- KMA2014Strata_2_Middle[-which(KMA2014Strata_2_Middle == "SALITC14_2_Middle")]  # Remove SALITC14_2_Middle as I want to run that out further
dput(x = Round2Mixtures_2014, file = "Objects/Round2Mixtures_2014.txt")

# Create rolling prior based on Round 1 estimates
Round2Mixtures_2014_Prior <- sapply(Round1Mixtures_2014_EstimatesStats[1:4], function(Mix) {Prior.GCL(groupvec = KMA473PopsGroupVec14, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round2Mixtures_2014_Prior) <- gsub(pattern = "1_Early", replacement = "2_Middle", x = names(Round2Mixtures_2014_Prior))  # This changes the names
dput(x = Round2Mixtures_2014_Prior, file = "Objects/Round2Mixtures_2014_Prior.txt")
str(Round2Mixtures_2014_Prior)


## Dumping Mixture files
sapply(Round2Mixtures_2014, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round2Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round2Mixtures_2014_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2014, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
                                                              maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
                                                              mixvec = Round2Mixtures_2014, prior = "",  
                                                              ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures_2014_Estimates, file = "Estimates objects/Round2Mixtures_2014_Estimates.txt")
dput(Round2Mixtures_2014_Estimates$Stats, file = "Estimates objects/Round2Mixtures_2014_EstimatesStats.txt")

Round2Mixtures_2014_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2014_Estimates.txt")
Round2Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2014_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SAYAKC14_2_Middle and SUYAKC14_2_
sapply(Round2Mixtures_2014[c(3, 4, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round2Mixtures_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures_2014_Estimates$Output)
Round2Mixtures_2014_Header <- setNames(object = c("Ayakulik Section Middle June 28-July 25, 2014",
                                                  "Karluk Section Middle June 28-July 25, 2014",
                                                  "Uganik Section Middle June 28-July 25, 2014",
                                                  "Uyak Section Middle June 28-July 25, 2014"), 
                                       nm = Round2Mixtures_2014)
dput(x = Round2Mixtures_2014_Header, file = "Objects/Round2Mixtures_2014_Header.txt")

PlotPosterior(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], output = Round2Mixtures_2014_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round2Mixtures_2014_Header, set.mfrow = c(5, 3), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round2Mixtures_2014_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round2Mixtures_2014_Header)

## Quick barplot
QuickBarplot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round2Mixtures_2014_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round2Mixtures_2014_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimates = Round2Mixtures_2014_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round2Mixtures_2014_Header)
rm(Round2Mixtures_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 Repeated MSA files for BAYES 80K Iterations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping out a new control file to re-analyze all Round 2 mixtures with a 80K posterior as opposed to our standard 40K
## Perhaps it will converge on one of the modes?

# Just use original mixture, no need to create a new one!

## Dumping a New Control file that goes out further!

sapply(Round2Mixtures_2014[c(1, 4)], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round2Mixtures_2014_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2014[c(1, 4)], function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})
# Will rename later, but leave as is for now for BayesCopyPaste.GCL


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures80K_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
                                                                 maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/80K Iterations", 
                                                                 mixvec = Round2Mixtures_2014[c(1, 4)], prior = "",  
                                                                 ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures80K_2014_Estimates, file = "Estimates objects/Round2Mixtures80K_2014_Estimates.txt")
dput(Round2Mixtures80K_2014_Estimates$Stats, file = "Estimates objects/Round2Mixtures80K_2014_EstimatesStats.txt")

Round2Mixtures80K_2014_Estimates <- dget(file = "Estimates objects/Round2Mixtures80K_2014_Estimates.txt")
Round2Mixtures80K_2014_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures80K_2014_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures80K_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures80K_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUYAKC14_2_
sapply(Round2Mixtures_2014[c(4, 1)], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round2Mixtures80K_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures80K_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures80K_2014_Estimates$Output)

PlotPosterior(mixvec = Round2Mixtures_2014[c(4, 1)], output = Round2Mixtures80K_2014_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round2Mixtures_2014_Header, set.mfrow = c(5, 3), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round2Mixtures_2014[c(4, 1)], estimatesstats = Round2Mixtures80K_2014_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round2Mixtures_2014_Header[c(4, 1)])

## Quick barplot
QuickBarplot(mixvec = Round2Mixtures_2014[c(4, 1)], estimatesstats = Round2Mixtures80K_2014_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round2Mixtures_2014_Header[c(4, 1)])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round2Mixtures_2014[c(4, 1)], estimates = Round2Mixtures80K_2014_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round2Mixtures_2014_Header[c(4, 1)])
rm(Round2Mixtures80K_2014_Estimates)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round2Mixtures_2014_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2014_Estimates.txt")
Round2Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2014_EstimatesStats.txt")


sapply(Round2Mixtures_2014[c(4, 1)], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round2Mixtures_2014_EstimatesStats[[Mix]][, "median"],
                                                 Round2Mixtures80K_2014_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round2Mixtures_2014_EstimatesStats[[Mix]][, "5%"],
                                               Round2Mixtures80K_2014_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round2Mixtures_2014_EstimatesStats[[Mix]][, "95%"],
                                               Round2Mixtures80K_2014_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round2Mixtures_2014_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA14GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 2 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SKARLC14_2_Middle and SUGANC14_2_Middle as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from SAYAKC14_2_Middle and SUYAKC14_2_Middle

sapply(Round2Mixtures_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SAYAKC14_2_Middle and SUYAKC14_2_Middle
sapply(Round2Mixtures80K_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUYAKC14_2_Middle still an issue, but going to live with it

round(cbind("40K" = Round2Mixtures_2014_Estimates$Stats$SUYAKC14_2_Middle[, "GR"],
            "80K" = Round2Mixtures80K_2014_Estimates$Stats$SUYAKC14_2_Middle[, "GR"]), 2)

str(Round2Mixtures_2014_Estimates)
Round2Mixtures_2014

Round2Mixtures_2014_Estimates_Final <- list(Stats = c(Round2Mixtures80K_2014_Estimates$Stats["SAYAKC14_2_Middle"],
                                                         Round2Mixtures_2014_Estimates$Stats[c("SKARLC14_2_Middle", "SUGANC14_2_Middle")],
                                                         Round2Mixtures80K_2014_Estimates$Stats["SUYAKC14_2_Middle"]),
                                            Output = c(Round2Mixtures80K_2014_Estimates$Output["SAYAKC14_2_Middle"],
                                                       Round2Mixtures_2014_Estimates$Output[c("SKARLC14_2_Middle", "SUGANC14_2_Middle")],
                                                       Round2Mixtures80K_2014_Estimates$Output["SUYAKC14_2_Middle"]))
str(Round2Mixtures_2014_Estimates_Final)
dput(x = Round2Mixtures_2014_Estimates_Final, file = "Estimates objects/Round2Mixtures_2014_Estimates_Final.txt")
dput(x = Round2Mixtures_2014_Estimates_Final$Stats, file = "Estimates objects/Round2Mixtures_2014_EstimatesStats_Final.txt")

Round2Mixtures_2014_Estimates_Final <- dget(file = "Estimates objects/Round2Mixtures_2014_Estimates_Final.txt")
Round2Mixtures_2014_EstimatesStats_Final <- dget(file = "Estimates objects/Round2Mixtures_2014_EstimatesStats_Final.txt")

sapply(Round2Mixtures_2014_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # South of Cape Suckling out for SUYAKC14_2_Middle


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 2 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2014_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/80K Iterations", 
                                                                   mixvec = Round2Mixtures_2014, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round2Mixtures_2014_31RG_Estimates)
QuickBarplot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round2Mixtures_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round2Mixtures_2014_Header[c(3, 4, 2, 1)])
QuickBarplot(mixvec = Round1Mixtures_2014[5], estimatesstats = Round1Mixtures_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round1Mixtures_2014_Header[5])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 MSA files for BAYES 2014 Late Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify Strata to Run
KMA2014Strata_3_Late <- grep(pattern = "3_Late", x = KMA2014Strata, value = TRUE)
Round3Mixtures_2014 <- KMA2014Strata_3_Late
dput(x = Round3Mixtures_2014, file = "Objects/Round3Mixtures_2014.txt")

# Create rolling prior based on Round 2 estimates
Round1Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")

Round3Mixtures_2014_Prior <- sapply(c(Round1Mixtures_2014_EstimatesStats["SALITC14_2_Middle"], Round2Mixtures_2014_EstimatesStats_Final[1:4]), function(Mix) {Prior.GCL(groupvec = KMA473PopsGroupVec14, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
str(Round3Mixtures_2014_Prior)
names(Round3Mixtures_2014_Prior) <- gsub(pattern = "2_Middle", replacement = "3_Late", x = names(Round3Mixtures_2014_Prior))  # This changes the names
dput(x = Round3Mixtures_2014_Prior, file = "Objects/Round3Mixtures_2014_Prior.txt")
str(Round3Mixtures_2014_Prior)


## Dumping Mixture files
sapply(Round3Mixtures_2014, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round3Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round3Mixtures_2014_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round3Mixtures_2014, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create Control file for SUGANC14_3_Late 80K as it is known offender from 15RG
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sapply(Round3Mixtures_2014[4], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round3Mixtures_2014_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

sapply(Round3Mixtures_2014[4], function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci 14RG/Output/80K Iterations/", Mix, sep = ""))})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 3 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2014_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/40K Iterations", 
  mixvec = Round3Mixtures_2014, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round3Mixtures_2014_Estimates, file = "Estimates objects/Round3Mixtures_2014_Estimates.txt")
dput(Round3Mixtures_2014_Estimates$Stats, file = "Estimates objects/Round3Mixtures_2014_EstimatesStats.txt")

Round3Mixtures_2014_Estimates <- dget(file = "Estimates objects/Round3Mixtures_2014_Estimates.txt")
Round3Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures_2014_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round3Mixtures_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round3Mixtures_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # No Issues!!!
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round3Mixtures_2014[c(4, 5, 3, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round3Mixtures_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round3Mixtures_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round3Mixtures_2014_Estimates$Output)
Round3Mixtures_2014_Header <- setNames(object = c("Alitak Section Middle July 26-August 29, 2014",
                                                  "Ayakulik Section Middle July 26-August 29, 2014",
                                                  "Karluk Section Middle July 26-August 25, 2014",
                                                  "Uganik Section Middle July 26-August 25, 2014",
                                                  "Uyak Section Middle July 26-August 25, 2014"), 
                                       nm = Round3Mixtures_2014)
dput(x = Round3Mixtures_2014_Header, file = "Objects/Round3Mixtures_2014_Header.txt")

PlotPosterior(mixvec = Round3Mixtures_2014[c(4, 5, 3, 2, 1)], output = Round3Mixtures_2014_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round3Mixtures_2014_Header, set.mfrow = c(5, 3), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 3 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round3Mixtures_2014[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2014_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round3Mixtures_2014_Header)

## Quick barplot
QuickBarplot(mixvec = Round3Mixtures_2014[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2014_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round3Mixtures_2014_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round3Mixtures_2014[c(4, 5, 3, 2, 1)], estimates = Round3Mixtures_2014_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round3Mixtures_2014_Header)
rm(Round3Mixtures_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 3 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2014_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/40K Iterations", 
                                                                   mixvec = Round3Mixtures_2014, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round3Mixtures_2014_31RG_Estimates)
QuickBarplot(mixvec = Round3Mixtures_2014[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round3Mixtures_2014_Header[c(4, 5, 3, 2, 1)])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 4 MSA files for BAYES 2014 LateLate ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Indetify Strata to Run
KMA2014Strata_4_LateLate <- grep(pattern = "4_LateLate", x = KMA2014Strata, value = TRUE)
Round4Mixtures_2014 <- c(KMA2014Strata_4_LateLate)
dput(x = Round4Mixtures_2014, file = "Objects/Round4Mixtures_2014.txt")

# Create rolling prior based on Round 3 estimates
str(Round3Mixtures_2014_EstimatesStats)
sapply(Round3Mixtures_2014_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})

str(Round3Mixtures_2014_EstimatesStats[c("SKARLC14_3_Late", "SUGANC14_3_Late", "SUYAKC14_3_Late")])

Round4Mixtures_2014_Prior <- sapply(Round3Mixtures_2014_EstimatesStats[c("SKARLC14_3_Late", "SUGANC14_3_Late", "SUYAKC14_3_Late")],
                                    function(Mix) {Prior.GCL(groupvec = KMA473PopsGroupVec14, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round4Mixtures_2014_Prior) <- gsub(pattern = "3_Late", replacement = "4_LateLate", x = names(Round4Mixtures_2014_Prior))  # This changes the names
dput(x = Round4Mixtures_2014_Prior, file = "Objects/Round4Mixtures_2014_Prior.txt")
str(Round4Mixtures_2014_Prior)


## Dumping Mixture files
sapply(Round4Mixtures_2014, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round4Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round4Mixtures_2014_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round4Mixtures_2014, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 4 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round4Mixtures_2014_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
  mixvec = Round4Mixtures_2014, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round4Mixtures_2014_Estimates, file = "Estimates objects/Round4Mixtures_2014_Estimates.txt")
dput(Round4Mixtures_2014_Estimates$Stats, file = "Estimates objects/Round4Mixtures_2014_EstimatesStats.txt")

Round4Mixtures_2014_Estimates <- dget(file = "Estimates objects/Round4Mixtures_2014_Estimates.txt")
Round4Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round4Mixtures_2014_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round4Mixtures_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round4Mixtures_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # No issues
sapply(Round4Mixtures_2014[c(2, 3, 1)], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round4Mixtures_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round4Mixtures_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round4Mixtures_2014_Estimates$Output)
Round4Mixtures_2014_Header <- setNames(object = c("Karluk Section LateLate August 25-29, 2014",
                                                  "Uganik Section LateLate August 25-29, 2014",
                                                  "Uyak Section LateLate August 25-29, 2014"), 
                                       nm = Round4Mixtures_2014)
dput(x = Round4Mixtures_2014_Header, file = "Objects/Round4Mixtures_2014_Header.txt")

PlotPosterior(mixvec = Round4Mixtures_2014, output = Round4Mixtures_2014_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round4Mixtures_2014_Header, set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 4 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round4Mixtures_2014[c(2, 3, 1)], estimatesstats = Round4Mixtures_2014_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round4Mixtures_2014_Header)

## Quick barplot
QuickBarplot(mixvec = Round4Mixtures_2014[c(2, 3, 1)], estimatesstats = Round4Mixtures_2014_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round4Mixtures_2014_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round4Mixtures_2014[c(2, 3, 1)], estimates = Round4Mixtures_2014_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round4Mixtures_2014_Header)
rm(Round4Mixtures_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 4 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round4Mixtures_2014_31RG_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/40K Iterations", 
  mixvec = Round4Mixtures_2014, prior = "",  
  ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round4Mixtures_2014_31RG_Estimates)
QuickBarplot(mixvec = Round4Mixtures_2014[c(2, 3, 1)], estimatesstats = Round4Mixtures_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round4Mixtures_2014_Header[c(2, 3, 1)])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Estimator for Uganik, Uyak and Karluk 3_Late + 4_LateLate ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Using KMA Sockeye Harvest data to apply stratified estimator to combine the
# 380 fish in 3_Late with the 285 fish in 4_LateLate to get at a "final"
# estimate for the "Late" strata from July 26 - August 29, 2014

HarvestByStrata2014 <- read.table(file = "Harvest/2014HarvestByStrata.txt", header = TRUE, sep = "\t", as.is = TRUE)
HarvestByStrata2014.mat <- data.matrix(HarvestByStrata2014[, -1])
dimnames(HarvestByStrata2014.mat) <- list(HarvestByStrata2014$location, c("1_Early", "2_Middle", "3_Late", "4_LateLate"))
dput(x = HarvestByStrata2014.mat, file = "Objects/HarvestByStrata2014.txt"); rm(HarvestByStrata2014.mat)
HarvestByStrata2014 <- dget(file = "Objects/HarvestByStrata2014.txt")

HarvestByStrata2015 <- read.table(file = "Harvest/2015HarvestByStrata.txt", header = TRUE, sep = "\t", as.is = TRUE)
HarvestByStrata2015.mat <- data.matrix(HarvestByStrata2015[, -1])
dimnames(HarvestByStrata2015.mat) <- list(HarvestByStrata2015$location, c("1_Early", "2_Middle", "3_Late"))
dput(x = HarvestByStrata2015.mat, file = "Objects/HarvestByStrata2015.txt"); rm(HarvestByStrata2015.mat)
HarvestByStrata2015 <- dget(file = "Objects/HarvestByStrata2015.txt")

HarvestByStrata2016 <- read.table(file = "Harvest/2016HarvestByStrata.txt", header = TRUE, sep = "\t", as.is = TRUE)
HarvestByStrata2016.mat <- data.matrix(HarvestByStrata2016[, -1])
dimnames(HarvestByStrata2016.mat) <- list(HarvestByStrata2016$location, c("1_Early", "2_Middle", "3_Late"))
dput(x = HarvestByStrata2016.mat, file = "Objects/HarvestByStrata2016.txt"); rm(HarvestByStrata2016.mat)
HarvestByStrata2016 <- dget(file = "Objects/HarvestByStrata2016.txt")


max(HarvestByStrata2014, HarvestByStrata2015, na.rm = TRUE)


# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
# Need to test first to make sure it will work with mixtures of different length (Uganik 3 had to be run out to 80K)
sapply(rownames(HarvestByStrata2014)[!is.na(HarvestByStrata2014[, 4])], function(geomix) {
  assign(x = paste(geomix, "_3_Late_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
           mixvec = c(paste(geomix, "_3_Late", sep = ''), paste(geomix, "_4_LateLate", sep = '')), 
           catchvec = HarvestByStrata2014[geomix, c("3_Late", "4_LateLate")], newname = paste(geomix, "_3_Late_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_3_Late_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_3_Late_Stratified.txt", sep = ''))
} ); beep(5)

str(SKARLC14_3_Late_Stratified)

dimnames(SKARLC14_3_Late_Stratified$Stats)
dimnames(Round4Mixtures_2014_Estimates$Stats$SKARLC14_4_LateLate)

# Look at Late and LateLate separately and then compare to stratified estimate
Round1Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")
Round2Mixtures_2014_EstimatesStats_Final <- dget("Estimates objects/Round2Mixtures_2014_EstimatesStats_Final.txt")
Round3Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round3Mixtures_2014_EstimatesStats.txt")
Round4Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round4Mixtures_2014_EstimatesStats.txt")

KMA2014Strata_EstimatesStats <- c(Round1Mixtures_2014_EstimatesStats, 
                                  Round2Mixtures_2014_EstimatesStats_Final, 
                                  Round3Mixtures_2014_EstimatesStats, 
                                  Round4Mixtures_2014_EstimatesStats)

TempMix <- c("_1_Early", "_2_Middle", "_3_Late", "_4_LateLate")
TempLegend14 <- c("June 1-27", "June 28-July 25", "July 26-August 25", "August 26-29")
GeoMix <- c("SKARLC", "SUGANC", "SUYAKC")
GeoHeader <- setNames(object = c(paste("Karluk/Sturgeon 255-10", "\u2013", "255-20; 256-40", sep = ''),
                                 paste("Uganik/Kupreanof 253", sep = ''),
                                 paste("Uyak Bay 254", sep = '')),
                      nm = GeoMix)
darkcolor <- "blue"
Estimates <- KMA2014Strata_EstimatesStats
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.5
cex.xaxis <- 0.7
cex.main <- 2
cex.leg <- 1.5
ci.lwd <- 2.5

par(mar = c(3.6, 5.1, 3.1, 1.1))
sapply(GeoMix, function(geomix) {
  Barplot1 <- barplot2(height = t(sapply(TempMix, function(tempmix) {Estimates[[paste(geomix, "14", tempmix, sep = '')]][, "median"]})) * 100, 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix, function(tempmix) {Estimates[[paste(geomix, "14", tempmix, sep = '')]][, "5%"]})) * 100, 
                       ci.u = t(sapply(TempMix, function(tempmix) {Estimates[[paste(geomix, "14", tempmix, sep = '')]][, "95%"]})) * 100, 
                       ylim = c(0, 100), col = colorpanel(length(TempMix), low = darkcolor, high = "white"), cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14, x = "topleft", fill = colorpanel(length(TempMix), low = darkcolor, high = "white"), border = "black", bty = "n", cex = cex.leg, title="2014")
  abline(h = 0)
  mtext(text = "Percentage of Catch", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[geomix], side = 3, cex = cex.main)
})

# Harvest percent of Late Stratified from Late and LateLate
round(HarvestByStrata2014[3:5, 3:4] / apply(HarvestByStrata2014[3:5, 3:4], 1, sum) * 100, 1)

sapply(GeoMix, function(geomix) {
  round(get(paste(geomix, "14_3_Late_Stratified", sep = ''))$Stats[, "median"] * 100, 1)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Final Estimates Objects for Each Temporal Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Round1Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")
Round2Mixtures_2014_EstimatesStats_Final <- dget("Estimates objects/Round2Mixtures_2014_EstimatesStats_Final.txt")
Round3Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round3Mixtures_2014_EstimatesStats.txt")
Round4Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round4Mixtures_2014_EstimatesStats.txt")
SKARLC14_3_Late_Stratified <- dget("Estimates objects/SKARLC14_3_Late_Stratified.txt")
SUGANC14_3_Late_Stratified <- dget("Estimates objects/SUGANC14_3_Late_Stratified.txt")
SUYAKC14_3_Late_Stratified <- dget("Estimates objects/SUYAKC14_3_Late_Stratified.txt")

# View as tables by year
require(reshape)
samp.df.2014 <- data.frame(t(sapply(KMA2014Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2014, location ~ temporal.strata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1_Early
str(Round1Mixtures_2014_EstimatesStats)
KMA2014Strata_1_Early_EstimatesStats <- Round1Mixtures_2014_EstimatesStats[1:4]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2_Middle
str(Round2Mixtures_2014_EstimatesStats_Final)
KMA2014Strata_2_Middle_EstimatesStats <- c(Round1Mixtures_2014_EstimatesStats[5], 
                                           Round2Mixtures_2014_EstimatesStats_Final)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3_Late
str(Round3Mixtures_2014_EstimatesStats)
str(Round4Mixtures_2014_EstimatesStats)
str(SKARLC14_3_Late_Stratified$Stats)

# dealing with different output from CustomCombineBAYESOutput and StratifiedEstimator
# colnames(SKARLC14_3_Late_Stratified$Summary) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")
# colnames(SUGANC14_3_Late_Stratified$Summary) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")
# colnames(SUYAKC14_3_Late_Stratified$Summary) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")
# 
# SKARLC14_3_Late_Stratified$Summary <- SKARLC14_3_Late_Stratified$Summary[, c("mean", "sd", "median", "5%", "95%", "P=0", "GR")]
# SUGANC14_3_Late_Stratified$Summary <- SUGANC14_3_Late_Stratified$Summary[, c("mean", "sd", "median", "5%", "95%", "P=0", "GR")]
# SUYAKC14_3_Late_Stratified$Summary <- SUYAKC14_3_Late_Stratified$Summary[, c("mean", "sd", "median", "5%", "95%", "P=0", "GR")]

KMA2014Strata_3_Late_EstimatesStats <- c(Round3Mixtures_2014_EstimatesStats[1:2],
                                         list("SKARLC14_3_Late" = SKARLC14_3_Late_Stratified$Stats),
                                         list("SUGANC14_3_Late" = SUGANC14_3_Late_Stratified$Stats),
                                         list("SUYAKC14_3_Late" = SUYAKC14_3_Late_Stratified$Stats))

str(KMA2014Strata_1_Early_EstimatesStats)
str(KMA2014Strata_2_Middle_EstimatesStats)
str(KMA2014Strata_3_Late_EstimatesStats)


dir.create("Estimates objects/Final")
dput(x = KMA2014Strata_1_Early_EstimatesStats, file = "Estimates objects/Final/KMA2014Strata_1_Early_EstimatesStats.txt")
dput(x = KMA2014Strata_2_Middle_EstimatesStats, file = "Estimates objects/Final/KMA2014Strata_2_Middle_EstimatesStats.txt")
dput(x = KMA2014Strata_3_Late_EstimatesStats, file = "Estimates objects/Final/KMA2014Strata_3_Late_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock composition results 2014 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA2014Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_1_Early_EstimatesStats.txt")
KMA2014Strata_2_Middle_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_2_Middle_EstimatesStats.txt")
KMA2014Strata_3_Late_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_3_Late_EstimatesStats.txt")

KMA2014Strata_EstimatesStats <- c(KMA2014Strata_1_Early_EstimatesStats, 
                                  KMA2014Strata_2_Middle_EstimatesStats, 
                                  KMA2014Strata_3_Late_EstimatesStats)
dput(x = KMA2014Strata_EstimatesStats, file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")
KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")

str(KMA2014Strata_EstimatesStats)


TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_EstimatesStats), value = TRUE)} )

Legend14 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend14 <- sapply(KMA2014, function(geo) {
  Legend14[sapply(TempMix14[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
} )

GeoHeader <- setNames(object = c(paste("Cape Alitak/Humpy Deadman 257-10,20,50,60,70", sep = ''),
                                 paste("Ayakulik/Halibut Bay 256-10", "\u2013", "256-30", sep = ''),
                                 paste("Karluk/Sturgeon 255-10", "\u2013", "255-20; 256-40", sep = ''),
                                 paste("Uganik/Kupreanof 253", sep = ''),
                                 paste("Uyak Bay 254", sep = '')),
                      nm = unlist(strsplit(x = KMA2014, split = "14")))

ProportionColors <- colorpanel(n = 3, low = "blue", high = "white")
TempProportionColors14 <- sapply(KMA2014, function(geo) {
  ProportionColors[sapply(TempMix14[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
} )

Estimates <- KMA2014Strata_EstimatesStats
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5

# dir.create("Figures/2014")
require(devEMF)

sapply(names(TempMix14[c(4, 5, 3, 2, 1)]), function(geomix) {
  emf(file = paste("Figures/2014/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))

  Barplot1 <- barplot2(height = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})) * 100, 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})) * 100, 
                       ci.u = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})) * 100, 
                       ylim = c(0, 100), col = TempProportionColors14[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14[[geomix]], x = "topleft", fill = TempProportionColors14[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2014")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Percentage of Catch", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "14"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock specific harvest results 2014 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Multiply by harvest
HarvestByStrata2014
# Collapse 3_Late and 4_LateLate into 3_Late Final
HarvestByStrata2014_Final <- cbind(HarvestByStrata2014[, 1:2], "3_Late" = rowSums(HarvestByStrata2014[, 3:4], na.rm = TRUE))
dput(x = HarvestByStrata2014_Final, file = "Objects/HarvestByStrata2014_Final.txt")


KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")
str(KMA2014Strata_EstimatesStats)
names(KMA2014Strata_EstimatesStats)
dimnames(KMA2014Strata_EstimatesStats[[1]])


KMA2014Strata_HarvestEstimatesStats <- sapply(names(KMA2014Strata_EstimatesStats), function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(KMA2014Strata_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2014_Final[strata.split[1], strata.split[2]],
        KMA2014Strata_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2014Strata_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
KMA2014Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
str(KMA2014Strata_HarvestEstimatesStats)

# What should ymax be?
max(sapply(KMA2014Strata_HarvestEstimatesStats, function(strata) strata[, "95%"]))


TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_EstimatesStats), value = TRUE)} )

Legend14 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend14 <- sapply(KMA2014, function(geo) {
  Legend14[sapply(TempMix14[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
} )

GeoHeader <- setNames(object = c(paste("Cape Alitak/Humpy Deadman 257-10,20,50,60,70", sep = ''),
                                 paste("Ayakulik/Halibut Bay 256-10", "\u2013", "256-30", sep = ''),
                                 paste("Karluk/Sturgeon 255-10", "\u2013", "255-20; 256-40", sep = ''),
                                 paste("Uganik/Kupreanof 253", sep = ''),
                                 paste("Uyak Bay 254", sep = '')),
                      nm = unlist(strsplit(x = KMA2014, split = "14")))

HarvestColors <- colorpanel(n = 3, low = "green", high = "white")
TempHarvestColors14 <- sapply(KMA2014, function(geo) {
  HarvestColors[sapply(TempMix14[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
} )

Estimates <- KMA2014Strata_HarvestEstimatesStats
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5
ymax <- 140000

sapply(names(TempMix14[c(4, 5, 3, 2, 1)]), function(geomix) {
  emf(file = paste("Figures/2014/", geomix, "Harvest.emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})), 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})), 
                       ci.u = t(sapply(TempMix14[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})), 
                       ylim = c(0, ymax), col = TempHarvestColors14[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 20000), labels = formatC(x = seq(0, ymax, 20000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14[[geomix]], x = "topleft", fill = TempHarvestColors14[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2014")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Number of Fish Harvested (Thousands)", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "14"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Regional Roll-Ups 2014 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HarvestByStrata2014

paste(geomix, colnames(HarvestByStrata2014)[!is.na(HarvestByStrata2014[geomix,])], sep = "_")


HarvestByStrata2014[geomix, !is.na(HarvestByStrata2014[geomix,])]

# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
# Need to test first to make sure it will work with mixtures of different length (Uganik 3 had to be run out to 80K)
sapply(KMA2014, function(geomix) {
  assign(x = paste(geomix, "_Annual_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
           mixvec = paste(geomix, colnames(HarvestByStrata2014)[!is.na(HarvestByStrata2014[geomix,])], sep = "_"), 
           catchvec = HarvestByStrata2014[geomix, !is.na(HarvestByStrata2014[geomix,])], 
           newname = paste(geomix, "_Annual_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_Annual_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_Annual_Stratified.txt", sep = ''))
} ); beep(5)

str(SALITC14_Annual_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2014_Annual_EstimatesStats <- sapply(KMA2014, function(geomix) {
  Stats <- get(paste(geomix, "_Annual_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2014_Annual_EstimatesStats)
dput(x = KMA2014_Annual_EstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_EstimatesStats.txt")




# Create a list object with all Stratified Annual Harvest
KMA2014_Annual_HarvestEstimatesStats <- sapply(KMA2014, function(strata) {
  cbind(KMA2014_Annual_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2014_Final[strata, ], na.rm = TRUE),
        KMA2014_Annual_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2014_Annual_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_HarvestEstimatesStats.txt")
str(KMA2014_Annual_HarvestEstimatesStats)






# Create a matrix of annual means
Annual2014_Stratified_Estimates <- sapply(KMA2014, function(geomix) {
  KMA2014_Annual_EstimatesStats[[geomix]][, "mean"]
})

# Create a matrix of early strata means
KMA2014Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_1_Early_EstimatesStats.txt")
EarlyStrata2014_Estimates <- sapply(KMA2014Strata_1_Early_EstimatesStats, function(geomix) {
  geomix[, "mean"]
})
colnames(EarlyStrata2014_Estimates) <- KMA2014[-1]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Annual vs. Early Means
par(mar = c(4.1, 5.1, 3.1, 1.1))
sapply(KMA2014[-1], function(geomix) {
  Barplot <- barplot2(height = rbind(Annual2014_Stratified_Estimates[, geomix], EarlyStrata2014_Estimates[, geomix]) * 100, 
           beside = TRUE, col = c("blue", "white"), yaxt = 'n', xaxt = 'n', main = geomix, ylab = "Precent of Mixture",
           cex.lab = 2, cex.main = 2, ylim = c(0, 100))
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = ",", digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  text(x = colMeans(Barplot), y = -1, labels = KMA14GroupsPC2Rows, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)
  legend("topleft", legend = c("Annual Means", "Early Strata Means"), fill = c("blue", "white"), bty = 'n', cex = 1.5)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Temporal Mixture ANOVA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA2014_2015Strata

# Pool 3_Late and 4_LateLate
sapply(KMA2014[3:5], function(silly) {
  PoolCollections.GCL(collections = c(paste(silly, "_3_Late", sep = ''), paste(silly, "_4_LateLate", sep = '')),
                      loci = loci46, newname = paste(silly, "_3_LateAll", sep = ''))
})
str(SKARLC14_3_LateAll.gcl)

KMA2014_2015Strata_All <- c(sapply(KMA2014[1:2], function(silly) {paste(silly, c("_1_Early", "_2_Middle", "_3_Late"), sep = '')} )[-1],
                            sapply(KMA2014[3:5], function(silly) {paste(silly, c("_1_Early", "_2_Middle", "_3_LateAll"), sep = '')} ))

dir.create("VarComps")
VarComps <- VarComp.GCL(sillyvec = KMA2014_2015Strata_All, loci = loci46, 
                        groupvec = unlist(sapply(KMA2014, function(silly) {rep(which(KMA2014 == silly), length(grep(pattern = silly, x = KMA2014_2015Strata_All)))})),
                        dir = "VarComps")
str(VarComps)
# Did not finish this, may come back later, not all that important



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA files for BAYES 2015 Early Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(reshape)
samp.df.2015 <- data.frame(t(sapply(KMA2015Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2015, location ~ temporal.strata)

# Indetify Strata to Run
KMA2015Strata_1_Early <- grep(pattern = "1_Early", x = KMA2015Strata, value = TRUE)
Round1Mixtures_2015 <- KMA2015Strata_1_Early
dput(x = Round1Mixtures_2015, file = "Objects/Round1Mixtures_2015.txt")

# Create rolling prior based on 2014 Round 1 estimates
Round1Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")

Round1Mixtures_2015_Prior <- sapply(Round1Mixtures_2014_EstimatesStats[1:4], function(Mix) {
  Prior.GCL(groupvec = KMA473PopsGroupVec14, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round1Mixtures_2015_Prior) <- gsub(pattern = "C14", replacement = "C15", x = names(Round1Mixtures_2015_Prior))  # This changes the names
# Use a flat prior for SALITC15_1_Early as there was no 2014 Early strata for Alitak
Round1Mixtures_2015_Prior <- c(list("SALITC15_1_Early" = KMA473Pops15FlatPrior), Round1Mixtures_2015_Prior)
Round1Mixtures_2015 %in% names(Round1Mixtures_2015_Prior)

dput(x = Round1Mixtures_2015_Prior, file = "Objects/Round1Mixtures_2015_Prior.txt")
str(Round1Mixtures_2015_Prior)

# Verify
sapply(Round1Mixtures_2015, function(geomix) {plot(as.vector(Round1Mixtures_2015_Prior[[geomix]]), type = "h", main = geomix)})

## Dumping Mixture files
sapply(Round1Mixtures_2015, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round1Mixtures_2015, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round1Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round1Mixtures_2015, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})


## Dumping Control files for 80K suspects
sapply(Round1Mixtures_2015[4:5], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round1Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 2015 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC ,
  maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/40K Iterations", 
  mixvec = Round1Mixtures_2015, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round1Mixtures_2015_Estimates, file = "Estimates objects/Round1Mixtures_2015_Estimates.txt")
dput(Round1Mixtures_2015_Estimates$Stats, file = "Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")

Round1Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2015_Estimates.txt")
Round1Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round1Mixtures_2015_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round1Mixtures_2015_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUYAKC15_1_Early
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round1Mixtures_2015[c(4, 5, 3, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round1Mixtures_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round1Mixtures_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round1Mixtures_2015_Estimates$Output)
Round1Mixtures_2015_Header <- setNames(object = c("Alitak Section Early June 1-July 3, 2015",
                                                  "Ayakulik Section Early June 1-July 3, 2015",
                                                  "Karluk Section Early June 1-July 3, 2015",
                                                  "Uganik Section Early June 1-July 3, 2015",
                                                  "Uyak Section Early June 1-July 3, 2015"), 
                                       nm = Round1Mixtures_2015)
dput(x = Round1Mixtures_2015_Header, file = "Objects/Round1Mixtures_2015_Header.txt")

PlotPosterior(mixvec = Round1Mixtures_2015[c(4, 5, 3, 2, 1)], output = Round1Mixtures_2015_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round1Mixtures_2015_Header, set.mfrow = c(5, 3), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 2015 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round1Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round1Mixtures_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round1Mixtures_2015_Header)

## Quick barplot
QuickBarplot(mixvec = Round1Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round1Mixtures_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round1Mixtures_2015_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round1Mixtures_2015[c(4, 5, 3, 2, 1)], estimates = Round1Mixtures_2015_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round1Mixtures_2015_Header)
rm(Round1Mixtures_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 1 2015 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2015_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
                                                                   mixvec = Round1Mixtures_2015, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round1Mixtures_2015_31RG_Estimates)
QuickBarplot(mixvec = Round1Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round1Mixtures_2015_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round1Mixtures_2015_Header)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures80K_2015_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
                               maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/80K Iterations", 
                               mixvec = Round1Mixtures_2015[5], prior = "",  
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round1Mixtures80K_2015_Estimates, file = "Estimates objects/Round1Mixtures80K_2015_Estimates.txt")
dput(Round1Mixtures80K_2015_Estimates$Stats, file = "Estimates objects/Round1Mixtures80K_2015_EstimatesStats.txt")

Round1Mixtures80K_2015_Estimates <- dget(file = "Estimates objects/Round1Mixtures80K_2015_Estimates.txt")
Round1Mixtures80K_2015_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures80K_2015_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round1Mixtures80K_2015_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(Round1Mixtures80K_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUGANC15_1_Early still has PWS SEAK issues
require(gplots)
sapply(Round1Mixtures_2015[5], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round1Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round1Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round1Mixtures80K_2015_Estimates$Output)

PlotPosterior(mixvec = Round1Mixtures_2015[5], output = Round1Mixtures80K_2015_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round1Mixtures_2015_Header[4:5], set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round1Mixtures_2015[5], estimatesstats = Round1Mixtures80K_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round1Mixtures_2015_Header[5])

## Quick barplot
QuickBarplot(mixvec = Round1Mixtures_2015[5], estimatesstats = Round1Mixtures80K_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round1Mixtures_2015_Header[5])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round1Mixtures_2015[5], estimates = Round1Mixtures80K_2015_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round1Mixtures_2015_Header)
rm(Round1Mixtures80K_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round1Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")

sapply(Round1Mixtures_2015[5], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round1Mixtures_2015_EstimatesStats[[Mix]][, "median"],
                                                 Round1Mixtures80K_2015_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round1Mixtures_2015_EstimatesStats[[Mix]][, "5%"],
                                               Round1Mixtures80K_2015_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round1Mixtures_2015_EstimatesStats[[Mix]][, "95%"],
                                               Round1Mixtures80K_2015_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round1Mixtures_2015_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA14GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 1 2015 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SALITC15_1_Early, SAYAKC15_1_Early, SKARLC15_1_Early, and SUGANC15_1_Early as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from SUYAKC15_1_Early
Round1Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2015_Estimates.txt")

sapply(Round1Mixtures_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUYAKC15_1_Early
sapply(Round1Mixtures80K_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUYAKC15_1_Early resolved

str(Round1Mixtures_2015_Estimates)
str(Round1Mixtures80K_2015_Estimates)

Round1Mixtures_2015_Estimates_Final <- list(Stats = c(Round1Mixtures_2015_Estimates$Stats[Round1Mixtures_2015[1:4]],
                                                      Round1Mixtures80K_2015_Estimates$Stats[Round1Mixtures_2015[5]]),
                                            Output = c(Round1Mixtures_2015_Estimates$Output[Round1Mixtures_2015[1:4]],
                                                       Round1Mixtures80K_2015_Estimates$Output[Round1Mixtures_2015[5]]))
str(Round1Mixtures_2015_Estimates_Final)
dput(x = Round1Mixtures_2015_Estimates_Final, file = "Estimates objects/Round1Mixtures_2015_Estimates_Final.txt")
dput(x = Round1Mixtures_2015_Estimates_Final$Stats, file = "Estimates objects/Round1Mixtures_2015_EstimatesStats_Final.txt")

Round1Mixtures_2015_Estimates_Final <- dget(file = "Estimates objects/Round1Mixtures_2015_Estimates_Final.txt")
Round1Mixtures_2015_EstimatesStats_Final <- dget(file = "Estimates objects/Round1Mixtures_2015_EstimatesStats_Final.txt")

sapply(Round1Mixtures_2015_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUGANC15_1_Early still has issues, but we are going to live with it

# Dput final Early Strata Estimates from 2015
dput(x = Round1Mixtures_2015_EstimatesStats_Final, file = "Estimates objects/Final/KMA2015Strata_1_Early_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 1 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2015_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
                                                                   mixvec = Round1Mixtures_2015, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round1Mixtures_2015_31RG_Estimates)
QuickBarplot(mixvec = Round1Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round1Mixtures_2015_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round1Mixtures_2015_Header[c(4, 5, 3, 2, 1)])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA files for BAYES 2015 Middle Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(reshape)
samp.df.2015 <- data.frame(t(sapply(KMA2015Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2015, location ~ temporal.strata)

# Indetify Strata to Run
KMA2015Strata_2_Middle <- grep(pattern = "2_Middle", x = KMA2015Strata, value = TRUE)
Round2Mixtures_2015 <- KMA2015Strata_2_Middle
dput(x = Round2Mixtures_2015, file = "Objects/Round2Mixtures_2015.txt")


# Create rolling prior based on Round 1 estimates
Round2Mixtures_2015_Prior <- sapply(Round1Mixtures_2015_EstimatesStats_Final, function(Mix) {
  Prior.GCL(groupvec = KMA473PopsGroupVec14, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round2Mixtures_2015_Prior) <- gsub(pattern = "1_Early", replacement = "2_Middle", 
                                         x = names(Round2Mixtures_2015_Prior))  # This changes the names
dput(x = Round2Mixtures_2015_Prior, file = "Objects/Round2Mixtures_2015_Prior.txt")
str(Round2Mixtures_2015_Prior)


## Dumping Mixture files
sapply(Round2Mixtures_2015, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round2Mixtures_2015, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round2Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2015, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})

## Dumping Control files for 80K suspects
sapply(Round2Mixtures_2015[5], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round2Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})


## Dumping Control files for 80K suspects, repeat
sapply(Round2Mixtures_2015[4], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round2Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 2015 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC ,
  maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/40K Iterations", 
  mixvec = Round2Mixtures_2015, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures_2015_Estimates, file = "Estimates objects/Round2Mixtures_2015_Estimates.txt")
dput(Round2Mixtures_2015_Estimates$Stats, file = "Estimates objects/Round2Mixtures_2015_EstimatesStats.txt")

Round2Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2015_Estimates.txt")
Round2Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2015_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures_2015_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures_2015_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC15_2_Middle
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round2Mixtures_2015[c(4, 5, 3, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round2Mixtures_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures_2015_Estimates$Output)
Round2Mixtures_2015_Header <- setNames(object = c("Alitak Section Middle July 4-August 1, 2015",
                                                  "Ayakulik Section Middle July 4-August 1, 2015",
                                                  "Karluk Section Middle July 4-August 1, 2015",
                                                  "Uganik Section Middle July 4-August 1, 2015",
                                                  "Uyak Section Middle July 4-August 1, 2015"), 
                                       nm = Round2Mixtures_2015)
dput(x = Round2Mixtures_2015_Header, file = "Objects/Round2Mixtures_2015_Header.txt")

PlotPosterior(mixvec = Round2Mixtures_2015[c(4, 5, 3, 2, 1)], output = Round2Mixtures_2015_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round2Mixtures_2015_Header, set.mfrow = c(5, 3), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 2015 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round2Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round2Mixtures_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round2Mixtures_2015_Header)

## Quick barplot
QuickBarplot(mixvec = Round2Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round2Mixtures_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round2Mixtures_2015_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round2Mixtures_2015[c(4, 5, 3, 2, 1)], estimates = Round2Mixtures_2015_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round2Mixtures_2015_Header)
rm(Round2Mixtures_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 2 2015 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2015_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
                                                                   mixvec = Round2Mixtures_2015, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round2Mixtures_2015_31RG_Estimates)
QuickBarplot(mixvec = Round2Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round2Mixtures_2015_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round2Mixtures_2015_Header)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures80K_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/80K Iterations", 
  mixvec = Round2Mixtures_2015[4], prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures80K_2015_Estimates, file = "Estimates objects/Round2Mixtures80K_2015_Estimates.txt")
dput(Round2Mixtures80K_2015_Estimates$Stats, file = "Estimates objects/Round2Mixtures80K_2015_EstimatesStats.txt")

Round2Mixtures80K_2015_Estimates <- dget(file = "Estimates objects/Round2Mixtures80K_2015_Estimates.txt")
Round2Mixtures80K_2015_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures80K_2015_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures80K_2015_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures80K_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUYAKC15_2_Middle resolved
require(gplots)
sapply(Round2Mixtures_2015[4], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round2Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures80K_2015_Estimates$Output)

PlotPosterior(mixvec = Round2Mixtures_2015[4], output = Round2Mixtures80K_2015_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round2Mixtures_2015_Header[4], set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round2Mixtures_2015[4], estimatesstats = Round2Mixtures80K_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round2Mixtures_2015_Header[3])

## Quick barplot
QuickBarplot(mixvec = Round2Mixtures_2015[4], estimatesstats = Round2Mixtures80K_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round2Mixtures_2015_Header[4])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round2Mixtures_2015[4], estimates = Round2Mixtures80K_2015_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round2Mixtures_2015_Header)
rm(Round2Mixtures80K_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round2Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2015_EstimatesStats.txt")

sapply(Round2Mixtures_2015[4], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round2Mixtures_2015_EstimatesStats[[Mix]][, "median"],
                                                 Round2Mixtures80K_2015_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round2Mixtures_2015_EstimatesStats[[Mix]][, "5%"],
                                               Round2Mixtures80K_2015_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round2Mixtures_2015_EstimatesStats[[Mix]][, "95%"],
                                               Round2Mixtures80K_2015_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round2Mixtures_2015_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA14GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 2 2015 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SALITC15_2_Middle, SAYAKC15_2_Middle, SKARLC15_2_Middle, and SUYAKC15_2_Middle as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from  and SUGANC15_2_Middle
Round2Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2015_Estimates.txt")

sapply(Round2Mixtures_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC15_2_Middle
sapply(Round2Mixtures80K_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUGANC15_2_Middle resolved

str(Round2Mixtures_2015_Estimates)
str(Round2Mixtures80K_2015_Estimates)

Round2Mixtures_2015_Estimates_Final <- list(Stats = c(Round2Mixtures_2015_Estimates$Stats[Round2Mixtures_2015[1:3]],
                                                      Round2Mixtures80K_2015_Estimates$Stats[Round2Mixtures_2015[4]],
                                                      Round2Mixtures_2015_Estimates$Stats[Round2Mixtures_2015[5]]),
                                            Output = c(Round2Mixtures_2015_Estimates$Output[Round2Mixtures_2015[1:3]],
                                                       Round2Mixtures80K_2015_Estimates$Output[Round2Mixtures_2015[4]],
                                                       Round2Mixtures_2015_Estimates$Output[Round2Mixtures_2015[5]]))
str(Round2Mixtures_2015_Estimates_Final)
dput(x = Round2Mixtures_2015_Estimates_Final, file = "Estimates objects/Round2Mixtures_2015_Estimates_Final.txt")
dput(x = Round2Mixtures_2015_Estimates_Final$Stats, file = "Estimates objects/Round2Mixtures_2015_EstimatesStats_Final.txt")

Round2Mixtures_2015_Estimates_Final <- dget(file = "Estimates objects/Round2Mixtures_2015_Estimates_Final.txt")
Round2Mixtures_2015_EstimatesStats_Final <- dget(file = "Estimates objects/Round2Mixtures_2015_EstimatesStats_Final.txt")

sapply(Round2Mixtures_2015_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})

# Dput final Middle Strata Estimates from 2015
dput(x = Round2Mixtures_2015_EstimatesStats_Final, file = "Estimates objects/Final/KMA2015Strata_2_Middle_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 2 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2015_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
                                                                   mixvec = Round2Mixtures_2015, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round2Mixtures_2015_31RG_Estimates)
QuickBarplot(mixvec = Round2Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round2Mixtures_2015_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round2Mixtures_2015_Header[c(4, 5, 3, 2, 1)])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 MSA files for BAYES 2015 Late Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(reshape)
samp.df.2015 <- data.frame(t(sapply(KMA2015Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2015, location ~ temporal.strata)

# Indetify Strata to Run
KMA2015Strata_3_Late <- grep(pattern = "3_Late", x = KMA2015Strata, value = TRUE)
Round3Mixtures_2015 <- KMA2015Strata_3_Late
dput(x = Round3Mixtures_2015, file = "Objects/Round3Mixtures_2015.txt")


# Create rolling prior based on Round 2 estimates
Round3Mixtures_2015_Prior <- sapply(Round2Mixtures_2015_EstimatesStats_Final, function(Mix) {
  Prior.GCL(groupvec = KMA473PopsGroupVec14, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round3Mixtures_2015_Prior) <- gsub(pattern = "2_Middle", replacement = "3_Late", 
                                         x = names(Round3Mixtures_2015_Prior))  # This changes the names
dput(x = Round3Mixtures_2015_Prior, file = "Objects/Round3Mixtures_2015_Prior.txt")
str(Round3Mixtures_2015_Prior)


## Dumping Mixture files
sapply(Round3Mixtures_2015, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round3Mixtures_2015, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round3Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round3Mixtures_2015, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})


## Dumping Control files for 80K suspects
sapply(Round3Mixtures_2015[c(1, 4, 5)], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round3Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory for 80K suspects
sapply(Round3Mixtures_2015[c(1, 4, 5)], function(Mix) {dir.create(paste("BAYES/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 3 2015 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC ,
  maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
  mixvec = Round3Mixtures_2015, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round3Mixtures_2015_Estimates, file = "Estimates objects/Round3Mixtures_2015_Estimates.txt")
dput(Round3Mixtures_2015_Estimates$Stats, file = "Estimates objects/Round3Mixtures_2015_EstimatesStats.txt")

Round3Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round3Mixtures_2015_Estimates.txt")
Round3Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures_2015_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round3Mixtures_2015_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round3Mixtures_2015_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC14_3_Late
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round3Mixtures_2015[c(4, 5, 3, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round3Mixtures_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round3Mixtures_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round3Mixtures_2015_Estimates$Output)
Round3Mixtures_2015_Header <- setNames(object = c("Alitak Section Late August 2-29, 2015",
                                                  "Ayakulik Section Late August 2-29, 2015",
                                                  "Karluk Section Late August 2-29, 2015",
                                                  "Uganik Section Late August 2-29, 2015",
                                                  "Uyak Section Late August 2-29, 2015"), 
                                       nm = Round3Mixtures_2015)
dput(x = Round3Mixtures_2015_Header, file = "Objects/Round3Mixtures_2015_Header.txt")

PlotPosterior(mixvec = Round3Mixtures_2015[c(4, 5, 3, 2, 1)], output = Round3Mixtures_2015_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round3Mixtures_2015_Header, set.mfrow = c(5, 3), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 3 2015 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round3Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round3Mixtures_2015_Header)

## Quick barplot
QuickBarplot(mixvec = Round3Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round3Mixtures_2015_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round3Mixtures_2015[c(4, 5, 3, 2, 1)], estimates = Round3Mixtures_2015_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round3Mixtures_2015_Header)
rm(Round3Mixtures_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 3 2015 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2015_31RG_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
  mixvec = Round3Mixtures_2015, prior = "",  
  ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round3Mixtures_2015_31RG_Estimates)
QuickBarplot(mixvec = Round3Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2015_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round3Mixtures_2015_Header)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 3 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures80K_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output/80K Iterations", 
  mixvec = Round3Mixtures_2015[4], prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round3Mixtures80K_2015_Estimates, file = "Estimates objects/Round3Mixtures80K_2015_Estimates.txt")
dput(Round3Mixtures80K_2015_Estimates$Stats, file = "Estimates objects/Round3Mixtures80K_2015_EstimatesStats.txt")

Round3Mixtures80K_2015_Estimates <- dget(file = "Estimates objects/Round3Mixtures80K_2015_Estimates.txt")
Round3Mixtures80K_2015_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures80K_2015_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round3Mixtures80K_2015_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(Round3Mixtures80K_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUYAKC15_3_Late resolved
require(gplots)
sapply(Round3Mixtures_2015[4], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round3Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round3Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round3Mixtures80K_2015_Estimates$Output)

PlotPosterior(mixvec = Round3Mixtures_2015[4], output = Round3Mixtures80K_2015_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round3Mixtures_2015_Header[4], set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 3 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round3Mixtures_2015[4], estimatesstats = Round3Mixtures80K_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round3Mixtures_2015_Header[3])

## Quick barplot
QuickBarplot(mixvec = Round3Mixtures_2015[4], estimatesstats = Round3Mixtures80K_2015_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round3Mixtures_2015_Header[4])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round3Mixtures_2015[4], estimates = Round3Mixtures80K_2015_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round3Mixtures_2015_Header)
rm(Round3Mixtures80K_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round3Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures_2015_EstimatesStats.txt")

sapply(Round3Mixtures_2015[4], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round3Mixtures_2015_EstimatesStats[[Mix]][, "median"],
                                                 Round3Mixtures80K_2015_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round3Mixtures_2015_EstimatesStats[[Mix]][, "5%"],
                                               Round3Mixtures80K_2015_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round3Mixtures_2015_EstimatesStats[[Mix]][, "95%"],
                                               Round3Mixtures80K_2015_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round3Mixtures_2015_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA14GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 3 2015 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SALITC15_3_Late, SAYAKC15_3_Late, SKARLC15_3_Late, and SUYAKC15_3_Late as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from SUGANC15_3_Late
Round3Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round3Mixtures_2015_Estimates.txt")

sapply(Round3Mixtures_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC15_3_Late
sapply(Round3Mixtures80K_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC15_3_Late remain...

str(Round3Mixtures_2015_Estimates)
str(Round3Mixtures80K_2015_Estimates)

Round3Mixtures_2015_Estimates_Final <- list(Stats = c(Round3Mixtures_2015_Estimates$Stats[Round3Mixtures_2015[1:3]],
                                                      Round3Mixtures80K_2015_Estimates$Stats[Round3Mixtures_2015[4]],
                                                      Round3Mixtures_2015_Estimates$Stats[Round3Mixtures_2015[5]]),
                                            Output = c(Round3Mixtures_2015_Estimates$Output[Round3Mixtures_2015[1:3]],
                                                       Round3Mixtures80K_2015_Estimates$Output[Round3Mixtures_2015[4]],
                                                       Round3Mixtures_2015_Estimates$Output[Round3Mixtures_2015[5]]))
str(Round3Mixtures_2015_Estimates_Final)
dput(x = Round3Mixtures_2015_Estimates_Final, file = "Estimates objects/Round3Mixtures_2015_Estimates_Final.txt")
dput(x = Round3Mixtures_2015_Estimates_Final$Stats, file = "Estimates objects/Round3Mixtures_2015_EstimatesStats_Final.txt")

Round3Mixtures_2015_Estimates_Final <- dget(file = "Estimates objects/Round3Mixtures_2015_Estimates_Final.txt")
Round3Mixtures_2015_EstimatesStats_Final <- dget(file = "Estimates objects/Round3Mixtures_2015_EstimatesStats_Final.txt")

sapply(Round3Mixtures_2015_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC15_3_Late remain...

# Dput final Middle Strata Estimates from 2015
dput(x = Round3Mixtures_2015_EstimatesStats_Final, file = "Estimates objects/Final/KMA2015Strata_3_Late_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 3 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2015_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci 14RG/Output", 
                                                                   mixvec = Round3Mixtures_2015, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round3Mixtures_2015_31RG_Estimates)
QuickBarplot(mixvec = Round3Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2015_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round3Mixtures_2015_Header[c(4, 5, 3, 2, 1)])




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Final Estimates Objects for Each Temporal Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# View as tables by year
require(reshape)
samp.df.2015 <- data.frame(t(sapply(KMA2015Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2015, location ~ temporal.strata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock composition results 2015 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA2015Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_1_Early_EstimatesStats.txt")
KMA2015Strata_2_Middle_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_2_Middle_EstimatesStats.txt")

#~~~~~~~~
# Add Igvak 2_Middle 2015
str(KMA2015Strata_2_Middle_EstimatesStats)
SIGVAC15_2_Middle_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2016_EstimatesStats.txt")[1]
str(SIGVAC15_2_Middle_EstimatesStats)

KMA2015Strata_2_Middle_EstimatesStats <- c(KMA2015Strata_2_Middle_EstimatesStats[1:2], SIGVAC15_2_Middle_EstimatesStats, KMA2015Strata_2_Middle_EstimatesStats[3:5])
str(KMA2015Strata_2_Middle_EstimatesStats)
dput(x = KMA2015Strata_2_Middle_EstimatesStats, file = "Estimates objects/Final/KMA2015Strata_2_Middle_EstimatesStats.txt")
#~~~~~~~~

KMA2015Strata_3_Late_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_3_Late_EstimatesStats.txt")

KMA2015Strata_EstimatesStats <- c(KMA2015Strata_1_Early_EstimatesStats, 
                                  KMA2015Strata_2_Middle_EstimatesStats, 
                                  KMA2015Strata_3_Late_EstimatesStats)
dput(x = KMA2015Strata_EstimatesStats, file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")
KMA2015Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")

str(KMA2015Strata_EstimatesStats)


TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend15 <- setNames(object = c("June 1-July 3", "July 4-August 1", "August 2-29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend15 <- sapply(KMA2015, function(geo) {
  Legend15[sapply(TempMix15[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)

GeoHeader <- setNames(object = c(paste("Cape Alitak/Humpy Deadman 257-10,20,50,60,70", sep = ''),
                                 paste("Ayakulik/Halibut Bay 256-10", "\u2013", "256-30", sep = ''),
                                 paste("Karluk/Sturgeon 255-10", "\u2013", "255-20; 256-40", sep = ''),
                                 paste("Uganik/Kupreanof 253", sep = ''),
                                 paste("Uyak Bay 254", sep = '')),
                      nm = unlist(strsplit(x = KMA2015, split = "15")))
dput(x = GeoHeader, file = "Objects/GeoHeader.txt")

ProportionColors <- colorpanel(n = 3, low = "blue", high = "white")
TempProportionColors15 <- sapply(KMA2015, function(geo) {
  ProportionColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates <- KMA2015Strata_EstimatesStats
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5

# dir.create("Figures/2015")

sapply(names(TempMix15[c(4, 5, 3, 2, 1)]), function(geomix) {
  emf(file = paste("Figures/2015/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})) * 100, 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})) * 100, 
                       ci.u = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})) * 100, 
                       ylim = c(0, 100), col = TempProportionColors15[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix]], x = "topleft", fill = TempProportionColors15[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Percentage of Catch", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "15"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock specific harvest results 2015 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Multiply by harvest
HarvestByStrata2015
# Nothing to collapse, final is same as preliminary
HarvestByStrata2015_Final <- HarvestByStrata2015
dput(x = HarvestByStrata2015_Final, file = "Objects/HarvestByStrata2015_Final.txt")


KMA2015Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")
str(KMA2015Strata_EstimatesStats)
names(KMA2015Strata_EstimatesStats)
dimnames(KMA2015Strata_EstimatesStats[[1]])


KMA2015Strata_HarvestEstimatesStats <- sapply(names(KMA2015Strata_EstimatesStats), function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(KMA2015Strata_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2015_Final[strata.split[1], strata.split[2]],
        KMA2015Strata_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2015Strata_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015Strata_HarvestEstimatesStats.txt")
KMA2015Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_HarvestEstimatesStats.txt")
str(KMA2015Strata_HarvestEstimatesStats)


TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend15 <- setNames(object = c("June 1-July 3", "July 4-August 1", "August 2-29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend15 <- sapply(KMA2015, function(geo) {
  Legend15[sapply(TempMix15[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)

GeoHeader <- setNames(object = c(paste("Cape Alitak/Humpy Deadman 257-10,20,50,60,70", sep = ''),
                                 paste("Ayakulik/Halibut Bay 256-10", "\u2013", "256-30", sep = ''),
                                 paste("Karluk/Sturgeon 255-10", "\u2013", "255-20; 256-40", sep = ''),
                                 paste("Uganik/Kupreanof 253", sep = ''),
                                 paste("Uyak Bay 254", sep = '')),
                      nm = unlist(strsplit(x = KMA2015, split = "15")))

HarvestColors <- colorpanel(n = 3, low = "green", high = "white")
TempHarvestColors15 <- sapply(KMA2015, function(geo) {
  HarvestColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)



# What should ymax be?
max(sapply(KMA2015Strata_HarvestEstimatesStats, function(strata) strata[, "95%"]))


Estimates <- KMA2015Strata_HarvestEstimatesStats
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5
ymax <- 200000

sapply(names(TempMix15[c(4, 5, 3, 2, 1)]), function(geomix) {
  emf(file = paste("Figures/2015/", geomix, "Harvest.emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})), 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})), 
                       ci.u = t(sapply(TempMix15[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})), 
                       ylim = c(0, ymax), col = TempHarvestColors15[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 20000), labels = formatC(x = seq(0, ymax, 20000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix]], x = "topleft", fill = TempHarvestColors15[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Number of Fish Harvested (Thousands)", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "15"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Regional Roll-Ups 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HarvestByStrata2015

paste(geomix, colnames(HarvestByStrata2015)[!is.na(HarvestByStrata2015[geomix,])], sep = "_")


HarvestByStrata2015[geomix, !is.na(HarvestByStrata2015[geomix,])]

# Which mixtures were 80K as opposed to 40K?
Round1Mixtures_2015_Estimates_Final <- dget(file = "Estimates objects/Round1Mixtures_2015_Estimates_Final.txt")
Round2Mixtures_2015_Estimates_Final <- dget(file = "Estimates objects/Round2Mixtures_2015_Estimates_Final.txt")
Round3Mixtures_2015_Estimates_Final <- dget(file = "Estimates objects/Round3Mixtures_2015_Estimates_Final.txt")

sort(sapply(c(Round1Mixtures_2015_Estimates_Final$Output,
              Round2Mixtures_2015_Estimates_Final$Output,
              Round3Mixtures_2015_Estimates_Final$Output),
            function(mixture) {dim(mixture)[1] == 200000} ), decreasing = TRUE)


# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(KMA2015, function(geomix) {
  assign(x = paste(geomix, "_Annual_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = paste(geomix, colnames(HarvestByStrata2015)[!is.na(HarvestByStrata2015[geomix,])], sep = "_"), 
           catchvec = HarvestByStrata2015[geomix, !is.na(HarvestByStrata2015[geomix,])], 
           newname = paste(geomix, "_Annual_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_Annual_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_Annual_Stratified.txt", sep = ''))
} ); beep(5)

str(SALITC15_Annual_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2015_Annual_EstimatesStats <- sapply(KMA2015, function(geomix) {
  get(paste(geomix, "_Annual_Stratified", sep = ''))$Stats
}, simplify = FALSE)
str(KMA2015_Annual_EstimatesStats)

# Add Igvak 2_Middle 2015
str(KMA2015Strata_2_Middle_EstimatesStats)
SIGVAC15_2_Middle_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2016_EstimatesStats.txt")[1]
str(SIGVAC15_2_Middle_EstimatesStats)

KMA2015_Annual_EstimatesStats <- c(KMA2015_Annual_EstimatesStats[1:2], list("SIGVAC15" = SIGVAC15_2_Middle_EstimatesStats[[1]]), KMA2015_Annual_EstimatesStats[3:5])
dput(x = KMA2015_Annual_EstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_EstimatesStats.txt")




# Create a list object with all Stratified Annual Harvest
KMA2015_Annual_HarvestEstimatesStats <- sapply(KMA2015, function(strata) {
  cbind(KMA2015_Annual_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2015_Final[strata, ], na.rm = TRUE),
        KMA2015_Annual_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

str(KMA2015_Annual_HarvestEstimatesStats)
dput(x = KMA2015_Annual_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_HarvestEstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a matrix of annual means
#~~~~~~~~~~~~~~~~~~
# 2014
KMA2014_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_EstimatesStats.txt")
Annual2014_Stratified_Estimates <- sapply(KMA2014, function(geomix) {
  KMA2014_Annual_EstimatesStats[[geomix]][, "mean"]
})


KMA2014_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_HarvestEstimatesStats.txt")
Annual2014_Stratified_HarvestEstimates <- sapply(KMA2014, function(geomix) {
  round(KMA2014_Annual_HarvestEstimatesStats[[geomix]][, "mean"])
})



#~~~~~~~~~~~~~~~~~~
# 2015
KMA2015_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_EstimatesStats.txt")
Annual2015_Stratified_Estimates <- sapply(KMA2015, function(geomix) {
  KMA2015_Annual_EstimatesStats[[geomix]][, "mean"]
})


KMA2015_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_HarvestEstimatesStats.txt")
Annual2015_Stratified_HarvestEstimates <- sapply(KMA2015, function(geomix) {
  round(KMA2015_Annual_HarvestEstimatesStats[[geomix]][, "mean"])
})




require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2015_Stratified_Estimates[14:1, c(4,5,3,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2015 Proportion", at = seq(0, 1, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2014_Stratified_Estimates[14:1, c(4,5,3,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2014 Proportion", at = seq(0, 1, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })


# Bubble chart
require(ggplot2)
require(reshape2)
Annual2015_Stratified_Estimates_df <- melt(Annual2015_Stratified_Estimates)
names(Annual2015_Stratified_Estimates_df) <- c("RG", "Fishery", "Proportion")
ggplot(data = Annual2015_Stratified_Estimates_df, aes(x = Fishery, y = RG, size = Proportion)) + geom_point()


# Create a matrix of annual median harvest
Annual2015_Stratified_HarvestEstimates <- sapply(KMA2015, function(geomix) {
  round(KMA2015_Annual_HarvestEstimatesStats[[geomix]][, "mean"])
})


require(lattice)
max(Annual2015_Stratified_HarvestEstimates)

new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2015_Stratified_HarvestEstimates[14:1, c(4,5,3,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2015 Harvest", at = seq(0, 250000, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2014_Stratified_HarvestEstimates[14:1, c(4,5,3,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2014 Harvest", at = seq(0, 250000, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

sort(rowSums(Annual2015_Stratified_HarvestEstimates))
sort(rowSums(Annual2014_Stratified_HarvestEstimates))



# Bubble chart example
apply(col2rgb(KMA14Colors), 2, function(col) {rgb(red = col[1], green = col[2], blue = col[3], maxColorValue = 255)} )

# 2015
Annual2015_Stratified_HarvestEstimates_df <- melt(Annual2015_Stratified_HarvestEstimates)
names(Annual2015_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2015_Stratified_HarvestEstimates_df$RG <- factor(Annual2015_Stratified_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC))
Annual2015_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2015_Stratified_HarvestEstimates_df$Fishery, levels = KMA2015[c(4,5,3,2,1)])
Annual2015_Stratified_HarvestEstimates_df$Color <- rep(rev(KMA14Colors), 5)
str(Annual2015_Stratified_HarvestEstimates_df)


ggplot(data = Annual2015_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + geom_point() + scale_size_area(max_size = 25) + scale_color_manual(values = rev(rep(KMA14Colors, 5)))


ggplot(data = Annual2015_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, 250000), breaks = seq(50000, 250000, 50000), range = c(0, 25)) + 
  scale_color_manual(values = rev(rep(KMA14Colors, 5))) +
  ggtitle("2015 Harvest")


# 2014
Annual2014_Stratified_HarvestEstimates_df <- melt(Annual2014_Stratified_HarvestEstimates)
names(Annual2014_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2014_Stratified_HarvestEstimates_df$RG <- factor(Annual2014_Stratified_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC))
Annual2014_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2014_Stratified_HarvestEstimates_df$Fishery, levels = KMA2014[c(4,5,3,2,1)])
Annual2014_Stratified_HarvestEstimates_df$Color <- rep(rev(KMA14Colors), 5)
str(Annual2014_Stratified_HarvestEstimates_df)

ggplot(data = Annual2014_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + geom_point() + scale_size_area(max_size = 25) + scale_color_manual(values = rev(rep(KMA14Colors, 5)))

ggplot(data = Annual2014_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, 250000), breaks = seq(50000, 250000, 50000), range = c(0, 25)) + 
  scale_color_manual(values = rev(rep(KMA14Colors, 5))) +
  ggtitle("2014 Harvest")


# Create a matrix of early strata means
KMA2015Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_1_Early_EstimatesStats.txt")
EarlyStrata2015_Estimates <- sapply(KMA2015Strata_1_Early_EstimatesStats, function(geomix) {
  geomix[, "mean"]
})
colnames(EarlyStrata2015_Estimates) <- KMA2015


# KMA2014Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_1_Early_EstimatesStats.txt")
# EarlyStrata2014_Estimates <- sapply(KMA2014Strata_1_Early_EstimatesStats, function(geomix) {
#   geomix[, "mean"]
# })
# colnames(EarlyStrata2014_Estimates) <- KMA2014

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Annual vs. Early Means
par(mar = c(4.1, 5.1, 3.1, 1.1))
sapply(KMA2015, function(geomix) {
  Barplot <- barplot2(height = rbind(Annual2015_Stratified_Estimates[, geomix], EarlyStrata2015_Estimates[, geomix]) * 100, 
                      beside = TRUE, col = c("blue", "white"), yaxt = 'n', xaxt = 'n', main = geomix, ylab = "Precent of Mixture",
                      cex.lab = 2, cex.main = 2, ylim = c(0, 100))
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = ",", digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  text(x = colMeans(Barplot), y = -1, labels = KMA14GroupsPC2Rows, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)
  legend("topleft", legend = c("Annual Means", "Early Strata Means"), fill = c("blue", "white"), bty = 'n', cex = 1.5)
})



# par(mar = c(4.1, 5.1, 3.1, 1.1))
# sapply(KMA2014, function(geomix) {
#   Barplot <- barplot2(height = rbind(Annual2014_Stratified_Estimates[, geomix], EarlyStrata2014_Estimates[, geomix]) * 100, 
#                       beside = TRUE, col = c("blue", "white"), yaxt = 'n', xaxt = 'n', main = geomix, ylab = "Precent of Mixture",
#                       cex.lab = 2, cex.main = 2, ylim = c(0, 100))
#   axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = ",", digits = 0, format = "f"), cex.axis = 1.5)
#   abline(h = 0, xpd = FALSE)
#   text(x = colMeans(Barplot), y = -1, labels = KMA14GroupsPC2Rows, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)
#   legend("topleft", legend = c("Annual Means", "Early Strata Means"), fill = c("blue", "white"), bty = 'n', cex = 1.5)
# })


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl48.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("EASSIP-LateAugust2014", "OriginalLocusControl96.txt", "OriginalLocusControl48.txt", "LocusControl98.txt", "LocusControl99.txt", "OLD 15RG")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get post-QC, stratified, combined loci, mixtures
invisible(sapply(KMA2014_2015Strata, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_Strata_PostQC_CombinedLoci/", silly, ".txt", sep = "")), pos = 1)} )); beep(4)
objects(pattern = "\\.gcl")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pull genotypes from LOKI 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get collection SILLYs
KMA2015 <- c("SALITC15", "SAYAKC15",  "SIGVAC15", "SKARLC15", "SUGANC15", "SUYAKC15")  # Good Igvak samples, but very low fishing fishing
dput(x = KMA2015, file = "Objects/KMA2015.txt")
KMA2016 <- c("SALITC16", "SAYAKC16", "SIGVAC16", "SKARLC16", "SUGANC16", "SUYAKC16")  # We'll see
dput(x = KMA2016, file = "Objects/KMA2016.txt")

#~~~~~~~~~~~~~~~~~~
## Pull all data for each silly code and create .gcl objects for each
LOKI2R.GCL(sillyvec = c("SIGVAC15", KMA2016), username = username, password = password)
rm(username, password)
objects(pattern = "\\.gcl")


## Fix SKARLC16 date from 6/26/16 to 6/28/16
str(SKARLC16.gcl$attributes$CAPTURE_DATE)
SKARLC16.gcl$attributes$CAPTURE_DATE <- as.POSIXct(gsub(pattern = "2016-06-26", replacement = "2016-06-28", x = SKARLC16.gcl$attributes$CAPTURE_DATE))


## Save unaltered .gcl's as back-up:
invisible(sapply(c("SIGVAC15", KMA2016), function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/" , silly, ".txt", sep = ''))} )); beep(8)

## Original sample sizes by SILLY
collection.size.original <- sapply(c("SIGVAC15", KMA2016), function(silly) get(paste(silly, ".gcl", sep = ""))$n)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

## Get objects
LocusControl <- dget(file = "Objects/OriginalLocusControl48.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("EASSIP-LateAugust2014", "OLD 15RG", "OriginalLocusControl96.txt", "OriginalLocusControl48.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get un-altered mixtures
invisible(sapply(c("SIGVAC15", KMA2016), function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Define strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## All fish have a capture date?
sapply(c("SIGVAC15", KMA2016), function(silly) {sum(is.na(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE))} )  # Zeros are good

## Confirming samples sizes by date
sapply(c("SIGVAC15", KMA2016), function(silly) {table(get(paste(silly, ".gcl", sep = ''))$attributes$CAPTURE_DATE)} )


## Get dataframes of strata dates
KMA.Strata.Dates.2015Igvak <- read.xlsx(file = "MixtureStrataDates.xlsx", sheetName = "2015Igvak", header = TRUE)
KMA.Strata.Dates.2016Igvak <- read.xlsx(file = "MixtureStrataDates.xlsx", sheetName = "2016Igvak", header = TRUE)
KMA.Strata.Dates.2016 <- read.xlsx(file = "MixtureStrataDates.xlsx", sheetName = "2016", header = TRUE)


## Function to define strata by dates (date.df)

# # Inputs
# silly <- "SKARLC14"
# date.df <- KMA.Strata.Dates.2014KarlUganUyak
# loci <- loci48

PoolCollectionsByDateDF <- function(silly, date.df, loci) {
  sapply(silly, function(mix) {
    mix.dates <- unique(as.Date(get(paste(mix, ".gcl", sep = ''))$attributes$CAPTURE_DATE))
    by(data = date.df, INDICES = date.df$Strata, function(x) {
      IDs <- AttributesToIDs.GCL(silly = mix, attribute = "CAPTURE_DATE", matching = mix.dates[mix.dates >= x$Begin & mix.dates <= x$End])
      IDs <- list(na.omit(IDs))
      names(IDs) <- mix
      PoolCollections.GCL(collections = mix, loci = loci, IDs = IDs, newname = paste(mix, as.character(x$Strata), sep = "_"))
      list("First Last Fish" = range(as.numeric(unlist(IDs))), "n" = get(paste(mix, "_", as.character(x$Strata), ".gcl", sep = ''))$n)
    } )
  }, simplify = FALSE, USE.NAMES = TRUE)
}

# # Example
# PoolCollectionsByDateDF(silly = LateAugustMixtures2014[c(1, 3)], date.df = KMA.Strata.Dates.2014KarlUganUyak, loci = loci48)

PoolCollectionsByDateDF(silly = "SIGVAC15", date.df = KMA.Strata.Dates.2015Igvak, loci = loci48)
PoolCollectionsByDateDF(silly = "SIGVAC16", date.df = KMA.Strata.Dates.2016Igvak, loci = loci48)
PoolCollectionsByDateDF(silly = KMA2016[-3], date.df = KMA.Strata.Dates.2016, loci = loci48)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Mixture sillyvec
KMA2015Strata <- unlist(strsplit(x = grep(pattern = "15_", x = objects(pattern = "\\.gcl"), value = TRUE), split = "\\.gcl"))
KMA2016Strata <- unlist(strsplit(x = grep(pattern = "16_", x = objects(pattern = "\\.gcl"), value = TRUE), split = "\\.gcl"))

KMA2014_2016Strata <- c(KMA2014Strata, KMA2015Strata, KMA2016Strata)


dput(x = KMA2015Strata, file = "Objects/KMA2015Strata.txt")
dput(x = KMA2016Strata, file = "Objects/KMA2016Strata.txt")
dput(x = KMA2014_2016Strata, file = "Objects/KMA2014_2016Strata.txt")

# Confirm sample sizes
sapply(KMA2014_2016Strata, function(silly) get(paste(silly, ".gcl", sep = ""))$n)

# View as tables by year
require(reshape)
samp.df.2014 <- data.frame(t(sapply(KMA2014Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2014, location ~ temporal.strata)

samp.df.2015 <- data.frame(t(sapply(KMA2015Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2015, location ~ temporal.strata)

samp.df.2016 <- data.frame(t(sapply(KMA2016Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2016, location ~ temporal.strata)


# dput mixture sillys
invisible(sapply(c("SIGVAC15_2_Middle", KMA2016Strata), function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_Strata/" , silly, ".txt", sep = ''))} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

c("SIGVAC15_2_Middle", KMA2016Strata)

KMA2016Strata_SampleSizes <- matrix(data = NA, nrow = length(c("SIGVAC15_2_Middle", KMA2016Strata)), ncol = 5, 
                                         dimnames = list(c("SIGVAC15_2_Middle", KMA2016Strata), c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_KMA2016Strata_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = c("SIGVAC15_2_Middle", KMA2016Strata), loci = loci48)
min(Original_KMA2016Strata_SampleSizebyLocus)  ## 210/223
apply(Original_KMA2016Strata_SampleSizebyLocus, 1, min) / apply(Original_KMA2016Strata_SampleSizebyLocus, 1, max)  ## Good, 0.913

Original_KMA2016Strata_PercentbyLocus <- apply(Original_KMA2016Strata_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_KMA2016Strata_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_KMA2016Strata_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


#### Check individuals
## View Histogram of Failure Rate by Strata
invisible(sapply(c("SIGVAC15_2_Middle", KMA2016Strata), function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  failure <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(ind == "0") / length(ind)} )
  hist(x = failure, main = mix, xlab = "Failure Rate", col = 8, xlim = c(0, 1), ylim = c(0, 20), breaks = seq(from = 0, to = 1, by = 0.02))
  abline(v = 0.2, lwd = 3)
}))



### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_KMA2016Strata_ColSize <- sapply(paste(c("SIGVAC15_2_Middle", KMA2016Strata), ".gcl", sep = ''), function(x) get(x)$n)
KMA2016Strata_SampleSizes[, "Genotyped"] <- Original_KMA2016Strata_ColSize


### Alternate
## Indentify alternate species individuals
KMA2016Strata_Alternate <- FindAlternateSpecies.GCL(sillyvec = c("SIGVAC15_2_Middle", KMA2016Strata), species = "sockeye")

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport = KMA2016Strata_Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)

## Get number of individuals per silly after removing alternate species individuals
ColSize_KMA2016Strata_PostAlternate <- sapply(paste(c("SIGVAC15_2_Middle", KMA2016Strata), ".gcl", sep = ''), function(x) get(x)$n)
KMA2016Strata_SampleSizes[, "Alternate"] <- Original_KMA2016Strata_ColSize-ColSize_KMA2016Strata_PostAlternate


### Missing
## Remove individuals with >20% missing data
KMA2016Strata_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = c("SIGVAC15_2_Middle", KMA2016Strata), proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_KMA2016Strata_PostMissLoci <- sapply(paste(c("SIGVAC15_2_Middle", KMA2016Strata), ".gcl", sep = ''), function(x) get(x)$n)
KMA2016Strata_SampleSizes[, "Missing"] <- ColSize_KMA2016Strata_PostAlternate-ColSize_KMA2016Strata_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
KMA2016Strata_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = c("SIGVAC15_2_Middle", KMA2016Strata), loci = loci48, quantile = NULL, minproportion = 0.95)
KMA2016Strata_DuplicateCheckReportSummary <- sapply(c("SIGVAC15_2_Middle", KMA2016Strata), function(x) KMA2016Strata_DuplicateCheck95MinProportion[[x]]$report)
KMA2016Strata_DuplicateCheckReportSummary

## Remove duplicate individuals
KMA2016Strata_RemovedDups <- RemoveDups.GCL(KMA2016Strata_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_KMA2016Strata_PostDuplicate <- sapply(paste(c("SIGVAC15_2_Middle", KMA2016Strata), ".gcl", sep = ''), function(x) get(x)$n)
KMA2016Strata_SampleSizes[, "Duplicate"] <- ColSize_KMA2016Strata_PostMissLoci-ColSize_KMA2016Strata_PostDuplicate


### Final
KMA2016Strata_SampleSizes[, "Final"] <- ColSize_KMA2016Strata_PostDuplicate
KMA2016Strata_SampleSizes

write.xlsx(KMA2016Strata_SampleSizes, file = "Output/KMA2016Strata_SampleSizes.xlsx")
dput(x = KMA2016Strata_SampleSizes, file = "Objects/KMA2016Strata_SampleSizes.txt")


KMA2014_2016Strata_SampleSizes <- rbind(KMA2014_2015Strata_SampleSizes, KMA2016Strata_SampleSizes)
write.xlsx(KMA2014_2016Strata_SampleSizes, file = "Output/KMA2014_2016Strata_SampleSizes.xlsx")
dput(x = KMA2014_2016Strata_SampleSizes, file = "Objects/KMA2014_2016Strata_SampleSizes.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Combine Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### loci46
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get the final locus set from baseline file to avoid any errors
loci46

## Combine loci
combined.loci.46 <- sapply(grep(pattern = "\\.", x = loci46, value = TRUE), function(locus) {unlist(strsplit(x = locus, split = "\\."))}, simplify = FALSE)
sapply(combined.loci.46, function(loci2combine) {CombineLoci.GCL(sillyvec = c("SIGVAC15_2_Middle", KMA2016Strata), markerset = loci2combine, update = TRUE)} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save PostQC/Combined loci .gcl's as back-up:
invisible(sapply(c("SIGVAC15_2_Middle", KMA2016Strata), function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections_Strata_PostQC_CombinedLoci/" , silly, ".txt", sep = ''))} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Geneop
# Kick out Genepop file to look for gross excesses of hets - using HWE probability test option with default settings to see P-values w/ Fis:
gcl2Genepop.GCL(sillyvec = c("SIGVAC15_2_Middle", KMA2016Strata), loci = loci48[-mito.loci48], path = "Genepop/KMA2016Strata_46nuclearloci.gen", VialNums = TRUE)

# Read in Genepop output .P file
HWE <- ReadGenepopHWE.GCL(file = "Genepop/KMA2016Strata_46nuclearloci.txt.P")

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
LocusControl <- dget(file = "Objects/LocusControl50.txt")

KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("EASSIP-LateAugust2014", "OriginalLocusControl96.txt", "OriginalLocusControl48.txt", "LocusControl98.txt", "LocusControl99.txt", "OLD 15RG")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get post-QC, stratified, combined loci, mixtures
invisible(sapply(KMA2014_2016Strata, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_Strata_PostQC_CombinedLoci/", silly, ".txt", sep = "")), pos = 1)} )); beep(4)
objects(pattern = "\\.gcl")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA files for BAYES 2016 Early Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(reshape)
samp.df.2016 <- data.frame(t(sapply(KMA2016Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2016, location ~ temporal.strata)

# Indetify Strata to Run
KMA2016Strata_1_Early <- grep(pattern = "1_Early", x = KMA2016Strata, value = TRUE)
Round1Mixtures_2016 <- c("SIGVAC15_2_Middle", KMA2016Strata_1_Early)
dput(x = Round1Mixtures_2016, file = "Objects/Round1Mixtures_2016.txt")

# Create rolling prior based on 2014 Round 1 estimates
Round1Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")

Round1Mixtures_2016_Prior <- sapply(Round1Mixtures_2015_EstimatesStats, function(Mix) {
  Prior.GCL(groupvec = KMA473PopsGroupVec14, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round1Mixtures_2016_Prior) <- gsub(pattern = "C15", replacement = "C16", x = names(Round1Mixtures_2016_Prior))  # This changes the names
# Use a flat prior for SIGVAC15_2_Middle and SIGVAC16_1_Early as these temporal strata have never yet been analyzed
Round1Mixtures_2016_Prior <- c(list("SIGVAC15_2_Middle" = KMA473Pops15FlatPrior), list("SIGVAC16_1_Early" = KMA473Pops15FlatPrior), Round1Mixtures_2016_Prior)
Round1Mixtures_2016 %in% names(Round1Mixtures_2016_Prior)

dput(x = Round1Mixtures_2016_Prior, file = "Objects/Round1Mixtures_2016_Prior.txt")
str(Round1Mixtures_2016_Prior)

# Verify
sapply(Round1Mixtures_2016, function(geomix) {plot(as.vector(Round1Mixtures_2016_Prior[[geomix]]), type = "h", main = geomix)})

## Dumping Mixture files
sapply(Round1Mixtures_2016, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2016 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round1Mixtures_2016, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round1Mixtures_2016_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2016 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round1Mixtures_2016, function(Mix) {dir.create(paste("BAYES/2014-2016 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})


## Dumping Control files for 80K suspects
sapply(Round1Mixtures_2016[5], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round1Mixtures_2016_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2016 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 2016 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2016_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC ,
  maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output/40K Iterations", 
  mixvec = Round1Mixtures_2016, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round1Mixtures_2016_Estimates, file = "Estimates objects/Round1Mixtures_2016_Estimates.txt")
dput(Round1Mixtures_2016_Estimates$Stats, file = "Estimates objects/Round1Mixtures_2016_EstimatesStats.txt")

Round1Mixtures_2016_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2016_Estimates.txt")
Round1Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2016_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round1Mixtures_2016_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round1Mixtures_2016_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SKARLC16_1_Early
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round1Mixtures_2016[c(4, 5, 3, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round1Mixtures_2016_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round1Mixtures_2016_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round1Mixtures_2016_Estimates$Output)
Round1Mixtures_2016_Header <- setNames(object = c("Cape Igvak Section Middle July 4-August 1, 2015",
                                                  "Alitak Section Early June 1-June 27, 2016",
                                                  "Ayakulik Section Early June 1-June 27, 2016",
                                                  "Cape Igvak Section Early June 1-June 27, 2016",
                                                  "Karluk Section Early June 1-June 27, 2016",
                                                  "Uganik Section Early June 1-June 27, 2016",
                                                  "Uyak Section Early June 1-June 27, 2016"), 
                                       nm = Round1Mixtures_2016)
dput(x = Round1Mixtures_2016_Header, file = "Objects/Round1Mixtures_2016_Header.txt")

PlotPosterior(mixvec = Round1Mixtures_2016[c(1, 4, 6, 7, 5, 3, 2)], output = Round1Mixtures_2016_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round1Mixtures_2016_Header, set.mfrow = c(5, 3), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 2016 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick barplot
QuickBarplot(mixvec = Round1Mixtures_2016[c(1, 4, 6, 7, 5, 3, 2)], estimatesstats = Round1Mixtures_2016_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round1Mixtures_2016_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round1Mixtures_2016[c(1, 4, 6, 7, 5, 3, 2)], estimates = Round1Mixtures_2016_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round1Mixtures_2016_Header)
rm(Round1Mixtures_2016_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 1 2016 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2016_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output/40K Iterations", 
                                                                   mixvec = Round1Mixtures_2016, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round1Mixtures_2016_31RG_Estimates)
QuickBarplot(mixvec = Round1Mixtures_2016[c(1, 4, 6, 7, 5, 3, 2)], estimatesstats = Round1Mixtures_2016_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round1Mixtures_2016_Header)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures80K_2016_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
                               maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output/80K Iterations", 
                               mixvec = Round1Mixtures_2016[5], prior = "",  
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round1Mixtures80K_2016_Estimates, file = "Estimates objects/Round1Mixtures80K_2016_Estimates.txt")
dput(Round1Mixtures80K_2016_Estimates$Stats, file = "Estimates objects/Round1Mixtures80K_2016_EstimatesStats.txt")

Round1Mixtures80K_2016_Estimates <- dget(file = "Estimates objects/Round1Mixtures80K_2016_Estimates.txt")
Round1Mixtures80K_2016_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures80K_2016_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round1Mixtures80K_2016_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(Round1Mixtures80K_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUGANC15_1_Early still has PWS SEAK issues
require(gplots)
sapply(Round1Mixtures_2016[5], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round1Mixtures80K_2016_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round1Mixtures80K_2016_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round1Mixtures80K_2016_Estimates$Output)

PlotPosterior(mixvec = Round1Mixtures_2016[5], output = Round1Mixtures80K_2016_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round1Mixtures_2016_Header[4:5], set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick barplot
QuickBarplot(mixvec = Round1Mixtures_2016[5], estimatesstats = Round1Mixtures80K_2016_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round1Mixtures_2016_Header[5])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round1Mixtures_2016[5], estimates = Round1Mixtures80K_2016_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round1Mixtures_2016_Header)
rm(Round1Mixtures80K_2016_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round1Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2016_EstimatesStats.txt")

sapply(Round1Mixtures_2016[5], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round1Mixtures_2016_EstimatesStats[[Mix]][, "median"],
                                                 Round1Mixtures80K_2016_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round1Mixtures_2016_EstimatesStats[[Mix]][, "5%"],
                                               Round1Mixtures80K_2016_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round1Mixtures_2016_EstimatesStats[[Mix]][, "95%"],
                                               Round1Mixtures80K_2016_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round1Mixtures_2016_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA14GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 1 2016 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SIGVAC15_2_Middle, SALITC16_1_Early, SAYAKC16_1_Early, SIGVAC_1_Early, SUGANC16_1_Early, and SUYAKC16_1_Early as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from SKARLC16_1_Early
Round1Mixtures_2016_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2016_Estimates.txt")

sapply(Round1Mixtures_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SKARLC16_1_Early
sapply(Round1Mixtures80K_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SKARLC16_1_Early still bad, but better

str(Round1Mixtures_2016_Estimates)
str(Round1Mixtures80K_2016_Estimates)

Round1Mixtures_2016_Estimates_Final <- list(Stats = c(Round1Mixtures_2016_Estimates$Stats[Round1Mixtures_2016[1:4]],
                                                      Round1Mixtures80K_2016_Estimates$Stats[Round1Mixtures_2016[5]],
                                                      Round1Mixtures_2016_Estimates$Stats[Round1Mixtures_2016[6:7]]),
                                            Output = c(Round1Mixtures_2016_Estimates$Output[Round1Mixtures_2016[1:4]],
                                                       Round1Mixtures80K_2016_Estimates$Output[Round1Mixtures_2016[5]],
                                                       Round1Mixtures_2016_Estimates$Output[Round1Mixtures_2016[6:7]]))
str(Round1Mixtures_2016_Estimates_Final)
dput(x = Round1Mixtures_2016_Estimates_Final, file = "Estimates objects/Round1Mixtures_2016_Estimates_Final.txt")
dput(x = Round1Mixtures_2016_Estimates_Final$Stats, file = "Estimates objects/Round1Mixtures_2016_EstimatesStats_Final.txt")

Round1Mixtures_2016_Estimates_Final <- dget(file = "Estimates objects/Round1Mixtures_2016_Estimates_Final.txt")
Round1Mixtures_2016_EstimatesStats_Final <- dget(file = "Estimates objects/Round1Mixtures_2016_EstimatesStats_Final.txt")

sapply(Round1Mixtures_2016_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SKARLC16_1_Early still has issues, but we are going to live with it

# Dput final Early Strata Estimates from 2016
dput(x = Round1Mixtures_2016_EstimatesStats_Final, file = "Estimates objects/Final/KMA2016Strata_1_Early_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA files for BAYES 2016 Middle Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(reshape)
samp.df.2016 <- data.frame(t(sapply(KMA2016Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2016, location ~ temporal.strata)

# Indetify Strata to Run
KMA2016Strata_2_Middle <- grep(pattern = "2_Middle", x = KMA2016Strata, value = TRUE)
Round2Mixtures_2016 <- KMA2016Strata_2_Middle
dput(x = Round2Mixtures_2016, file = "Objects/Round2Mixtures_2016.txt")


# Create rolling prior based on Round 1 estimates
Round2Mixtures_2016_Prior <- sapply(Round1Mixtures_2016_EstimatesStats_Final, function(Mix) {
  Prior.GCL(groupvec = KMA473PopsGroupVec14, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round2Mixtures_2016_Prior) <- gsub(pattern = "1_Early", replacement = "2_Middle", 
                                         x = names(Round2Mixtures_2016_Prior))  # This changes the names
Round2Mixtures_2016_Prior <- Round2Mixtures_2016_Prior[-1]
dput(x = Round2Mixtures_2016_Prior, file = "Objects/Round2Mixtures_2016_Prior.txt")
str(Round2Mixtures_2016_Prior)


## Dumping Mixture files
sapply(Round2Mixtures_2016, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2016 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round2Mixtures_2016, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round2Mixtures_2016_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2016 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2016, function(Mix) {dir.create(paste("BAYES/2014-2016 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})

## Dumping Control files for 80K
sapply(Round2Mixtures_2016[2], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round2Mixtures_2016_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2016 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 2016 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2016_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC ,
  maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output/40K Iterations", 
  mixvec = Round2Mixtures_2016, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures_2016_Estimates, file = "Estimates objects/Round2Mixtures_2016_Estimates.txt")
dput(Round2Mixtures_2016_Estimates$Stats, file = "Estimates objects/Round2Mixtures_2016_EstimatesStats.txt")

Round2Mixtures_2016_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2016_Estimates.txt")
Round2Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2016_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures_2016_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures_2016_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC15_2_Middle
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round2Mixtures_2016[c(4, 5, 3, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round2Mixtures_2016_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures_2016_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures_2016_Estimates$Output)
Round2Mixtures_2016_Header <- setNames(object = c("Alitak Section Middle June 28-July 25, 2016",
                                                  "Ayakulik Section Middle June 28-July 25, 2016",
                                                  "Igvak Section Middle June 28-July 25, 2016",
                                                  "Karluk Section Middle June 28-July 25, 2016",
                                                  "Uganik Section Middle June 28-July 25, 2016",
                                                  "Uyak Section Middle June 28-July 25, 2016"), 
                                       nm = Round2Mixtures_2016)
dput(x = Round2Mixtures_2016_Header, file = "Objects/Round2Mixtures_2016_Header.txt")

PlotPosterior(mixvec = Round2Mixtures_2016[c(3, 5, 6, 4, 2, 1)], output = Round2Mixtures_2016_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round2Mixtures_2016_Header, set.mfrow = c(5, 3), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 2016 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick barplot
QuickBarplot(mixvec = Round2Mixtures_2016[c(3, 5, 6, 4, 2, 1)], estimatesstats = Round2Mixtures_2016_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round2Mixtures_2016_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round2Mixtures_2016[c(3, 5, 6, 4, 2, 1)], estimates = Round2Mixtures_2016_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round2Mixtures_2016_Header)
rm(Round2Mixtures_2016_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 2 2016 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2016_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
                                                                   mixvec = Round2Mixtures_2016, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round2Mixtures_2016_31RG_Estimates)
QuickBarplot(mixvec = Round2Mixtures_2016[c(4, 5, 3, 2, 1)], estimatesstats = Round2Mixtures_2016_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round2Mixtures_2016_Header)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures80K_2016_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
                               maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output/80K Iterations", 
                               mixvec = Round2Mixtures_2016[2], prior = "",  
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures80K_2016_Estimates, file = "Estimates objects/Round2Mixtures80K_2016_Estimates.txt")
dput(Round2Mixtures80K_2016_Estimates$Stats, file = "Estimates objects/Round2Mixtures80K_2016_EstimatesStats.txt")

Round2Mixtures80K_2016_Estimates <- dget(file = "Estimates objects/Round2Mixtures80K_2016_Estimates.txt")
Round2Mixtures80K_2016_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures80K_2016_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures80K_2016_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures80K_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUGANC15_1_Early still has PWS SEAK issues
require(gplots)
sapply(Round2Mixtures_2016[2], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round2Mixtures80K_2016_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures80K_2016_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures80K_2016_Estimates$Output)

PlotPosterior(mixvec = Round2Mixtures_2016[2], output = Round2Mixtures80K_2016_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round2Mixtures_2016_Header[2], set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick barplot
QuickBarplot(mixvec = Round2Mixtures_2016[5], estimatesstats = Round2Mixtures80K_2016_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round2Mixtures_2016_Header[5])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round2Mixtures_2016[5], estimates = Round2Mixtures80K_2016_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round2Mixtures_2016_Header)
rm(Round2Mixtures80K_2016_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round2Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2016_EstimatesStats.txt")

sapply(Round2Mixtures_2016[2], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round2Mixtures_2016_EstimatesStats[[Mix]][, "median"],
                                                 Round2Mixtures80K_2016_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round2Mixtures_2016_EstimatesStats[[Mix]][, "5%"],
                                               Round2Mixtures80K_2016_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round2Mixtures_2016_EstimatesStats[[Mix]][, "95%"],
                                               Round2Mixtures80K_2016_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round2Mixtures_2016_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA14GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 2 2016 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SALITC16_2_Middle, SKARLC16_2_Middle, SIGVAC_2_Middle, SUGANC16_2_Middle, and SUYAKC16_2_Middle as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from SAYAKC16_2_Middle
Round2Mixtures_2016_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2016_Estimates.txt")

sapply(Round2Mixtures_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SAYAKC16_2_Middle
sapply(Round2Mixtures80K_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SAYAKC16_2_Middle resolved

str(Round2Mixtures_2016_Estimates)
str(Round2Mixtures80K_2016_Estimates)

Round2Mixtures_2016_Estimates_Final <- list(Stats = c(Round2Mixtures_2016_Estimates$Stats[Round2Mixtures_2016[1]],
                                                      Round2Mixtures80K_2016_Estimates$Stats[Round2Mixtures_2016[2]],
                                                      Round2Mixtures_2016_Estimates$Stats[Round2Mixtures_2016[3:6]]),
                                            Output = c(Round2Mixtures_2016_Estimates$Output[Round2Mixtures_2016[1]],
                                                       Round2Mixtures80K_2016_Estimates$Output[Round2Mixtures_2016[2]],
                                                       Round2Mixtures_2016_Estimates$Output[Round2Mixtures_2016[3:6]]))
str(Round2Mixtures_2016_Estimates_Final)
dput(x = Round2Mixtures_2016_Estimates_Final, file = "Estimates objects/Round2Mixtures_2016_Estimates_Final.txt")
dput(x = Round2Mixtures_2016_Estimates_Final$Stats, file = "Estimates objects/Round2Mixtures_2016_EstimatesStats_Final.txt")

Round2Mixtures_2016_Estimates_Final <- dget(file = "Estimates objects/Round2Mixtures_2016_Estimates_Final.txt")
Round2Mixtures_2016_EstimatesStats_Final <- dget(file = "Estimates objects/Round2Mixtures_2016_EstimatesStats_Final.txt")

sapply(Round2Mixtures_2016_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SAYAKC16_2_Middle resolved

# Dput final Early Strata Estimates from 2016
dput(x = Round2Mixtures_2016_EstimatesStats_Final, file = "Estimates objects/Final/KMA2016Strata_2_Middle_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 MSA files for BAYES 2016 Middle Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(reshape)
samp.df.2016 <- data.frame(t(sapply(KMA2016Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2016, location ~ temporal.strata)

# Indetify Strata to Run
KMA2016Strata_3_Late <- grep(pattern = "3_Late", x = KMA2016Strata, value = TRUE)
Round3Mixtures_2016 <- KMA2016Strata_3_Late
dput(x = Round3Mixtures_2016, file = "Objects/Round3Mixtures_2016.txt")


# Create rolling prior based on Round 2 estimates
Round3Mixtures_2016_Prior <- sapply(Round2Mixtures_2016_EstimatesStats_Final, function(Mix) {
  Prior.GCL(groupvec = KMA473PopsGroupVec14, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round3Mixtures_2016_Prior) <- gsub(pattern = "2_Middle", replacement = "3_Late", 
                                         x = names(Round3Mixtures_2016_Prior))  # This changes the names
Round3Mixtures_2016_Prior <- Round3Mixtures_2016_Prior[-3]
dput(x = Round3Mixtures_2016_Prior, file = "Objects/Round3Mixtures_2016_Prior.txt")
str(Round3Mixtures_2016_Prior)


## Dumping Mixture files
sapply(Round3Mixtures_2016, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2016 Mixtures 46loci 14RG/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round3Mixtures_2016, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round3Mixtures_2016_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2016 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round3Mixtures_2016, function(Mix) {dir.create(paste("BAYES/2014-2016 Mixtures 46loci 14RG/Output/", Mix, sep = ""))})

## Dumping Control files for 80K
sapply(Round3Mixtures_2016[5], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec14, priorvec = Round3Mixtures_2016_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2016 Mixtures 46loci 14RG/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 3 2016 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2016_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC ,
  maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output/40K Iterations", 
  mixvec = Round3Mixtures_2016, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round3Mixtures_2016_Estimates, file = "Estimates objects/Round3Mixtures_2016_Estimates.txt")
dput(Round3Mixtures_2016_Estimates$Stats, file = "Estimates objects/Round3Mixtures_2016_EstimatesStats.txt")

Round3Mixtures_2016_Estimates <- dget(file = "Estimates objects/Round3Mixtures_2016_Estimates.txt")
Round3Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures_2016_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round3Mixtures_2016_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round3Mixtures_2016_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUYAKC15_3_Late
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round3Mixtures_2016[c(4, 5, 3, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round3Mixtures_2016_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round3Mixtures_2016_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round3Mixtures_2016_Estimates$Output)
Round3Mixtures_2016_Header <- setNames(object = c("Alitak Section Middle July 26-August 29, 2016",
                                                  "Ayakulik Section Middle July 26-August 29, 2016",
                                                  "Karluk Section Middle July 26-August 29, 2016",
                                                  "Uganik Section Middle July 26-August 29, 2016",
                                                  "Uyak Section Middle July 26-August 29, 2016"), 
                                       nm = Round3Mixtures_2016)
dput(x = Round3Mixtures_2016_Header, file = "Objects/Round3Mixtures_2016_Header.txt")

PlotPosterior(mixvec = Round3Mixtures_2016[c(4, 5, 3, 2, 1)], output = Round3Mixtures_2016_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round3Mixtures_2016_Header, set.mfrow = c(5, 3), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 3 2016 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick barplot
QuickBarplot(mixvec = Round3Mixtures_2016[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2016_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round3Mixtures_2016_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round3Mixtures_2016[c(4, 5, 3, 2, 1)], estimates = Round3Mixtures_2016_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round3Mixtures_2016_Header)
rm(Round3Mixtures_2016_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 3 2016 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2016_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output/40K Iterations", 
                                                                   mixvec = Round3Mixtures_2016, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round3Mixtures_2016_31RG_Estimates)
QuickBarplot(mixvec = Round3Mixtures_2016[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2016_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round3Mixtures_2016_Header)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 3 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures80K_2016_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, 
                               maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output/80K Iterations", 
                               mixvec = Round3Mixtures_2016[5], prior = "",  
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round3Mixtures80K_2016_Estimates, file = "Estimates objects/Round3Mixtures80K_2016_Estimates.txt")
dput(Round3Mixtures80K_2016_Estimates$Stats, file = "Estimates objects/Round3Mixtures80K_2016_EstimatesStats.txt")

Round3Mixtures80K_2016_Estimates <- dget(file = "Estimates objects/Round3Mixtures80K_2016_Estimates.txt")
Round3Mixtures80K_2016_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures80K_2016_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round3Mixtures80K_2016_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(Round3Mixtures80K_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUGANC15_1_Early still has PWS SEAK issues
require(gplots)
sapply(Round3Mixtures_2016[2], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round3Mixtures80K_2016_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round3Mixtures80K_2016_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA14GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round3Mixtures80K_2016_Estimates$Output)

PlotPosterior(mixvec = Round3Mixtures_2016[5], output = Round3Mixtures80K_2016_Estimates$Output, 
              groups = KMA14GroupsPC, colors = KMA14Colors, 
              header = Round3Mixtures_2016_Header[2], set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 3 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick barplot
QuickBarplot(mixvec = Round3Mixtures_2016[5], estimatesstats = Round3Mixtures80K_2016_Estimates, groups = KMA14GroupsPC, groups2rows = KMA14GroupsPC2Rows, header = Round3Mixtures_2016_Header[5])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round3Mixtures_2016[5], estimates = Round3Mixtures80K_2016_Estimates, groups = KMA14GroupsPC2Rows, colors = KMA14Colors, header = Round3Mixtures_2016_Header)
rm(Round3Mixtures80K_2016_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round3Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures_2016_EstimatesStats.txt")

sapply(Round3Mixtures_2016[5], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round3Mixtures_2016_EstimatesStats[[Mix]][, "median"],
                                                 Round3Mixtures80K_2016_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round3Mixtures_2016_EstimatesStats[[Mix]][, "5%"],
                                               Round3Mixtures80K_2016_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round3Mixtures_2016_EstimatesStats[[Mix]][, "95%"],
                                               Round3Mixtures80K_2016_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA14GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round3Mixtures_2016_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA14GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 3 2016 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SALITC16_3_Late, SAYAKC16_3_Late, SKARLC16_3_Late, and SUGANC16_3_Late as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from SUYAKC16_3_Late
Round3Mixtures_2016_Estimates <- dget(file = "Estimates objects/Round3Mixtures_2016_Estimates.txt")

sapply(Round3Mixtures_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUYAKC16_3_Late
sapply(Round3Mixtures80K_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUYAKC16_3_Late resolved

str(Round3Mixtures_2016_Estimates)
str(Round3Mixtures80K_2016_Estimates)

Round3Mixtures_2016_Estimates_Final <- list(Stats = c(Round3Mixtures_2016_Estimates$Stats[Round3Mixtures_2016[1:4]],
                                                      Round3Mixtures80K_2016_Estimates$Stats[Round3Mixtures_2016[5]]),
                                            Output = c(Round3Mixtures_2016_Estimates$Output[Round3Mixtures_2016[1:4]],
                                                       Round3Mixtures80K_2016_Estimates$Output[Round3Mixtures_2016[5]]))
str(Round3Mixtures_2016_Estimates_Final)
dput(x = Round3Mixtures_2016_Estimates_Final, file = "Estimates objects/Round3Mixtures_2016_Estimates_Final.txt")
dput(x = Round3Mixtures_2016_Estimates_Final$Stats, file = "Estimates objects/Round3Mixtures_2016_EstimatesStats_Final.txt")

Round3Mixtures_2016_Estimates_Final <- dget(file = "Estimates objects/Round3Mixtures_2016_Estimates_Final.txt")
Round3Mixtures_2016_EstimatesStats_Final <- dget(file = "Estimates objects/Round3Mixtures_2016_EstimatesStats_Final.txt")

sapply(Round3Mixtures_2016_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUYAKC16_3_Late resolved

# Dput final Early Strata Estimates from 2016
dput(x = Round3Mixtures_2016_EstimatesStats_Final, file = "Estimates objects/Final/KMA2016Strata_3_Late_EstimatesStats.txt")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Final Estimates Objects for Each Temporal Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# View as tables by year
require(reshape)
samp.df.2016 <- data.frame(t(sapply(KMA2016Strata, function(strata) {
  location <- unlist(strsplit(x = strata, split = "_"))[1]
  temporal.strata <- as.numeric(unlist(strsplit(x = strata, split = "_"))[2])
  n <- get(paste(strata, ".gcl", sep = ""))$n
  c(location = location, temporal.strata = temporal.strata, n = n)
})))
cast(data = samp.df.2016, location ~ temporal.strata)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock composition results 2016 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA2016Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_1_Early_EstimatesStats.txt")
KMA2016Strata_2_Middle_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_2_Middle_EstimatesStats.txt")
KMA2016Strata_3_Late_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_3_Late_EstimatesStats.txt")

KMA2016Strata_EstimatesStats <- c(KMA2016Strata_1_Early_EstimatesStats, 
                                  KMA2016Strata_2_Middle_EstimatesStats, 
                                  KMA2016Strata_3_Late_EstimatesStats)
str(KMA2016Strata_EstimatesStats)

# Add Igvak 2015 to KMA2015Strata_EstimatesStats
KMA2015Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")
str(KMA2015Strata_EstimatesStats)
KMA2015Strata_EstimatesStats <- c(KMA2015Strata_EstimatesStats[1:7], KMA2016Strata_EstimatesStats[1], KMA2015Strata_EstimatesStats[8:15])
str(KMA2015Strata_EstimatesStats)
dput(x = KMA2015Strata_EstimatesStats, file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")


# Dput KMA2016Strata_EstimatesStats
KMA2016Strata_EstimatesStats <- KMA2016Strata_EstimatesStats[-1]
str(KMA2016Strata_EstimatesStats)
dput(x = KMA2016Strata_EstimatesStats, file = "Estimates objects/Final/KMA2016Strata_EstimatesStats.txt")
KMA2016Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_EstimatesStats.txt")

str(KMA2016Strata_EstimatesStats)




TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend16 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend16 <- sapply(KMA2016, function(geo) {
  Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)

GeoHeader <- setNames(object = c(paste("Cape Alitak/Humpy Deadman 257-10,20,50,60,70", sep = ''),
                                 paste("Ayakulik/Halibut Bay 256-10", "\u2013", "256-30", sep = ''),
                                 paste("Cape Igvak 262-75,80,90,95", sep = ''),
                                 paste("Karluk/Sturgeon 255-10", "\u2013", "255-20; 256-40", sep = ''),
                                 paste("Uganik/Kupreanof 253", sep = ''),
                                 paste("Uyak Bay 254", sep = '')),
                      nm = unlist(strsplit(x = KMA2016, split = "16")))
dput(x = GeoHeader, file = "Objects/GeoHeader.txt")

ProportionColors <- colorpanel(n = 3, low = "blue", high = "white")
TempProportionColors16 <- sapply(KMA2016, function(geo) {
  ProportionColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates <- KMA2016Strata_EstimatesStats
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5

# dir.create("Figures/2016")

sapply(names(TempMix16[c(3, 5, 6, 4, 2, 1)]), function(geomix) {
  emf(file = paste("Figures/2016/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})) * 100, 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})) * 100, 
                       ci.u = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})) * 100, 
                       ylim = c(0, 100), col = TempProportionColors16[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend16[[geomix]], x = "topleft", fill = TempProportionColors16[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Percentage of Catch", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "16"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot stock specific harvest results 2016 KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Multiply by harvest
HarvestByStrata2016
# Nothing to collapse, final is same as preliminary
HarvestByStrata2016_Final <- HarvestByStrata2016
dput(x = HarvestByStrata2016_Final, file = "Objects/HarvestByStrata2016_Final.txt")


KMA2016Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_EstimatesStats.txt")
str(KMA2016Strata_EstimatesStats)
names(KMA2016Strata_EstimatesStats)
dimnames(KMA2016Strata_EstimatesStats[[1]])


KMA2016Strata_HarvestEstimatesStats <- sapply(names(KMA2016Strata_EstimatesStats), function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(KMA2016Strata_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2016_Final[strata.split[1], strata.split[2]],
        KMA2016Strata_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2016Strata_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016Strata_HarvestEstimatesStats.txt")
KMA2016Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_HarvestEstimatesStats.txt")
str(KMA2016Strata_HarvestEstimatesStats)


TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend16 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend16 <- sapply(KMA2016, function(geo) {
  Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)

GeoHeader <- setNames(object = c(paste("Cape Alitak/Humpy Deadman 257-10,20,50,60,70", sep = ''),
                                 paste("Ayakulik/Halibut Bay 256-10", "\u2013", "256-30", sep = ''),
                                 paste("Cape Igvak 262-75,80,90,95", sep = ''),
                                 paste("Karluk/Sturgeon 255-10", "\u2013", "255-20; 256-40", sep = ''),
                                 paste("Uganik/Kupreanof 253", sep = ''),
                                 paste("Uyak Bay 254", sep = '')),
                      nm = unlist(strsplit(x = KMA2016, split = "16")))

HarvestColors <- colorpanel(n = 3, low = "green", high = "white")
TempHarvestColors16 <- sapply(KMA2016, function(geo) {
  HarvestColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)



# What should ymax be?
max(sapply(KMA2016Strata_HarvestEstimatesStats, function(strata) strata[, "95%"]))


Estimates <- KMA2016Strata_HarvestEstimatesStats
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5
ymax <- 200000

sapply(names(TempMix16[c(3, 5, 6, 4, 2, 1)]), function(geomix) {
  emf(file = paste("Figures/2016/", geomix, "Harvest.emf", sep = ''), width = 8.5, height = 6.5, family = "serif", bg = "white")
  par(mar = c(2.1, 4.1, 2.6, 0.6))
  
  Barplot1 <- barplot2(height = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "median"]})), 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "5%"]})), 
                       ci.u = t(sapply(TempMix16[[geomix]], function(tempmix) {Estimates[[tempmix]][, "95%"]})), 
                       ylim = c(0, ymax), col = TempHarvestColors16[[geomix]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 20000), labels = formatC(x = seq(0, ymax, 20000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend16[[geomix]], x = "topleft", fill = TempHarvestColors16[[geomix]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = "Number of Fish Harvested (Thousands)", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[unlist(strsplit(geomix, split = "16"))], side = 3, cex = cex.main, line = 1)
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Regional Roll-Ups 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HarvestByStrata2016

paste(geomix, colnames(HarvestByStrata2016)[!is.na(HarvestByStrata2016[geomix,])], sep = "_")


HarvestByStrata2016[geomix, !is.na(HarvestByStrata2016[geomix,])]

# Which mixtures were 80K as opposed to 40K?
Round1Mixtures_2016_Estimates_Final <- dget(file = "Estimates objects/Round1Mixtures_2016_Estimates_Final.txt")
Round2Mixtures_2016_Estimates_Final <- dget(file = "Estimates objects/Round2Mixtures_2016_Estimates_Final.txt")
Round3Mixtures_2016_Estimates_Final <- dget(file = "Estimates objects/Round3Mixtures_2016_Estimates_Final.txt")

sort(sapply(c(Round1Mixtures_2016_Estimates_Final$Output,
              Round2Mixtures_2016_Estimates_Final$Output,
              Round3Mixtures_2016_Estimates_Final$Output),
            function(mixture) {dim(mixture)[1] == 200000} ), decreasing = TRUE)


# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(KMA2016, function(geomix) {
  assign(x = paste(geomix, "_Annual_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = paste(geomix, colnames(HarvestByStrata2016)[!is.na(HarvestByStrata2016[geomix,])], sep = "_"), 
           catchvec = HarvestByStrata2016[geomix, !is.na(HarvestByStrata2016[geomix,])], 
           newname = paste(geomix, "_Annual_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_Annual_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_Annual_Stratified.txt", sep = ''))
} ); beep(5)

str(SALITC16_Annual_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2016_Annual_EstimatesStats <- sapply(KMA2016, function(geomix) {
  get(paste(geomix, "_Annual_Stratified", sep = ''))$Stats
}, simplify = FALSE)
str(KMA2016_Annual_EstimatesStats)
dput(x = KMA2016_Annual_EstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_EstimatesStats.txt")




# Create a list object with all Stratified Annual Harvest
KMA2016_Annual_HarvestEstimatesStats <- sapply(KMA2016, function(strata) {
  cbind(KMA2016_Annual_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2016_Final[strata, ], na.rm = TRUE),
        KMA2016_Annual_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

str(KMA2016_Annual_HarvestEstimatesStats)
dput(x = KMA2016_Annual_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_HarvestEstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Bubble Plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a matrix of annual medians

KMA14GroupsPC2 <- c(KMA14GroupsPC[1:4], "Ayakulik / Frazer", KMA14GroupsPC[6:14])
dput(x = KMA14GroupsPC2, file = "Objects/KMA14GroupsPC2.txt")

#~~~~~~~~~~~~~~~~~~
# 2014
KMA2014_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_EstimatesStats.txt")
Annual2014_Stratified_Estimates <- sapply(KMA2014, function(geomix) {
  KMA2014_Annual_EstimatesStats[[geomix]][, "mean"]
})
Annual2014_Stratified_Estimates <- cbind(Annual2014_Stratified_Estimates[, 1:2], "SIGVAC14" = rep(0, 14), Annual2014_Stratified_Estimates[, 3:5])


KMA2014_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_HarvestEstimatesStats.txt")
Annual2014_Stratified_HarvestEstimates <- sapply(KMA2014, function(geomix) {
  round(KMA2014_Annual_HarvestEstimatesStats[[geomix]][, "median"])
})
Annual2014_Stratified_HarvestEstimates <- cbind(Annual2014_Stratified_HarvestEstimates[, 1:2], "SIGVAC14" = rep(0, 14), Annual2014_Stratified_HarvestEstimates[, 3:5])
dimnames(Annual2014_Stratified_HarvestEstimates) <- list(KMA14GroupsPC2,
                                                         c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak"))





KMA2014_Temporal_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
Early2014_Stratified_HarvestEstimates <- sapply(grep(pattern = "Early", x = KMA2014Strata, value = TRUE), function(geomix) {
  round(KMA2014_Temporal_HarvestEstimatesStats[[geomix]][, "median"])
})
Early2014_Stratified_HarvestEstimates <- cbind("SALITC14_1_Early" = rep(0, 14), Early2014_Stratified_HarvestEstimates[, 1, drop = FALSE], "SIGVAC14" = rep(0, 14), Early2014_Stratified_HarvestEstimates[, 2:4])
dimnames(Early2014_Stratified_HarvestEstimates) <- list(KMA14GroupsPC2,
                                                         c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak"))




#~~~~~~~~~~~~~~~~~~
# 2015
KMA2015_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_EstimatesStats.txt")
Annual2015_Stratified_Estimates <- sapply(KMA2015, function(geomix) {
  KMA2015_Annual_EstimatesStats[[geomix]][, "mean"]
})


KMA2015_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_HarvestEstimatesStats.txt")
Annual2015_Stratified_HarvestEstimates <- sapply(KMA2015, function(geomix) {
  round(KMA2015_Annual_HarvestEstimatesStats[[geomix]][, "median"])
})
dimnames(Annual2015_Stratified_HarvestEstimates) <- list(KMA14GroupsPC2,
                                                         c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak"))



#~~~~~~~~~~~~~~~~~~
# 2016
KMA2016_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_EstimatesStats.txt")
Annual2016_Stratified_Estimates <- sapply(KMA2016, function(geomix) {
  KMA2016_Annual_EstimatesStats[[geomix]][, "mean"]
})


KMA2016_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_HarvestEstimatesStats.txt")
Annual2016_Stratified_HarvestEstimates <- sapply(KMA2016, function(geomix) {
  round(KMA2016_Annual_HarvestEstimatesStats[[geomix]][, "median"])
})
dimnames(Annual2016_Stratified_HarvestEstimates) <- list(KMA14GroupsPC2,
                                                         c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak"))

zmax <- max(Annual2014_Stratified_HarvestEstimates, 
            Annual2015_Stratified_HarvestEstimates,
            Annual2016_Stratified_HarvestEstimates)


require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2016_Stratified_Estimates[14:1, c(3,5,6,4,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2016 Proportion", at = seq(0, 1, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2015_Stratified_Estimates[14:1, c(3,5,6,4,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2015 Proportion", at = seq(0, 1, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2014_Stratified_Estimates[14:1, c(3,5,6,4,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2014 Proportion", at = seq(0, 1, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })


# Bubble chart
require(ggplot2)
require(reshape2)
require(devEMF)
Annual2015_Stratified_Estimates_df <- melt(Annual2015_Stratified_Estimates)
names(Annual2015_Stratified_Estimates_df) <- c("RG", "Fishery", "Proportion")
ggplot(data = Annual2015_Stratified_Estimates_df, aes(x = Fishery, y = RG, size = Proportion)) + geom_point()


# Create a matrix of annual median harvest
Annual2015_Stratified_HarvestEstimates <- sapply(KMA2015, function(geomix) {
  round(KMA2015_Annual_HarvestEstimatesStats[[geomix]][, "mean"])
})


require(lattice)

new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2016_Stratified_HarvestEstimates[14:1, c(3,5,6,4,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2016 Harvest", at = seq(0, 250000, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2015_Stratified_HarvestEstimates[14:1, c(3,5,6,4,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2015 Harvest", at = seq(0, 250000, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2014_Stratified_HarvestEstimates[14:1, c(3,5,6,4,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2014 Harvest", at = seq(0, 250000, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

sort(rowSums(Annual2016_Stratified_HarvestEstimates))
sort(rowSums(Annual2015_Stratified_HarvestEstimates))
sort(rowSums(Annual2014_Stratified_HarvestEstimates))


KMA_AnnualHarvest_14RG <- cbind("2014" = rowSums(Annual2014_Stratified_HarvestEstimates),
                                "2015" = rowSums(Annual2015_Stratified_HarvestEstimates),
                                '2016' = rowSums(Annual2016_Stratified_HarvestEstimates))
KMA_AnnualHarvest_14RG

apply(KMA_AnnualHarvest_14RG, 2, function(yr) {round(yr["Cook Inlet"] / sum(yr) * 100, 1)})

# Bubble chart example
apply(col2rgb(KMA14Colors), 2, function(col) {rgb(red = col[1], green = col[2], blue = col[3], maxColorValue = 255)} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2016
Annual2016_Stratified_HarvestEstimates_df <- melt(Annual2016_Stratified_HarvestEstimates)
names(Annual2016_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2016_Stratified_HarvestEstimates_df$RG <- factor(Annual2016_Stratified_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC2))
Annual2016_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2016_Stratified_HarvestEstimates_df$Fishery, levels = c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak"))
Annual2016_Stratified_HarvestEstimates_df$Color <- rep(rev(KMA14Colors), 5)
str(Annual2016_Stratified_HarvestEstimates_df)


ggplot(data = Annual2016_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + geom_point() + scale_size_area(max_size = 25) + scale_color_manual(values = rev(rep(KMA14Colors, 5)))


ggplot(data = Annual2016_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, zmax), breaks = seq(50000, 250000, 50000), range = c(0, 25)) + 
  scale_color_manual(values = rev(rep(KMA14Colors, 5))) +
  ggtitle("2016 Harvest")



# Figure for Report
Annual2016_Stratified_HarvestEstimates_df <- melt(Annual2016_Stratified_HarvestEstimates)
names(Annual2016_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2016_Stratified_HarvestEstimates_df$RG <- factor(Annual2016_Stratified_HarvestEstimates_df$RG, levels = KMA14GroupsPC2)
Annual2016_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2016_Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak")))
Annual2016_Stratified_HarvestEstimates_df$Color <- rep(KMA14Colors, 5)
Annual2016_Stratified_HarvestEstimates_df$Harvest[Annual2016_Stratified_HarvestEstimates_df$Harvest == 0] <- NA
str(Annual2016_Stratified_HarvestEstimates_df)

KMA14GroupsPC2RowsBubble <- KMA14GroupsPC2Rows
KMA14GroupsPC2RowsBubble[c(6,7,9,11,13)] <- gsub(pattern = "\n", replacement = "", x = KMA14GroupsPC2RowsBubble[c(6,7,9,11,13)])

emf(file ="Figures/All Years/2016 Harvest Bubble Plot.emf", width = 9, height = 5.75, family = "serif", bg = "white")

ggplot(data = Annual2016_Stratified_HarvestEstimates_df, aes(x = RG, y = Fishery, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(name = "Harvest\n(Thousands)", limits = c(0, zmax), breaks = c(1000, 5000, 10000, 20000, seq(50000, 250000, 50000)), range = c(0, 20), labels = c(1, 5, 10, 20, 50, 100, 150, 200, 250)) + 
  scale_color_manual(values = rep(KMA14Colors, 5), guide = FALSE) +
  scale_x_discrete(name = "Reporting Group", labels = KMA14GroupsPC2RowsBubble) +
  scale_y_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak"))) +
  theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = rel(1.3))) +
  theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,0.2,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(-0.5,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.7), angle = 00)) +
  theme(text = element_text(family = "times"))

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2015
Annual2015_Stratified_HarvestEstimates_df <- melt(Annual2015_Stratified_HarvestEstimates)
names(Annual2015_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2015_Stratified_HarvestEstimates_df$RG <- factor(Annual2015_Stratified_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC2))
Annual2015_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2015_Stratified_HarvestEstimates_df$Fishery, levels = c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak"))
Annual2015_Stratified_HarvestEstimates_df$Color <- rep(rev(KMA14Colors), 5)
str(Annual2015_Stratified_HarvestEstimates_df)


ggplot(data = Annual2015_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + geom_point() + scale_size_area(max_size = 25) + scale_color_manual(values = rev(rep(KMA14Colors, 5)))


ggplot(data = Annual2015_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, zmax), breaks = seq(50000, 250000, 50000), range = c(0, 25)) + 
  scale_color_manual(values = rev(rep(KMA14Colors, 5))) +
  ggtitle("2015 Harvest")




# Figure for Report
Annual2015_Stratified_HarvestEstimates_df <- melt(Annual2015_Stratified_HarvestEstimates)
names(Annual2015_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2015_Stratified_HarvestEstimates_df$RG <- factor(Annual2015_Stratified_HarvestEstimates_df$RG, levels = KMA14GroupsPC2)
Annual2015_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2015_Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak")))
Annual2015_Stratified_HarvestEstimates_df$Color <- rep(KMA14Colors, 5)
Annual2015_Stratified_HarvestEstimates_df$Harvest[Annual2015_Stratified_HarvestEstimates_df$Harvest == 0] <- NA
str(Annual2015_Stratified_HarvestEstimates_df)


emf(file ="Figures/All Years/2015 Harvest Bubble Plot.emf", width = 9, height = 5.75, family = "serif", bg = "white")

ggplot(data = Annual2015_Stratified_HarvestEstimates_df, aes(x = RG, y = Fishery, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(name = "Harvest\n(Thousands)", limits = c(0, zmax), breaks = c(1000, 5000, 10000, 20000, seq(50000, 250000, 50000)), range = c(0, 20), labels = c(1, 5, 10, 20, 50, 100, 150, 200, 250)) + 
  scale_color_manual(values = rep(KMA14Colors, 5), guide = FALSE) +
  scale_x_discrete(name = "Reporting Group", labels = KMA14GroupsPC2RowsBubble) +
  scale_y_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak"))) +
  theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = rel(1.3))) +
  theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,0.2,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(-0.5,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.7), angle = 00)) +
  theme(text = element_text(family = "times"))

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2014
Annual2014_Stratified_HarvestEstimates_df <- melt(Annual2014_Stratified_HarvestEstimates)
names(Annual2014_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2014_Stratified_HarvestEstimates_df$RG <- factor(Annual2014_Stratified_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC2))
Annual2014_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2014_Stratified_HarvestEstimates_df$Fishery, levels = c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak"))
Annual2014_Stratified_HarvestEstimates_df$Color <- rep(rev(KMA14Colors), 5)
str(Annual2014_Stratified_HarvestEstimates_df)

ggplot(data = Annual2014_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + geom_point() + scale_size_area(max_size = 25) + scale_color_manual(values = rev(rep(KMA14Colors, 5)))

ggplot(data = Annual2014_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, zmax), breaks = seq(50000, 250000, 50000), range = c(0, 25)) + 
  scale_color_manual(values = rev(rep(KMA14Colors, 5))) +
  ggtitle("2014 Harvest")






# Figure for Report
Annual2014_Stratified_HarvestEstimates_df <- melt(Annual2014_Stratified_HarvestEstimates)
names(Annual2014_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2014_Stratified_HarvestEstimates_df$RG <- factor(Annual2014_Stratified_HarvestEstimates_df$RG, levels = KMA14GroupsPC2)
Annual2014_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2014_Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak")))
Annual2014_Stratified_HarvestEstimates_df$Color <- rep(KMA14Colors, 5)
Annual2014_Stratified_HarvestEstimates_df$Harvest[Annual2014_Stratified_HarvestEstimates_df$Harvest == 0] <- NA
str(Annual2014_Stratified_HarvestEstimates_df)


emf(file ="Figures/All Years/2014 Harvest Bubble Plot.emf", width = 9, height = 5.75, family = "serif", bg = "white")

ggplot(data = Annual2014_Stratified_HarvestEstimates_df, aes(x = RG, y = Fishery, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(name = "Harvest\n(Thousands)", limits = c(0, zmax), breaks = c(1000, 5000, 10000, 20000, seq(50000, 250000, 50000)), range = c(0, 20), labels = c(1, 5, 10, 20, 50, 100, 150, 200, 250)) + 
  scale_color_manual(values = rep(KMA14Colors, 5), guide = FALSE) +
  scale_x_discrete(name = "Reporting Group", labels = KMA14GroupsPC2RowsBubble) +
  scale_y_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak"))) +
  theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = rel(1.3))) +
  theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,0.2,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(-0.5,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.7), angle = 00)) +
  theme(text = element_text(family = "times"))

dev.off()




Annual2016_Stratified_HarvestEstimates_df$Harvest <- NA
Annual2016_Stratified_HarvestEstimates_df$Harvest[1] <- 0

emf(file ="Figures/All Years/Blank Harvest Bubble Plot.emf", width = 9, height = 5.75, family = "serif", bg = "white")

ggplot(data = Annual2016_Stratified_HarvestEstimates_df, aes(x = RG, y = Fishery, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(name = "Harvest\n(Thousands)", limits = c(0, zmax), breaks = c(1000, 5000, 10000, 20000, seq(50000, 250000, 50000)), range = c(0, 20), labels = c(1, 5, 10, 20, 50, 100, 150, 200, 250)) + 
  scale_color_manual(values = rep(KMA14Colors, 5), guide = FALSE) +
  scale_x_discrete(name = "Reporting Group", labels = KMA14GroupsPC2RowsBubble) +
  scale_y_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak"))) +
  theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = rel(1.3))) +
  theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,0.2,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(-0.5,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.7), angle = 00)) +
  theme(text = element_text(family = "times"))

dev.off()




Annual2016_Stratified_HarvestEstimates_df$Harvest <- 100000


emf(file ="Figures/All Years/Color Harvest Bubble Plot.emf", width = 9, height = 5.75, family = "serif", bg = "white")

ggplot(data = Annual2016_Stratified_HarvestEstimates_df, aes(x = RG, y = Fishery, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(name = "Harvest\n(Thousands)", limits = c(0, zmax), breaks = c(1000, 5000, 10000, 20000, seq(50000, 250000, 50000)), range = c(0, 20), labels = c(1, 5, 10, 20, 50, 100, 150, 200, 250)) + 
  scale_color_manual(values = rep(KMA14Colors, 5), guide = FALSE) +
  scale_x_discrete(name = "Reporting Group", labels = KMA14GroupsPC2RowsBubble) +
  scale_y_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak"))) +
  theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = rel(1.3))) +
  theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,0.2,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(-0.5,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.7), angle = 00)) +
  theme(text = element_text(family = "times"))

dev.off()






# Data for Eric Dieters 12/20/16
t(rbind("2014" = Annual2014_Stratified_HarvestEstimates["Upper Station / Akalura", c(5, 6, 4, 2, 1)],
        "2015" = Annual2015_Stratified_HarvestEstimates["Upper Station / Akalura", c(5, 6, 4, 2, 1)],
        "2016" = Annual2016_Stratified_HarvestEstimates["Upper Station / Akalura", c(5, 6, 4, 2, 1)]))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Temporal Strata Bubble Plots
require(ggplot2)
require(reshape2)
require(devEMF)
require(abind)
require(gridExtra)
require(grid)

# Get data
KMA2014Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
KMA2015Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_HarvestEstimatesStats.txt")
KMA2016Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_HarvestEstimatesStats.txt")

KMAStrata_HarvestEstimatesStats <- c(KMA2014Strata_HarvestEstimatesStats,
                                     KMA2015Strata_HarvestEstimatesStats,
                                     KMA2016Strata_HarvestEstimatesStats)
max(sapply(KMAStrata_HarvestEstimatesStats, function(strata) {strata[, "median"]}))
max(sapply(KMAStrata_HarvestEstimatesStats, function(strata) {strata[, "95%"]}))

max(HarvestByStrata2014_Final, HarvestByStrata2015_Final, HarvestByStrata2016_Final, na.rm = TRUE)


KMA2014_Temporal_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Temporal_Stratified_HarvestEstimatesStats.txt")
KMA2015_Temporal_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Temporal_Stratified_HarvestEstimatesStats.txt")
KMA2016_Temporal_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Temporal_Stratified_HarvestEstimatesStats.txt")

KMATemporalStrata_HarvestEstimatesStats <- c(KMA2014_Temporal_Stratified_HarvestEstimatesStats,
                                             KMA2015_Temporal_Stratified_HarvestEstimatesStats,
                                             KMA2016_Temporal_Stratified_HarvestEstimatesStats)
max(sapply(KMATemporalStrata_HarvestEstimatesStats, function(strata) {strata[, "median"]}))
max(sapply(KMATemporalStrata_HarvestEstimatesStats, function(strata) {strata[, "95%"]}))


KMA14GroupsPC2RowsBubble <- KMA14GroupsPC2Rows
KMA14GroupsPC2RowsBubble[c(6,7,9,11,13)] <- gsub(pattern = "\n", replacement = "", x = KMA14GroupsPC2RowsBubble[c(6,7,9,11,13)])


str(KMA2016Strata_HarvestEstimatesStats)

Strata_Bubbleplot.f <- function(yr, strata, zmax = 185118, bubrange = 30) {
  harvest_lst <- get(paste0("KMA", yr, "Strata_HarvestEstimatesStats"))
  strata_nms <- grep(pattern = strata, x = names(harvest_lst), value = TRUE)
  strata_silly <- c("SALIT", "SAYAK", "SIGVA", "SKARL", "SUGAN", "SUYAK")
  
  harvest_median_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
  harvest_median_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "median"])})
  harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
  harvest_median_mat[, harvest_col_index] <- harvest_median_mat.temp
  
  Stratified_HarvestEstimates_df <- melt(harvest_median_mat)
  names(Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
  Stratified_HarvestEstimates_df$RG <- factor(Stratified_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC2))
  Stratified_HarvestEstimates_df$Fishery <- factor(Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak")))
  Stratified_HarvestEstimates_df$Color <- rep(KMA14Colors, 6)
  Stratified_HarvestEstimates_df$Harvest[Stratified_HarvestEstimates_df$Harvest == 0] <- NA
  # str(Stratified_HarvestEstimates_df)
  
  ggplot(data = Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
    geom_point() + 
    scale_size_continuous(name = "Harvest\n(1,000s)", limits = c(0, zmax), breaks = c(1000, 5000, 10000, 20000, seq(50000, 200000, 50000)), range = c(0, bubrange), labels = c(1, 5, 10, 20, 50, 100, 150, 200)) + 
    scale_color_manual(values = rev(KMA14Colors), guide = FALSE) +
    scale_y_discrete(name = "Reporting Group", labels = rev(KMA14GroupsPC2RowsBubble)) +
    scale_x_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak"))) +
    theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = rel(1.3))) +
    theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,-0.5,0,0), "cm"))) +
    theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
    theme(legend.title = element_text(size = rel(1.7), angle = 00)) +
    theme(text = element_text(family = "times"))
}


emf(file ="Figures/Harvest Maps/2014 Early Color Harvest Bubble Plot.emf", width = 6, height = 6.75, family = "serif", bg = "white")
Strata_Bubbleplot.f(yr = 2014, strata = "1_Early", zmax = 185118, bubrange = 20); dev.off()

emf(file ="Figures/Harvest Maps/2014 Middle Color Harvest Bubble Plot.emf", width = 6, height = 6.75, family = "serif", bg = "white")
Strata_Bubbleplot.f(yr = 2014, strata = "2_Middle", zmax = 185118, bubrange = 20); dev.off()

emf(file ="Figures/Harvest Maps/2014 Late Color Harvest Bubble Plot.emf", width = 6, height = 6.75, family = "serif", bg = "white")
Strata_Bubbleplot.f(yr = 2014, strata = "3_Late", zmax = 185118, bubrange = 20); dev.off()

emf(file ="Figures/Harvest Maps/2015 Early Color Harvest Bubble Plot.emf", width = 6, height = 6.75, family = "serif", bg = "white")
Strata_Bubbleplot.f(yr = 2015, strata = "1_Early", zmax = 185118, bubrange = 20); dev.off()

emf(file ="Figures/Harvest Maps/2015 Middle Color Harvest Bubble Plot.emf", width = 6, height = 6.75, family = "serif", bg = "white")
Strata_Bubbleplot.f(yr = 2015, strata = "2_Middle", zmax = 185118, bubrange = 20); dev.off()

emf(file ="Figures/Harvest Maps/2015 Late Color Harvest Bubble Plot.emf", width = 6, height = 6.75, family = "serif", bg = "white")
Strata_Bubbleplot.f(yr = 2015, strata = "3_Late", zmax = 185118, bubrange = 20); dev.off()

emf(file ="Figures/Harvest Maps/2016 Early Color Harvest Bubble Plot.emf", width = 6, height = 6.75, family = "serif", bg = "white")
Strata_Bubbleplot.f(yr = 2016, strata = "1_Early", zmax = 185118, bubrange = 20); dev.off()

emf(file ="Figures/Harvest Maps/2016 Middle Color Harvest Bubble Plot.emf", width = 6, height = 6.75, family = "serif", bg = "white")
Strata_Bubbleplot.f(yr = 2016, strata = "2_Middle", zmax = 185118, bubrange = 20); dev.off()

emf(file ="Figures/Harvest Maps/2016 Late Color Harvest Bubble Plot.emf", width = 6, height = 6.75, family = "serif", bg = "white")
Strata_Bubbleplot.f(yr = 2016, strata = "3_Late", zmax = 185118, bubrange = 20); dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Add marginal totals
# Tyler suggested ggExtra::ggMarginal
ggExtra::ggMarginal(p, type = "density")  # fail

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Make a blank grid

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Credibility intervals ? Try plotting the 5%, median, and 95% all at once with alpha = 0.33 for transparency.
yr = 2015; strata = "2_Middle"; zmax = 202725; bubrange = 20

# Read in stock-specific harvest data
harvest_lst <- get(paste0("KMA", yr, "Strata_HarvestEstimatesStats"))
strata_nms <- grep(pattern = strata, x = names(harvest_lst), value = TRUE)
strata_silly <- c("SALIT", "SAYAK", "SIGVA", "SKARL", "SUGAN", "SUYAK")

# Matrix of medians
harvest_median_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
harvest_median_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "median"])})
harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
harvest_median_mat[, harvest_col_index] <- harvest_median_mat.temp

# Matrix of 5%
harvest_5_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
harvest_5_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "5%"])})
harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
harvest_5_mat[, harvest_col_index] <- harvest_5_mat.temp

# Matrix of 95%
harvest_95_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
harvest_95_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "95%"])})
harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
harvest_95_mat[, harvest_col_index] <- harvest_95_mat.temp

# Create array of 5%, median, and 95% estimates
harvest_array <- abind("lower" = harvest_5_mat, "median" = harvest_median_mat, "upper" = harvest_95_mat, along = 3)

# Melt array into a dataframe for ggplot
Stratified_HarvestEstimates_df <- melt(harvest_array)
names(Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Estimator", "Harvest")
Stratified_HarvestEstimates_df$RG <- factor(Stratified_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC2))
Stratified_HarvestEstimates_df$Fishery <- factor(Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak")))
Stratified_HarvestEstimates_df$Harvest[Stratified_HarvestEstimates_df$Harvest == 0] <- NA
# str(Stratified_HarvestEstimates_df)

bubble_chart <- ggplot(data = Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point(alpha = 0.33) + 
  scale_size_area(name = "Harvest\n(1,000s)", limits = c(0, max(Stratified_HarvestEstimates_df$Harvest)), breaks = c(1000, 5000, 10000, 20000, seq(50000, 200000, 50000)), max_size = bubrange, labels = c(1, 5, 10, 20, 50, 100, 150, 200)) + 
  scale_color_manual(values = rev(KMA14Colors), guide = FALSE) +
  scale_y_discrete(name = "Reporting Group", labels = rev(KMA14GroupsPC2RowsBubble)) +
  scale_x_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak"))) +
  theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = rel(1.3))) +
  theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,-0.5,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.7), angle = 00)) +
  theme(text = element_text(family = "times"))

bubble_chart


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Credibility intervals + marginal totals as bars? Try plotting the 5%, median, and 95% all at once with alpha = 0.33 for transparency.
yr = 2015; strata = "2_Middle"; zmax = 202725; bubrange = 20

Strata_Marginal_Credibility_Bubbleplot.f <- function(yr, strata, zmax = 202725, bubrange = 20) {
  
  # Read in stock-specific harvest data
  harvest_lst <- get(paste0("KMA", yr, "Strata_HarvestEstimatesStats"))
  strata_nms <- grep(pattern = strata, x = names(harvest_lst), value = TRUE)
  strata_silly <- c("SALIT", "SAYAK", "SIGVA", "SKARL", "SUGAN", "SUYAK")
  
  # Matrix of medians
  harvest_median_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
  harvest_median_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "median"])})
  harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
  harvest_median_mat[, harvest_col_index] <- harvest_median_mat.temp
  
  # Matrix of 5%
  harvest_5_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
  harvest_5_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "5%"])})
  harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
  harvest_5_mat[, harvest_col_index] <- harvest_5_mat.temp
  
  # Matrix of 95%
  harvest_95_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
  harvest_95_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "95%"])})
  harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
  harvest_95_mat[, harvest_col_index] <- harvest_95_mat.temp
  
  # Create array of 5%, median, and 95% estimates
  harvest_array <- abind("lower" = harvest_5_mat, "median" = harvest_median_mat, "upper" = harvest_95_mat, along = 3)
  
  # Melt array into a dataframe for ggplot
  Stratified_HarvestEstimates_df <- melt(harvest_array)
  names(Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Estimator", "Harvest")
  Stratified_HarvestEstimates_df$RG <- factor(Stratified_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC2))
  Stratified_HarvestEstimates_df$Fishery <- factor(Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak")))
  Stratified_HarvestEstimates_df$Harvest[Stratified_HarvestEstimates_df$Harvest == 0] <- NA
  # str(Stratified_HarvestEstimates_df)
  
  # Matrix of means
  harvest_mean_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
  harvest_mean_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "mean"])})
  harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
  harvest_mean_mat[, harvest_col_index] <- harvest_mean_mat.temp
  harvest_mean_mat <- addmargins(harvest_mean_mat)
  
  # Melt means into dataframe for barplot
  Total_HarvestEstimates_df <- melt(harvest_mean_mat["Sum", rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak")), drop = FALSE])
  names(Total_HarvestEstimates_df) <- c("Estimator", "Fishery", "Harvest")
  Total_HarvestEstimates_df$Harvest[Total_HarvestEstimates_df$Harvest == 0] <- NA
  
  # Matrix of temporal medians
  temporal_harvest_mat <- get(paste0("KMA", yr, "_Temporal_Stratified_HarvestEstimatesStats"))[[strata]]
  Temporal_HarvestEstimates_df <- data.frame("RG" = KMA14GroupsPC2, temporal_harvest_mat)[, c(1, 4:6)]
  names(Temporal_HarvestEstimates_df) <- c("RG", "Median", "Lower5", "Upper95")
  Temporal_HarvestEstimates_df$RG <- factor(Temporal_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC2))
  
  bar_top <- ggplot(data = Total_HarvestEstimates_df, aes(x = Fishery, y = Harvest)) + 
    geom_bar(stat = "identity") +
    scale_y_continuous(name = "Harvest (1,000s of fish)", breaks = seq(0, 400000, by = 100000), labels = c(seq(0, 400, by = 100)), limits = c(0, 400000)) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    theme(axis.text.y = element_text(size = rel(1)), axis.title.y = element_text(size = rel(1.3), margin = unit(c(0,-3,0,0), "cm")))
  
  # empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  #   theme(axis.ticks=element_blank(), 
  #         panel.background=element_blank(), 
  #         axis.text.x=element_blank(), axis.text.y=element_blank(),           
  #         axis.title.x=element_blank(), axis.title.y=element_blank())
  
  bubble_chart <- ggplot(data = Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
    geom_point(alpha = 0.33) + 
    scale_size_area(name = "Harvest\n(1,000s)", limits = c(0, max(Stratified_HarvestEstimates_df$Harvest)), breaks = c(1000, 5000, 10000, 20000, seq(50000, 200000, 50000)), max_size = bubrange, labels = c(1, 5, 10, 20, 50, 100, 150, 200), guide = "legend") + 
    scale_color_manual(values = rev(KMA14Colors), guide = FALSE) +
    scale_y_discrete(name = "Reporting Group", labels = rev(KMA14GroupsPC2RowsBubble)) +
    scale_x_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak"))) +
    theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = rel(1.3))) +
    theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,-0.5,0,0), "cm"))) +
    theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
    theme(legend.title = element_text(size = rel(1.7), angle = 00), legend.position = "none") +
    guides(size = guide_legend(ncol = 2, byrow = TRUE)) +
    theme(text = element_text(family = "times"))
  
  bar_right <- ggplot(data = Temporal_HarvestEstimates_df, aes(x = RG, y = Median, ymin = Lower5, ymax = Upper95)) + 
    geom_bar(stat = "identity", show.legend = FALSE, fill = rev(KMA14Colors)) +
    geom_errorbar(width = 0.5) +
    scale_y_continuous(name = "Harvest (1,000s of fish)", breaks = seq(0, 500000, by = 100000), labels = c(seq(0, 500, by = 100)), limits = c(0, 500000)) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = rel(1)), axis.title.x = element_text(size = rel(1.3), vjust = 20)) +
    coord_flip()
  
  legend_topright_feeder <- ggplot(data = Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
    geom_point(alpha = 0.33) + 
    scale_size_area(name = "Harvest\n(1,000s)", limits = c(0, zmax), breaks = c(1000, 5000, 10000, 20000, seq(50000, 200000, 50000)), max_size = bubrange, labels = c(1, 5, 10, 20, 50, 100, 150, 200), guide = "legend") + 
    scale_color_manual(values = rev(KMA14Colors), guide = FALSE) +
    scale_y_discrete(name = "Reporting Group", labels = rev(KMA14GroupsPC2RowsBubble)) +
    scale_x_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak"))) +
    theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = rel(1.3))) +
    theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,-0.5,0,0), "cm"))) +
    theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
    theme(legend.title = element_text(size = rel(1.7), angle = 00), legend.position = "right") +
    guides(size = guide_legend(ncol = 2, byrow = TRUE)) +
    theme(text = element_text(family = "times"))
  
  g_legend<-function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)} 
  
  legend_topright <- g_legend(legend_topright_feeder)
  
  # Plot all four, but not aligned
  # grid.arrange(bar_top, legend_topright, bubble_chart, bar_right, ncol=2, nrow=2, widths=c(4, 1.5), heights=c(1.5, 4))
  
  #~~~~~~~~~~~~~~~~~~
  # Convert to tables and then 
  p1 <- ggplot_gtable(ggplot_build(bar_top))
  # p2 <- ggplot_gtable(ggplot_build(legend_topright))
  p3 <- ggplot_gtable(ggplot_build(bubble_chart))
  p4 <- ggplot_gtable(ggplot_build(bar_right))
  
  maxWidth <- unit.pmax(p1$widths[2:3], p3$widths[2:3])
  p1$widths[2:3] <- maxWidth
  p3$widths[2:3] <- maxWidth
  
  maxHeight <- unit.pmax(p3$heights[7], p4$heights[7])
  p3$heights[7] <- maxHeight
  p4$heights[7] <- maxHeight
  
  grid.arrange(p1, legend_topright, p3, p4, ncol=2, nrow=2, widths=c(4, 1.35), heights=c(1.35, 4))
}


png(file ="Figures/Harvest Maps/2014 1_Early Marginal Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Bubbleplot.f(yr = 2014, strata = "1_Early"); dev.off()
png(file ="Figures/Harvest Maps/2015 1_Early Marginal Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Bubbleplot.f(yr = 2015, strata = "1_Early"); dev.off()
png(file ="Figures/Harvest Maps/2016 1_Early Marginal Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Bubbleplot.f(yr = 2016, strata = "1_Early"); dev.off()

png(file ="Figures/Harvest Maps/2014 2_Middle Marginal Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Bubbleplot.f(yr = 2014, strata = "2_Middle"); dev.off()
png(file ="Figures/Harvest Maps/2015 2_Middle Marginal Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Bubbleplot.f(yr = 2015, strata = "2_Middle"); dev.off()
png(file ="Figures/Harvest Maps/2016 2_Middle Marginal Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Bubbleplot.f(yr = 2016, strata = "2_Middle"); dev.off()

png(file ="Figures/Harvest Maps/2014 3_Late Marginal Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Bubbleplot.f(yr = 2014, strata = "3_Late"); dev.off()
png(file ="Figures/Harvest Maps/2015 3_Late Marginal Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Bubbleplot.f(yr = 2015, strata = "3_Late"); dev.off()
png(file ="Figures/Harvest Maps/2016 3_Late Marginal Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Bubbleplot.f(yr = 2016, strata = "3_Late"); dev.off()
#~~~~~~~~~~~~~~~~~~
# Extract Legend 
# g_legend<-function(a.gplot){ 
#   tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
#   legend <- tmp$grobs[[leg]] 
#   return(legend)} 
# 
# legend_topright <- g_legend(bubble_chart)
# grid.draw(legend_topright)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Credibility intervals + marginal totals as bars, including unsampled areas and post-sampling? Try plotting the 5%, median, and 95% all at once with alpha = 0.33 for transparency.
# Need to run KMA_Harvest_Simple.f first to get s14, s15, and s16 objects.

topmax <- max(s14, s15, s16)

yr = 2015; strata = "2_Middle"; zmax = 202725; bubrange = 20

Strata_Marginal_Credibility_Unsamp_Bubbleplot.f <- function(yr, strata, zmax = 202725, bubrange = 20) {
  
  if(strata == "4_Post") {
    # Get harvest data
    Total_HarvestEstimates_df <- get(paste0("s", yr - 2000))[, c("Geo", unlist(strsplit(x = strata, split = "_"))[2])]
    names(Total_HarvestEstimates_df) <- c("Fishery", "Harvest")
    levels(Total_HarvestEstimates_df$Fishery)
    Total_HarvestEstimates_df$Harvest[Total_HarvestEstimates_df$Harvest == 0] <- NA
    
    # Create blank bubble plot
    harvest_array <- array(data = c(0, rep(NA, 293)), dim = c(14, 7, 3), dimnames = list(KMA14GroupsPC2, 
                                                                        rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak", "Unsampled")),
                                                                        c("lower", "median", "upper")))
    Stratified_HarvestEstimates_df <- melt(harvest_array)
    names(Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Estimator", "Harvest")
    Stratified_HarvestEstimates_df$RG <- factor(Stratified_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC2))
    Stratified_HarvestEstimates_df$Fishery <- factor(Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak", "Unsampled")))

    # Create blank right plot
    Temporal_HarvestEstimates_df <- data.frame("RG" = factor(x = KMA14GroupsPC2, levels = KMA14GroupsPC2),
                                               "Median" = rep(0, 14),
                                               "Lower5" = rep(0, 14),
                                               "Upper95" = rep(0, 14))
  } else {
    # Read in stock-specific harvest data
    harvest_lst <- get(paste0("KMA", yr, "Strata_HarvestEstimatesStats"))
    strata_nms <- grep(pattern = strata, x = names(harvest_lst), value = TRUE)
    strata_silly <- c("SALIT", "SAYAK", "SIGVA", "SKARL", "SUGAN", "SUYAK")
    
    # Matrix of medians
    harvest_median_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
    harvest_median_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "median"])})
    harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
    harvest_median_mat[, harvest_col_index] <- harvest_median_mat.temp
    harvest_median_mat <- cbind("Unsampled" = rep(0, 14), harvest_median_mat)
    
    # Matrix of 5%
    harvest_5_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
    harvest_5_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "5%"])})
    harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
    harvest_5_mat[, harvest_col_index] <- harvest_5_mat.temp
    harvest_5_mat <- cbind("Unsampled" = rep(0, 14), harvest_5_mat)
    
    # Matrix of 95%
    harvest_95_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
    harvest_95_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "95%"])})
    harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
    harvest_95_mat[, harvest_col_index] <- harvest_95_mat.temp
    harvest_95_mat <- cbind("Unsampled" = rep(0, 14), harvest_95_mat)
    
    # Create array of 5%, median, and 95% estimates
    harvest_array <- abind("lower" = harvest_5_mat, "median" = harvest_median_mat, "upper" = harvest_95_mat, along = 3)
    
    # Melt array into a dataframe for ggplot
    Stratified_HarvestEstimates_df <- melt(harvest_array)
    names(Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Estimator", "Harvest")
    Stratified_HarvestEstimates_df$RG <- factor(Stratified_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC2))
    Stratified_HarvestEstimates_df$Fishery <- factor(Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak", "Unsampled")))
    Stratified_HarvestEstimates_df$Harvest[Stratified_HarvestEstimates_df$Harvest == 0] <- NA
    # str(Stratified_HarvestEstimates_df)
    
    # Get harvest data
    Total_HarvestEstimates_df <- get(paste0("s", yr - 2000))[, c("Geo", unlist(strsplit(x = strata, split = "_"))[2])]
    names(Total_HarvestEstimates_df) <- c("Fishery", "Harvest")
    levels(Total_HarvestEstimates_df$Fishery)
    Total_HarvestEstimates_df$Harvest[Total_HarvestEstimates_df$Harvest == 0] <- NA
    
    # Matrix of temporal medians
    temporal_harvest_mat <- get(paste0("KMA", yr, "_Temporal_Stratified_HarvestEstimatesStats"))[[strata]]
    Temporal_HarvestEstimates_df <- data.frame("RG" = KMA14GroupsPC2, temporal_harvest_mat)[, c(1, 4:6)]
    names(Temporal_HarvestEstimates_df) <- c("RG", "Median", "Lower5", "Upper95")
    Temporal_HarvestEstimates_df$RG <- factor(Temporal_HarvestEstimates_df$RG, levels = rev(KMA14GroupsPC2))
  }
  
  bar_top <- ggplot(data = Total_HarvestEstimates_df, aes(x = Fishery, y = Harvest)) + 
    geom_bar(stat = "identity") +
    scale_y_continuous(name = "Harvest (1,000s of fish)", breaks = seq(0, 500000, by = 100000), labels = c(seq(0, 500, by = 100)), limits = c(0, topmax)) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    theme(axis.text.y = element_text(size = rel(1)), axis.title.y = element_text(size = rel(1.3), margin = unit(c(0,-3,0,0), "cm")))
  
  # empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  #   theme(axis.ticks=element_blank(), 
  #         panel.background=element_blank(), 
  #         axis.text.x=element_blank(), axis.text.y=element_blank(),           
  #         axis.title.x=element_blank(), axis.title.y=element_blank())
  
  bubble_chart <- ggplot(data = Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
    geom_point(alpha = 0.33) + 
    scale_size_area(name = "Harvest\n(1,000s)", limits = c(0, max(Stratified_HarvestEstimates_df$Harvest)), breaks = c(1000, 5000, 10000, 20000, seq(50000, 200000, 50000)), max_size = bubrange, labels = c(1, 5, 10, 20, 50, 100, 150, 200), guide = "legend") + 
    scale_color_manual(values = rev(KMA14Colors), guide = FALSE) +
    scale_y_discrete(name = "Reporting Group", labels = rev(KMA14GroupsPC2RowsBubble)) +
    scale_x_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak", "Unsampled\nAreas"))) +
    theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = rel(1.3))) +
    theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,-0.5,0,0), "cm"))) +
    theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
    theme(legend.title = element_text(size = rel(1.7), angle = 00), legend.position = "none") +
    guides(size = guide_legend(ncol = 2, byrow = TRUE)) +
    theme(text = element_text(family = "times"))
  
  bar_right <- ggplot(data = Temporal_HarvestEstimates_df, aes(x = RG, y = Median, ymin = Lower5, ymax = Upper95)) + 
    geom_bar(stat = "identity", show.legend = FALSE, fill = rev(KMA14Colors)) +
    geom_errorbar(width = 0.5) +
    scale_y_continuous(name = "Harvest (1,000s of fish)", breaks = seq(0, 500000, by = 100000), labels = c(seq(0, 500, by = 100)), limits = c(0, 500000)) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = rel(1)), axis.title.x = element_text(size = rel(1.3), vjust = 20)) +
    coord_flip()
  
  legend_topright_feeder <- ggplot(data = Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
    geom_point(alpha = 0.33) + 
    scale_size_area(name = "Harvest\n(1,000s)", limits = c(0, zmax), breaks = c(1000, 5000, 10000, 20000, seq(50000, 200000, 50000)), max_size = bubrange, labels = c(1, 5, 10, 20, 50, 100, 150, 200), guide = "legend") + 
    scale_color_manual(values = rev(KMA14Colors), guide = FALSE) +
    scale_y_discrete(name = "Reporting Group", labels = rev(KMA14GroupsPC2RowsBubble)) +
    scale_x_discrete(name = "Sampling Area", labels = rev(c("Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak", "Unsampled\nAreas"))) +
    theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = rel(1.3))) +
    theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,-0.5,0,0), "cm"))) +
    theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
    theme(legend.title = element_text(size = rel(1.7), angle = 00), legend.position = "right") +
    guides(size = guide_legend(ncol = 2, byrow = TRUE)) +
    theme(text = element_text(family = "times"))
  
  g_legend<-function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)} 
  
  legend_topright <- g_legend(legend_topright_feeder)
  
  # Plot all four, but not aligned
  # grid.arrange(bar_top, legend_topright, bubble_chart, bar_right, ncol=2, nrow=2, widths=c(4, 1.5), heights=c(1.5, 4))
  
  #~~~~~~~~~~~~~~~~~~
  # Convert to tables and then 
  p1 <- ggplot_gtable(ggplot_build(bar_top))
  # p2 <- ggplot_gtable(ggplot_build(legend_topright))
  p3 <- ggplot_gtable(ggplot_build(bubble_chart))
  p4 <- ggplot_gtable(ggplot_build(bar_right))
  
  maxWidth <- unit.pmax(p1$widths[2:3], p3$widths[2:3])
  p1$widths[2:3] <- maxWidth
  p3$widths[2:3] <- maxWidth
  
  maxHeight <- unit.pmax(p3$heights[7], p4$heights[7])
  p3$heights[7] <- maxHeight
  p4$heights[7] <- maxHeight
  
  grid.arrange(p1, legend_topright, p3, p4, ncol=2, nrow=2, widths=c(4, 1.35), heights=c(1.35, 4))
}

png(file ="Figures/Harvest Maps/2014 1_Early Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2014, strata = "1_Early"); dev.off()
png(file ="Figures/Harvest Maps/2014 2_Middle Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2014, strata = "2_Middle"); dev.off()
png(file ="Figures/Harvest Maps/2014 3_Late Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2014, strata = "3_Late"); dev.off()
png(file ="Figures/Harvest Maps/2014 4_Post Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2014, strata = "4_Post"); dev.off()

png(file ="Figures/Harvest Maps/2015 1_Early Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2015, strata = "1_Early"); dev.off()
png(file ="Figures/Harvest Maps/2015 2_Middle Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2015, strata = "2_Middle"); dev.off()
png(file ="Figures/Harvest Maps/2015 3_Late Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2015, strata = "3_Late"); dev.off()
png(file ="Figures/Harvest Maps/2015 4_Post Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2015, strata = "4_Post"); dev.off()

png(file ="Figures/Harvest Maps/2016 1_Early Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2016, strata = "1_Early"); dev.off()
png(file ="Figures/Harvest Maps/2016 2_Middle Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2016, strata = "2_Middle"); dev.off()
png(file ="Figures/Harvest Maps/2016 3_Late Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2016, strata = "3_Late"); dev.off()
png(file ="Figures/Harvest Maps/2016 4_Post Marginal Unsamp Bubble.png", width = 811, height = 957, units = "px", res = 96, family = "serif", bg = "white")
Strata_Marginal_Credibility_Unsamp_Bubbleplot.f(yr = 2016, strata = "4_Post"); dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Credibility intervals + marginal totals as bubbles? Try plotting the 5%, median, and 95% all at once with alpha = 0.33 for transparency.
yr = 2015; strata = "2_Middle"; zmax = 202725; bubrange = 20

harvest_lst <- get(paste0("KMA", yr, "Strata_HarvestEstimatesStats"))
strata_nms <- grep(pattern = strata, x = names(harvest_lst), value = TRUE)
strata_silly <- c("SALIT", "SAYAK", "SIGVA", "SKARL", "SUGAN", "SUYAK")

harvest_median_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
harvest_median_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "median"])})
harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
harvest_median_mat[, harvest_col_index] <- harvest_median_mat.temp
harvest_median_mat <- addmargins(harvest_median_mat)

harvest_5_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
harvest_5_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "5%"])})
harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
harvest_5_mat[, harvest_col_index] <- harvest_5_mat.temp
harvest_5_mat <- addmargins(harvest_5_mat)

harvest_95_mat <- matrix(data = 0, nrow = 14, ncol = 6, dimnames = list(KMA14GroupsPC2, c("Alitak", "Ayakulik", "Igvak", "Karluk", "Uganik", "Uyak")))
harvest_95_mat.temp <- sapply(strata_nms, function(strat) {round(harvest_lst[[strat]][, "95%"])})
harvest_col_index <- sapply(strata_nms, function(strat) {which(strata_silly == paste(unlist(strsplit(x = strat, split = ''))[1:5], collapse = ''))})
harvest_95_mat[, harvest_col_index] <- harvest_95_mat.temp
harvest_95_mat <- addmargins(harvest_95_mat)


harvest_array <- abind("lower" = harvest_5_mat, "median" = harvest_median_mat, "upper" = harvest_95_mat, along = 3)


Stratified_HarvestEstimates_df <- melt(harvest_array)
names(Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Estimator", "Harvest")
Stratified_HarvestEstimates_df$RG <- factor(Stratified_HarvestEstimates_df$RG, levels = c("Sum", rev(KMA14GroupsPC2)))
Stratified_HarvestEstimates_df$Fishery <- factor(Stratified_HarvestEstimates_df$Fishery, levels = rev(c("Sum", "Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak")))
Stratified_HarvestEstimates_df$Harvest[Stratified_HarvestEstimates_df$Harvest == 0] <- NA
# str(Stratified_HarvestEstimates_df)

ggplot(data = Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point(alpha = 0.33) + 
  scale_size_continuous(name = "Harvest\n(1,000s)", limits = c(0, max(Stratified_HarvestEstimates_df$Harvest)), breaks = c(1000, 5000, 10000, 20000, seq(50000, 200000, 50000)), range = c(0, bubrange), labels = c(1, 5, 10, 20, 50, 100, 150, 200)) + 
  scale_color_manual(values = c("grey80", rev(KMA14Colors)), guide = FALSE) +
  scale_y_discrete(name = "Reporting Group", labels = c("Total", rev(KMA14GroupsPC2RowsBubble))) +
  scale_x_discrete(name = "Sampling Area", labels = rev(c("Total", "Uganik\nKupreanof", "Uyak", "Karluk\nSturgeon", "Ayakulik\nHalibut Bay", "Alitak", "Igvak"))) +
  theme(axis.text.x = element_text(size = rel(1.2), angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = rel(1.3))) +
  theme(axis.title.y = element_text(size = rel(1.7), angle = 90, margin = unit(c(0,-0.5,0,0), "cm"))) +
  theme(axis.title.x = element_text(size = rel(1.7), angle = 00, margin = unit(c(0,0,0,0), "cm"))) +
  theme(legend.title = element_text(size = rel(1.7), angle = 00)) +
  theme(text = element_text(family = "times"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Exvessel value


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Percentages for KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")
KMA2015Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")
KMA2016Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_EstimatesStats.txt")

str(KMA2014Strata_EstimatesStats)

sapply(c(KMA2014Strata_EstimatesStats, KMA2015Strata_EstimatesStats, KMA2016Strata_EstimatesStats), function(mix) {table(mix[, "GR"] > 1.2)})
round(KMA2014Strata_EstimatesStats$SUYAKC14_2_Middle[, c("mean", "GR")], 3)
round(KMA2015Strata_EstimatesStats$SUGANC15_3_Late[, c("mean", "GR")], 3)
round(KMA2016Strata_EstimatesStats$SKARLC16_1_Early[, c("mean", "GR")], 3)



# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             1, 3,
                             1, 4,
                             5, 6), nrow = 4, ncol = 2, byrow = TRUE)



GeoMix <- c("SUGANC", "SUYAKC", "SKARLC", "SAYAKC", "SALITC", "SIGVAC")

filenames <- setNames(object = c("Uganik Proportions 2014-2016", 
                                 "Uyak Proportions 2014-2016",
                                 "Karluk Proportions 2014-2016",
                                 "Ayakulik Proportions 2014-2016",
                                 "Alitak Proportions 2014-2016",
                                 "Igvak Proportions 2014-2016"), nm = GeoMix)

# If showing proportions (percetages) use blue, otherwise green as "low"
ProportionColors <- colorpanel(n = 3, low = "blue", high = "white")


#~~~~~~~~~~~~~~~~~~
# 2014
TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend14 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend14 <- sapply(KMA2014, function(geo) {
  Legend14[sapply(TempMix14[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempProportionColors14 <- sapply(KMA2014, function(geo) {
  ProportionColors[sapply(TempMix14[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates14 <- KMA2014Strata_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2015
TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend15 <- setNames(object = c("June 1-July 3", "July 4-August 1", "August 2-29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend15 <- sapply(KMA2015, function(geo) {
  Legend15[sapply(TempMix15[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempProportionColors15 <- sapply(KMA2015, function(geo) {
  ProportionColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates15 <- KMA2015Strata_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2016
TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend16 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"),
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend16 <- sapply(KMA2016, function(geo) {
  Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempProportionColors16 <- sapply(KMA2016, function(geo) {
  ProportionColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates16 <- KMA2016Strata_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)

sapply(GeoMix, function(geomix) {
  
  emf(file = paste("Figures/All Years/", filenames[geomix], ".emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")
  
  
  layout(mat = layoutmat, widths = c(0.075, 1, 1), heights = c(0.9, 0.9, 0.9, 0.15))
  par(mar = rep(0, 4))
  par(family = "times")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Y-axis label
  plot.new()
  text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2014 Barplot
  if(length(grep(pattern = geomix, x = names(TempMix14), value = TRUE)) == 0) {
    par(mar = c(1, 1, 1, 1))
    Barplot14 <- barplot2(height = rep(0, 14), beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd, ci.l = rep(0, 14),  ci.u = rep(0, 14), 
                          ylim = c(0, 100), col = 1, yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
    abline(h = 0, xpd = FALSE)
    text(x = 8.5, y = 50, labels = "No estimates avaialable for 2014", cex = cex.leg)
  } else {
    geomix14 <- grep(pattern = geomix, x = names(TempMix14), value = TRUE)
    par(mar = c(1, 1, 1, 1))
    Barplot14 <- barplot2(height = t(sapply(TempMix14[[geomix14]], function(tempmix) {Estimates14[[tempmix]][, "median"]})) * 100, 
                          beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                          ci.l = t(sapply(TempMix14[[geomix14]], function(tempmix) {Estimates14[[tempmix]][, "5%"]})) * 100, 
                          ci.u = t(sapply(TempMix14[[geomix14]], function(tempmix) {Estimates14[[tempmix]][, "95%"]})) * 100, 
                          ylim = c(0, 100), col = TempProportionColors14[[geomix14]], yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
    legend(legend = TempLegend14[[geomix14]], x = "topleft", fill = TempProportionColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2014")
    abline(h = 0, xpd = FALSE)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2015 Barplot
  geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
  par(mar = c(1, 1, 1, 1))
  Barplot15 <- barplot2(height = t(sapply(TempMix15[[geomix15]], function(tempmix) {Estimates15[[tempmix]][, "median"]})) * 100, 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix15[[geomix15]], function(tempmix) {Estimates15[[tempmix]][, "5%"]})) * 100, 
                        ci.u = t(sapply(TempMix15[[geomix15]], function(tempmix) {Estimates15[[tempmix]][, "95%"]})) * 100, 
                        ylim = c(0, 100), col = TempProportionColors15[[geomix15]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix15]], x = "topleft", fill = TempProportionColors15[[geomix15]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2016 Barplot
  geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
  par(mar = c(1, 1, 1, 1))
  Barplot16 <- barplot2(height = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][, "median"]})) * 100,
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][, "5%"]})) * 100,
                        ci.u = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][, "95%"]})) * 100,
                        ylim = c(0, 100), col = TempProportionColors16[[geomix16]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  if(geomix == "SIGVAC") {
    legend(legend = TempLegend16[[geomix16]], x = "top", fill = TempProportionColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")      
  } else {
    legend(legend = TempLegend16[[geomix16]], x = "topleft", fill = TempProportionColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")  
  }
  abline(h = 0, xpd = FALSE)
  
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot16, 2, mean), adj = 0.5, cex = cex.xaxis)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Blank Corner
  par(mar = rep(0, 4))
  plot.new()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## x-axis label
  par(mar = rep(0, 4))
  plot.new()
  text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)
  
  
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Harvests for KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_HarvestEstimatesStats.txt")
KMA2015_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_HarvestEstimatesStats.txt")
KMA2016_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_HarvestEstimatesStats.txt")

# What is the total annual harvest per RG?
KMA2014_HarvestByRG <- apply(sapply(KMA2014_Annual_HarvestEstimatesStats, function(geo) {geo[, "median"]}), 1, sum)
KMA2015_HarvestByRG <- apply(sapply(KMA2015_Annual_HarvestEstimatesStats, function(geo) {geo[, "median"]}), 1, sum)
KMA2016_HarvestByRG <- apply(sapply(KMA2016_Annual_HarvestEstimatesStats, function(geo) {geo[, "median"]}), 1, sum)

KMA2014_HarvestByRG_LCI <- apply(sapply(KMA2014_Annual_HarvestEstimatesStats, function(geo) {geo[, "5%"]}), 1, sum)
KMA2015_HarvestByRG_LCI <- apply(sapply(KMA2015_Annual_HarvestEstimatesStats, function(geo) {geo[, "5%"]}), 1, sum)
KMA2016_HarvestByRG_LCI <- apply(sapply(KMA2016_Annual_HarvestEstimatesStats, function(geo) {geo[, "5%"]}), 1, sum)

KMA2014_HarvestByRG_UCI <- apply(sapply(KMA2014_Annual_HarvestEstimatesStats, function(geo) {geo[, "95%"]}), 1, sum)
KMA2015_HarvestByRG_UCI <- apply(sapply(KMA2015_Annual_HarvestEstimatesStats, function(geo) {geo[, "95%"]}), 1, sum)
KMA2016_HarvestByRG_UCI <- apply(sapply(KMA2016_Annual_HarvestEstimatesStats, function(geo) {geo[, "95%"]}), 1, sum)

# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.1
cex.leg <- 1.1
ci.lwd <- 2.5
ymax <- max(c(KMA2014_HarvestByRG_UCI, KMA2015_HarvestByRG_UCI, KMA2016_HarvestByRG_UCI))


par(mar = c(4.6, 4.6, 1, 1))
Barplot <- barplot2(height = rbind(KMA2014_HarvestByRG, KMA2015_HarvestByRG, KMA2016_HarvestByRG), 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = rbind(KMA2014_HarvestByRG_LCI, KMA2015_HarvestByRG_LCI, KMA2016_HarvestByRG_LCI), 
                    ci.u = rbind(KMA2014_HarvestByRG_UCI, KMA2015_HarvestByRG_UCI, KMA2016_HarvestByRG_UCI), 
                    ylim = c(0, ymax), col = c("white", "cyan", "blue"), yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 50000), labels = formatC(x = seq(0, ymax, 50000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
mtext(text = "Number of Fish Harvested (Thousands)", side = 2, line = 3, cex = cex.lab)
mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = cex.xaxis)
mtext(text = "Reporting Group", side = 1, line = 3, cex = cex.lab)
mtext(text = "Annual Stock-Specific Harvest\nfor Select KMA Fisheries", side = 3, cex = cex.lab, line = -2)
legend(legend = 2014:2016, x = "topleft", fill = c("white", "cyan", "blue"), border = "black", bty = "n", cex = cex.leg, title="")
abline(h = 0, xpd = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get data
KMA2014Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
KMA2015Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_HarvestEstimatesStats.txt")
KMA2016Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_HarvestEstimatesStats.txt")

str(KMA2014Strata_EstimatesStats)




# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             1, 3,
                             1, 4,
                             5, 6), nrow = 4, ncol = 2, byrow = TRUE)



GeoMix <- c("SUGANC", "SUYAKC", "SKARLC", "SAYAKC", "SALITC", "SIGVAC")

filenames <- setNames(object = c("Uganik Harvest 2014-2016", 
                                 "Uyak Harvest 2014-2016",
                                 "Karluk Harvest 2014-2016",
                                 "Ayakulik Harvest 2014-2016",
                                 "Alitak Harvest 2014-2016",
                                 "Igvak Harvest 2014-2016"), nm = GeoMix)

# If showing proportions (percetages) use blue, otherwise green as "low"
HarvestColors <- colorpanel(n = 3, low = "green", high = "white")


#~~~~~~~~~~~~~~~~~~
# 2014
TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_HarvestEstimatesStats), value = TRUE)}, simplify = FALSE)

Legend14 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend14 <- sapply(KMA2014, function(geo) {
  Legend14[sapply(TempMix14[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempHarvestColors14 <- sapply(KMA2014, function(geo) {
  HarvestColors[sapply(TempMix14[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

HarvestEstimates14 <- KMA2014Strata_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2015
TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_HarvestEstimatesStats), value = TRUE)}, simplify = FALSE)

Legend15 <- setNames(object = c("June 1-July 3", "July 4-August 1", "August 2-29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend15 <- sapply(KMA2015, function(geo) {
  Legend15[sapply(TempMix15[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempHarvestColors15 <- sapply(KMA2015, function(geo) {
  HarvestColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

HarvestEstimates15 <- KMA2015Strata_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2016
TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend16 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"),
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend16 <- sapply(KMA2016, function(geo) {
  Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempHarvestColors16 <- sapply(KMA2016, function(geo) {
  HarvestColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

HarvestEstimates16 <- KMA2016Strata_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5
ymax <- 200000  # max(sapply(c(HarvestEstimates14, HarvestEstimates15, HarvestEstimates16), function(strata) {strata[, "95%"]}))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)

sapply(GeoMix, function(geomix) {
  
  emf(file = paste("Figures/All Years/", filenames[geomix], ".emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")
  
  
  layout(mat = layoutmat, widths = c(0.075, 1, 1), heights = c(0.9, 0.9, 0.9, 0.15))
  par(mar = rep(0, 4))
  par(family = "times")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Y-axis label
  plot.new()
  text(x = 0.25, y = 0.5, labels = "Number of Fish Harvested (Thousands)", srt = 90, cex = cex.lab)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2014 Barplot
  if(length(grep(pattern = geomix, x = names(TempMix14), value = TRUE)) == 0) {
    par(mar = c(1, 1, 1, 1))
    Barplot14 <- barplot2(height = rep(0, 14), beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd, ci.l = rep(0, 14),  ci.u = rep(0, 14), 
                          ylim = c(0, ymax), col = 1, yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, ymax, 50000), labels = formatC(x = seq(0, ymax, 50000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
    abline(h = 0, xpd = FALSE)
    text(x = 8.5, y = 100000, labels = "No estimates avaialable for 2014", cex = cex.leg)
  } else {
    geomix14 <- grep(pattern = geomix, x = names(TempMix14), value = TRUE)
    par(mar = c(1, 1, 1, 1))
    Barplot14 <- barplot2(height = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][, "median"]})), 
                          beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                          ci.l = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][, "5%"]})), 
                          ci.u = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][, "95%"]})), 
                          ylim = c(0, ymax), col = TempHarvestColors14[[geomix14]], yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, ymax, 50000), labels = formatC(x = seq(0, ymax, 50000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
    legend(legend = TempLegend14[[geomix14]], x = "topleft", fill = TempHarvestColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2014")
    abline(h = 0, xpd = FALSE)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2015 Barplot
  geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
  par(mar = c(1, 1, 1, 1))
  Barplot15 <- barplot2(height = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][, "median"]})), 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][, "5%"]})), 
                        ci.u = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][, "95%"]})), 
                        ylim = c(0, ymax), col = TempHarvestColors15[[geomix15]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 50000), labels = formatC(x = seq(0, ymax, 50000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix15]], x = "topleft", fill = TempHarvestColors15[[geomix15]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2016 Barplot
  geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
  par(mar = c(1, 1, 1, 1))
  Barplot16 <- barplot2(height = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][, "median"]})),
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][, "5%"]})),
                        ci.u = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][, "95%"]})),
                        ylim = c(0, ymax), col = TempHarvestColors16[[geomix16]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 50000), labels = formatC(x = seq(0, ymax, 50000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend16[[geomix16]], x = "topleft", fill = TempHarvestColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot16, 2, mean), adj = 0.5, cex = cex.xaxis)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Blank Corner
  par(mar = rep(0, 4))
  plot.new()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## x-axis label
  par(mar = rep(0, 4))
  plot.new()
  text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)
  
  
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dates
# DatesStrata2014 <- read.table(file = "Harvest/2014DatesByStrata.txt", header = TRUE, sep = "\t", as.is = TRUE)
# DatesStrata2014.mat <- as.matrix(DatesStrata2014[-1])
# dimnames(DatesStrata2014.mat) <- list(DatesStrata2014$location, c("1_Early", "2_Middle", "3_Late"))
# dput(x = DatesStrata2014.mat, file = "Objects/DatesStrata2014_Final.txt"); rm(DatesStrata2014.mat)
DatesStrata2014_Final <- dget(file = "Objects/DatesStrata2014_Final.txt")


# DatesStrata2015 <- read.table(file = "Harvest/2015DatesByStrata.txt", header = TRUE, sep = "\t", as.is = TRUE)
# DatesStrata2015.mat <- as.matrix(DatesStrata2015[-1])
# dimnames(DatesStrata2015.mat) <- list(DatesStrata2015$location, c("1_Early", "2_Middle", "3_Late"))
# dput(x = DatesStrata2015.mat, file = "Objects/DatesStrata2015_Final.txt"); rm(DatesStrata2015.mat)
DatesStrata2015_Final <- dget(file = "Objects/DatesStrata2015_Final.txt")


# DatesStrata2016 <- read.table(file = "Harvest/2016DatesByStrata.txt", header = TRUE, sep = "\t", as.is = TRUE)
# DatesStrata2016.mat <- as.matrix(DatesStrata2016[-1])
# dimnames(DatesStrata2016.mat) <- list(DatesStrata2016$location, c("1_Early", "2_Middle", "3_Late"))
# dput(x = DatesStrata2016.mat, file = "Objects/DatesStrata2016_Final.txt"); rm(DatesStrata2016.mat)
DatesStrata2016_Final <- dget(file = "Objects/DatesStrata2016_Final.txt")


## Sample sizes
# KMA2014_2016Strata_SampleSizes_Final <- KMA2014_2016Strata_SampleSizes[, "Final"]
# LateLateStrata <- grep(pattern = "LateLate", x = names(KMA2014_2016Strata_SampleSizes_Final))
# KMA2014_2016Strata_SampleSizes_Final[LateLateStrata-1] <- KMA2014_2016Strata_SampleSizes_Final[LateLateStrata-1] + KMA2014_2016Strata_SampleSizes_Final[LateLateStrata]
# KMA2014_2016Strata_SampleSizes_Final_Condense <- KMA2014_2016Strata_SampleSizes_Final[-LateLateStrata]
# dput(x = KMA2014_2016Strata_SampleSizes_Final_Condense, file = "Objects/KMA2014_2016Strata_SampleSizes_Final_Condense.txt")
KMA2014_2016Strata_SampleSizes_Final_Condense <- dget(file = "Objects/KMA2014_2016Strata_SampleSizes_Final_Condense.txt")


## Geographic headers
GeoHeader <- setNames(object = c(paste0("Alitak (257-10, 20, 50, 60, 70)"),
                                 paste0("Ayakulik-Halibut Bay (256-10", "\u2013", "256-30)"),
                                 paste0("Igvak (262-75, 80, 90, 95)"),
                                 paste0("Karluk-Sturgeon (255-10, 20; 256-40)"),
                                 paste0("Uganik-Kupreanof (253)"),
                                 paste0("Uyak (254)")),
                      nm = unlist(strsplit(x = KMA2016, split = "16")))
dput(x = GeoHeader, file = "Objects/GeoHeader.txt")
GeoHeader <- dget(file = "Objects/GeoHeader.txt")


# Get Final Estimates Objects
KMAfinalestimatesobjects <- list.files(path = "Estimates objects/Final", recursive = FALSE)
invisible(sapply(KMAfinalestimatesobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste("Estimates objects/Final", objct, sep = "/")), pos = 1) })); beep(2)
KMAfinalestimatesobjects; rm(KMAfinalestimatesobjects)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Defining caption variables

EstimatesStats <- c(KMA2014Strata_EstimatesStats, KMA2014_Annual_EstimatesStats,
                    KMA2015Strata_EstimatesStats, KMA2015_Annual_EstimatesStats)

HarvestEstimatesStats <- c(KMA2014Strata_HarvestEstimatesStats, KMA2014_Annual_HarvestEstimatesStats,
                           KMA2015Strata_HarvestEstimatesStats, KMA2015_Annual_HarvestEstimatesStats)

SheetNames <- sort(names(EstimatesStats))
names(SheetNames) <- SheetNames
mixvec <- SheetNames


# mix <- SheetNames[2]
harvest <- rbind(HarvestByStrata2014_Final, HarvestByStrata2015_Final)
dates <- rbind(DatesStrata2014_Final, DatesStrata2015_Final)
sampsize <- KMA2014_2015Strata_SampleSizes_Final_Condense



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

for(mix in SheetNames) {
  
  # String split by "_
  SheetNames.split <- unlist(strsplit(x = mix, split = "_"))
  
  
  # The first element is the geographic area + year
  geomix <- SheetNames.split[1]
  yr <- unlist(strsplit(x = geomix, split = "C"))[2]
  geo <- unlist(strsplit(x = geomix, split = "1"))[1]
  
  
  # If it is not an annual roll-up, then get the strata number + strata name
  if(length(SheetNames.split) > 1) {
    tempmix <- paste(c(SheetNames.split[2], SheetNames.split[3]), collapse = "_")
    
    Caption <- paste("Table X.-Estimates of stock composition (%) and stock-specific harvest for temporal stratum ",
                     SheetNames.split[2], " (", dates[geomix, tempmix],
                     "; Harvest=", formatC(x = harvest[geomix, tempmix], format = "f", digits = 0, big.mark = ","),
                     "; n=", sampsize[mix], ")", " of ", GeoHeader[geo], ", 20", yr,
                     ". Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.",
                     sep = '')
  } else {
    Caption <- paste("Table X.-Annual estimates of stock composition (%) and stock-specific harvest for ", GeoHeader[geo], ", 20", yr,
                     ". Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.",
                     sep = '')
  }
  
  Disclaimer <- "Note: Stock composition estimates may not sum to 100% and stock-specific harvest estimates may not sum to the total harvest due to rounding error."
  
  
  TableX <- matrix(data = "", nrow = 20, ncol = 13)
  
  TableX[1, 1] <- Caption
  TableX[2, c(2, 9)] <- c("Stock Composition", "Stock-specific Harvest")
  TableX[3, c(3, 10)] <- rep("90% CI", 2)
  TableX[4, c(1, 2:4, 6:7, 9:13, 5)] <- c("Reporting Group", rep(c("Median", "5%", "95%", "Mean", "SD"), 2), "P=0")
  TableX[5:18, 1] <- KMA14GroupsPC
  TableX[5:18, c(2:4, 6:7)] <- formatC(x = EstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
  TableX[5:18, 5] <- formatC(x = EstimatesStats[[mix]][, "P=0"], digits = 2, format = "f")
  TableX[5:18, 9:13] <- formatC(x = HarvestEstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")
  TableX[19, 11:12] <- c("Total", formatC(x = sum(HarvestEstimatesStats[[mix]][, "mean"]), digits = 0, format = "f", big.mark = ","))
  TableX[20, 1] <- Disclaimer
  
  
  write.xlsx(x = as.data.frame(TableX), 
             file = "Estimates tables/KMA Sockeye Estimates Tables.xlsx",
             col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = mix)
}; beep(5)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2014 Late Late Harvest for Appendix
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Defining caption variables

EstimatesStats <- dget(file = "Estimates objects/Round4Mixtures_2014_EstimatesStats.txt")
str(EstimatesStats)

HarvestEstimatesStats <- sapply(names(EstimatesStats), function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2014[strata.split[1], strata.split[2]],
        EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )



HarvestEstimatesStats <- c(KMA2014Strata_HarvestEstimatesStats, KMA2014_Annual_HarvestEstimatesStats,
                           KMA2015Strata_HarvestEstimatesStats, KMA2015_Annual_HarvestEstimatesStats)

SheetNames <- sort(names(EstimatesStats))
names(SheetNames) <- SheetNames
mixvec <- SheetNames


# mix <- SheetNames[2]
harvest <- HarvestByStrata2014
dates <- matrix(data = rep("August 25-29", 3), nrow = 3, ncol = 1, 
                dimnames = list(c("SKARLC14", "SUGANC14", "SUYAKC14"),
                                "4_LateLate"))
sampsize <- KMA2014_2015Strata_SampleSizes[, "Final"]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

for(mix in SheetNames) {
  
  # String split by "_
  SheetNames.split <- unlist(strsplit(x = mix, split = "_"))
  
  
  # The first element is the geographic area + year
  geomix <- SheetNames.split[1]
  yr <- unlist(strsplit(x = geomix, split = "C"))[2]
  geo <- unlist(strsplit(x = geomix, split = "1"))[1]
  
  
  # If it is not an annual roll-up, then get the strata number + strata name
  if(length(SheetNames.split) > 1) {
    tempmix <- paste(c(SheetNames.split[2], SheetNames.split[3]), collapse = "_")
    
    Caption <- paste("Table X.-Estimates of stock composition (%) and stock-specific harvest for temporal stratum ",
                     SheetNames.split[2], " (", dates[geomix, tempmix],
                     "; Harvest=", formatC(x = harvest[geomix, tempmix], format = "f", digits = 0, big.mark = ","),
                     "; n=", sampsize[mix], ")", " of the ", GeoHeader[geo], ", 20", yr,
                     ". Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.",
                     sep = '')
  } else {
    Caption <- paste("Table X.-Annual estimates of stock composition (%) and stock-specific harvest for the ", GeoHeader[geo], ", 20", yr,
                     ". Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and SD.",
                     sep = '')
  }
  
  Disclaimer <- "Note: Stock composition estimates may not sum to 100% and stock-specific harvest estimates may not sum to the total harvest due to rounding error."
  
  
  TableX <- matrix(data = "", nrow = 20, ncol = 13)
  
  TableX[1, 1] <- Caption
  TableX[2, c(2, 9)] <- c("Stock Composition", "Stock-specific Harvest")
  TableX[3, c(3, 10)] <- rep("90% CI", 2)
  TableX[4, c(1, 2:4, 6:7, 9:13, 5)] <- c("Reporting Group", rep(c("Median", "5%", "95%", "Mean", "SD"), 2), "P=0")
  TableX[5:18, 1] <- KMA14GroupsPC
  TableX[5:18, c(2:4, 6:7)] <- formatC(x = EstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
  TableX[5:18, 5] <- formatC(x = EstimatesStats[[mix]][, "P=0"], digits = 2, format = "f")
  TableX[5:18, 9:13] <- formatC(x = HarvestEstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")
  TableX[19, 11:12] <- c("Total", formatC(x = sum(HarvestEstimatesStats[[mix]][, "mean"]), digits = 0, format = "f", big.mark = ","))
  TableX[20, 1] <- Disclaimer
  
  
  write.xlsx(x = as.data.frame(TableX), 
             file = "Estimates tables/KMA Sockeye Estimates Tables.xlsx",
             col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = mix)
}; beep(5)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Annual KMA Percentages ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_Stratified_EstimatesStats.txt")
KMA2015_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_Stratified_EstimatesStats.txt")
KMA2016_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_Stratified_EstimatesStats.txt")



# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             3, 4), nrow = 2, ncol = 2, byrow = TRUE)

ProportionColors <- colorpanel(n = 3, low = "blue", high = "white")

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)



emf(file = "Figures/All Years/KMA Proportions 2014-2016.emf", width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 1), heights = c(2.7, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of KMA Harvest", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplot
par(mar = c(1, 1, 1.5, 1))
Barplot <- barplot2(height = t(cbind(KMA2014_Annual_Stratified_EstimatesStats[, "median"],
                                     KMA2015_Annual_Stratified_EstimatesStats[, "median"],
                                     KMA2016_Annual_Stratified_EstimatesStats[, "median"])) * 100, 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = t(cbind(KMA2014_Annual_Stratified_EstimatesStats[, "5%"],
                                   KMA2015_Annual_Stratified_EstimatesStats[, "5%"],
                                   KMA2016_Annual_Stratified_EstimatesStats[, "5%"])) * 100, 
                    ci.u = t(cbind(KMA2014_Annual_Stratified_EstimatesStats[, "95%"],
                                   KMA2015_Annual_Stratified_EstimatesStats[, "95%"],
                                   KMA2016_Annual_Stratified_EstimatesStats[, "95%"])) * 100, 
                    ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = 2014:2016, x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="")
abline(h = 0, xpd = FALSE)

mtext(text = Groups2Rows, side = 1, line = 0.66, at = apply(Barplot, 2, mean), adj = 0.5, cex = cex.xaxis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)


dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Annual KMA Harvest ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_Stratified_HarvestEstimatesStats.txt")
KMA2015_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_Stratified_HarvestEstimatesStats.txt")
KMA2016_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_Stratified_HarvestEstimatesStats.txt")



# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             3, 4), nrow = 2, ncol = 2, byrow = TRUE)

ProportionColors <- colorpanel(n = 3, low = "green", high = "white")

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5
ymax <- 700000  # max(sapply(list(KMA2014_Annual_Stratified_HarvestEstimatesStats, KMA2015_Annual_Stratified_HarvestEstimatesStats, KMA2016_Annual_Stratified_HarvestEstimatesStats), function(strata) {strata[, "95%"]}))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)



emf(file = "Figures/All Years/KMA Harvest 2014-2016.emf", width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 1), heights = c(2.7, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Number of Fish Harvested", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplot
par(mar = c(1, 1, 1.5, 1))
Barplot <- barplot2(height = t(cbind(KMA2014_Annual_Stratified_HarvestEstimatesStats[, "median"],
                                     KMA2015_Annual_Stratified_HarvestEstimatesStats[, "median"],
                                     KMA2016_Annual_Stratified_HarvestEstimatesStats[, "median"])), 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = t(cbind(KMA2014_Annual_Stratified_HarvestEstimatesStats[, "5%"],
                                   KMA2015_Annual_Stratified_HarvestEstimatesStats[, "5%"],
                                   KMA2016_Annual_Stratified_HarvestEstimatesStats[, "5%"])), 
                    ci.u = t(cbind(KMA2014_Annual_Stratified_HarvestEstimatesStats[, "95%"],
                                   KMA2015_Annual_Stratified_HarvestEstimatesStats[, "95%"],
                                   KMA2016_Annual_Stratified_HarvestEstimatesStats[, "95%"])), 
                    ylim = c(0, ymax), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 100000), labels = formatC(x = seq(0, ymax, 100000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
legend(legend = 2014:2016, x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="")
abline(h = 0, xpd = FALSE)

mtext(text = Groups2Rows, side = 1, line = 0.66, at = apply(Barplot, 2, mean), adj = 0.5, cex = cex.xaxis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Reporting Group", cex = cex.lab)


dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Re-Summarize All Estimates for Regional Rollups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA6GroupsPC <- c(KMA14GroupsPC[1], "Chignik", "Kodiak", KMA14GroupsPC[12:14])

KMA473PopsGroupVec6 <- KMA473PopsGroupVec14
KMA473PopsGroupVec6[KMA473PopsGroupVec6 %in% 2:3] <- 2
KMA473PopsGroupVec6[KMA473PopsGroupVec6 %in% 4:11] <- 3
KMA473PopsGroupVec6[KMA473PopsGroupVec6 %in% 12] <- 4
KMA473PopsGroupVec6[KMA473PopsGroupVec6 %in% 13] <- 5
KMA473PopsGroupVec6[KMA473PopsGroupVec6 %in% 14] <- 6

KMA473PopsGroupVec6_14RG <- c(1, rep(2, 2), rep(3, 8), 4, 5, 6)

dput(x = KMA6GroupsPC, file = "Objects/KMA6GroupsPC.txt")
dput(x = KMA473PopsGroupVec6, file = "Objects/KMA473PopsGroupVec6.txt")
dput(x = KMA473PopsGroupVec6_14RG, file = "Objects/KMA473PopsGroupVec6_14RG.txt")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2014
Round1Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")
sapply(Round1Mixtures_2014_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})
Round2Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2014_EstimatesStats.txt")
sapply(Round2Mixtures_2014_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})
Round3Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures_2014_EstimatesStats.txt")
sapply(Round3Mixtures_2014_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})
Round4Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round4Mixtures_2014_EstimatesStats.txt")
sapply(Round4Mixtures_2014_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})



KMA2014Strata_Regional_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, 
                               maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
                               mixvec = KMA2014Strata, prior = "",  
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(KMA2014Strata_Regional_Estimates)


# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(KMA2014Strata_Regional_Estimates, file = "Estimates objects/KMA2014Strata_Regional_Estimates.txt")
dput(KMA2014Strata_Regional_Estimates$Stats, file = "Estimates objects/KMA2014Strata_Regional_EstimatesStats.txt")

KMA2014Strata_Regional_EstimatesStats <- dget(file = "Estimates objects/KMA2014Strata_Regional_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(KMA2014Strata_Regional_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(KMA2014Strata_Regional_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUYAKC14_2_Middle still has SEAK issues
KMA2014Strata_Regional_EstimatesStats["SUYAKC14_2_Middle"]

which(sapply(KMA2014Strata_Regional_Estimates$Output, function(mix) {dim(mix)[1]}) == 200000)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2015
KMA2015Strata_Regional_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, 
                               maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
                               mixvec = KMA2015Strata, prior = "",  
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(KMA2015Strata_Regional_Estimates)


# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(KMA2015Strata_Regional_Estimates, file = "Estimates objects/KMA2015Strata_Regional_Estimates.txt")
dput(KMA2015Strata_Regional_Estimates$Stats, file = "Estimates objects/KMA2015Strata_Regional_EstimatesStats.txt")

KMA2015Strata_Regional_EstimatesStats <- dget(file = "Estimates objects/KMA2015Strata_Regional_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(KMA2015Strata_Regional_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(KMA2015Strata_Regional_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUGANC15_3_Late still has SEAK issues
KMA2015Strata_Regional_EstimatesStats["SUGANC15_3_Late"]

which(sapply(KMA2015Strata_Regional_Estimates$Output, function(mix) {dim(mix)[1]}) == 200000)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2016
Round1Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2016_EstimatesStats.txt")
sapply(Round1Mixtures_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})
Round2Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2016_EstimatesStats.txt")
sapply(Round2Mixtures_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})
Round3Mixtures_2016_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures_2016_EstimatesStats.txt")
sapply(Round3Mixtures_2016_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})



KMA2016Strata_Regional_Estimates <- 
  CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, 
                               maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
                               mixvec = KMA2016Strata, prior = "",  
                               ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(KMA2016Strata_Regional_Estimates)


# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(KMA2016Strata_Regional_Estimates, file = "Estimates objects/KMA2016Strata_Regional_Estimates.txt")
dput(KMA2016Strata_Regional_Estimates$Stats, file = "Estimates objects/KMA2016Strata_Regional_EstimatesStats.txt")

KMA2016Strata_Regional_EstimatesStats <- dget(file = "Estimates objects/KMA2016Strata_Regional_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(KMA2016Strata_Regional_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(KMA2016Strata_Regional_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  #  SKARLC16_1_Early still has PWS SEAK issues
KMA2016Strata_Regional_EstimatesStats["SKARLC16_1_Early"]

which(sapply(KMA2016Strata_Regional_Estimates$Output, function(mix) {dim(mix)[1]}) == 200000)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Regional Roll-Ups 2014 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
names(KMA2014Strata_Regional_EstimatesStats)

HarvestByStrata2014_Final

# Regional stratified estiamtes for 2014 Uganik, Uyak, and Karluk 3_Late (combining 3_Late with 4_LateLate)
sapply(rownames(HarvestByStrata2014)[!is.na(HarvestByStrata2014[, 4])], function(geomix) {
  assign(x = paste(geomix, "_3_Late_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = c(paste(geomix, "_3_Late", sep = ''), paste(geomix, "_4_LateLate", sep = '')), 
           catchvec = HarvestByStrata2014[geomix, c("3_Late", "4_LateLate")], newname = paste(geomix, "_3_Late_Stratified_Regional", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_3_Late_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_3_Late_Stratified_Regional.txt", sep = ''))
} ); beep(5)

str(SKARLC14_3_Late_Stratified)

# Update 2014 Strata results
KMA2014Strata_Regional_EstimatesStats[["SKARLC14_3_Late"]] <- SKARLC14_3_Late_Stratified$Stats
KMA2014Strata_Regional_EstimatesStats[["SUGANC14_3_Late"]] <- SUGANC14_3_Late_Stratified$Stats
KMA2014Strata_Regional_EstimatesStats[["SUYAKC14_3_Late"]] <- SUYAKC14_3_Late_Stratified$Stats

# Dput final Strata restuls
dput(x = KMA2014Strata_Regional_EstimatesStats, file = "Estimates objects/Final/KMA2014Strata_Regional_EstimatesStats.txt")
dput(x = KMA2015Strata_Regional_EstimatesStats, file = "Estimates objects/Final/KMA2015Strata_Regional_EstimatesStats.txt")
dput(x = KMA2016Strata_Regional_EstimatesStats, file = "Estimates objects/Final/KMA2016Strata_Regional_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(KMA2014, function(geomix) {
  assign(x = paste(geomix, "_Annual_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = paste(geomix, colnames(HarvestByStrata2014)[!is.na(HarvestByStrata2014[geomix,])], sep = "_"), 
           catchvec = HarvestByStrata2014[geomix, !is.na(HarvestByStrata2014[geomix,])], 
           newname = paste(geomix, "_Annual_Stratified_Regional", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_Annual_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_Annual_Stratified_Regional.txt", sep = ''))
} ); beep(5)

str(SALITC14_Annual_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2014_Annual_Regional_EstimatesStats <- sapply(KMA2014, function(geomix) {
  Stats <- get(paste(geomix, "_Annual_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2014_Annual_Regional_EstimatesStats)
dput(x = KMA2014_Annual_Regional_EstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_Regional_EstimatesStats.txt")


# Create a list object with all Stratified Annual Harvest
KMA2014_Annual_Regional_HarvestEstimatesStats <- sapply(KMA2014, function(strata) {
  cbind(KMA2014_Annual_Regional_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2014_Final[strata, ], na.rm = TRUE),
        KMA2014_Annual_Regional_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2014_Annual_Regional_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_Regional_HarvestEstimatesStats.txt")
str(KMA2014_Annual_Regional_HarvestEstimatesStats)




# Create Strata_Regional_HarvestEstimatesStats
KMA2014Strata_Regional_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_Regional_EstimatesStats.txt")
KMA2014Strata_Regional_HarvestEstimatesStats <- sapply(names(KMA2014Strata_Regional_EstimatesStats)[-grep(pattern = "LateLate", x = names(KMA2014Strata_Regional_EstimatesStats))], function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(KMA2014Strata_Regional_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2014_Final[strata.split[1], strata.split[2]],
        KMA2014Strata_Regional_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )
str(KMA2014Strata_Regional_HarvestEstimatesStats)

dput(x = KMA2014Strata_Regional_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014Strata_Regional_HarvestEstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create Temporal Stratified rollup
# 14RG
#~~~~~~~~~~~~~~~~~~
# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(c("1_Early", "2_Middle"), function(tempmix) {
  assign(x = paste0("KMA2014_", tempmix, "_Temporal_Stratified"), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC2, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = paste(rownames(HarvestByStrata2014)[!is.na(HarvestByStrata2014[,tempmix])], tempmix, sep = "_"), 
           catchvec = HarvestByStrata2014[!is.na(HarvestByStrata2014[, tempmix]), tempmix], 
           newname = paste0("KMA2014_", tempmix, "_Temporal_Stratified"), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste0("KMA2014_", tempmix, "_Temporal_Stratified")), file = paste("Estimates objects/KMA2014_", tempmix, "_Temporal_Stratified.txt", sep = ''))
} ); beep(5)

KMA2014_3_Late_Temporal_Stratified <- StratifiedEstimator.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC2, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
  mixvec = c(paste(rownames(HarvestByStrata2014)[!is.na(HarvestByStrata2014[, "3_Late"])], "3_Late", sep = "_"), 
             paste(rownames(HarvestByStrata2014)[!is.na(HarvestByStrata2014[, "4_LateLate"])], "4_LateLate", sep = "_")), 
  catchvec = c(HarvestByStrata2014[!is.na(HarvestByStrata2014[, "3_Late"]), "3_Late"], 
               HarvestByStrata2014[!is.na(HarvestByStrata2014[, "4_LateLate"]), "4_LateLate"]), 
  newname = "KMA2014_3_Late_Temporal_Stratified", priorname = '',
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1)
dput(x = KMA2014_3_Late_Temporal_Stratified, file = paste("Estimates objects/KMA2014_3_Late_Temporal_Stratified.txt", sep = ''))


# Create a list object with all Stratified Temporal Rollups for the "Estimates objects/Final" folder
KMA2014_Temporal_Stratified_EstimatesStats <- sapply(c("1_Early", "2_Middle", "3_Late"), function(tempmix) {
  Stats <- get(paste0("KMA2014_", tempmix, "_Temporal_Stratified"))$Stats
  Stats
}, simplify = FALSE)
str(KMA2014_Temporal_Stratified_EstimatesStats)
dput(x = KMA2014_Temporal_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2014_Temporal_Stratified_EstimatesStats.txt")


# Create a list object with all Stratified Annual Harvest
KMA2014_Temporal_Stratified_HarvestEstimatesStats <- sapply(c("1_Early", "2_Middle", "3_Late"), function(tempmix) {
  cbind(KMA2014_Temporal_Stratified_EstimatesStats[[tempmix]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2014_Final[, tempmix], na.rm = TRUE),
        KMA2014_Temporal_Stratified_EstimatesStats[[tempmix]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2014_Temporal_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014_Temporal_Stratified_HarvestEstimatesStats.txt")
str(KMA2014_Temporal_Stratified_HarvestEstimatesStats)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create KMA-wide annual rollup
# 14RG
#~~~~~~~~~~~~~~~~~~
KMA2014_Annual_Stratified <- StratifiedEstimator.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
  mixvec = KMA2014Strata, 
  catchvec = as.vector(na.omit(as.vector(t(HarvestByStrata2014[sort(KMA2014), ])))), 
  newname = "KMA2014_Annual_Stratified", priorname = '',
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1
)
str(KMA2014_Annual_Stratified)
dput(x = KMA2014_Annual_Stratified, file = "Estimates objects/KMA2014_Annual_Stratified.txt")
KMA2014_Annual_Stratified_EstimatesStats <- KMA2014_Annual_Stratified$Stats
dput(x = KMA2014_Annual_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_Stratified_EstimatesStats.txt")

# Create a list object with all Stratified Annual Harvest
#~~~~~~~~~~~~~~~~~~
KMA2014_Annual_Stratified_HarvestEstimatesStats <-
  cbind(KMA2014_Annual_Stratified_EstimatesStats[, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2014_Final, na.rm = TRUE),
        KMA2014_Annual_Stratified_EstimatesStats[, c("P=0", "GR")])
dput(x = KMA2014_Annual_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_Stratified_HarvestEstimatesStats.txt")
str(KMA2014_Annual_Stratified_HarvestEstimatesStats)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6RG
#~~~~~~~~~~~~~~~~~~
KMA2014_Annual_Regional_Stratified <- StratifiedEstimator.GCL(
  groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
  mixvec = KMA2014Strata, 
  catchvec = as.vector(na.omit(as.vector(t(HarvestByStrata2014[sort(KMA2014), ])))), 
  newname = "KMA2014_Annual_Regional_Stratified", priorname = '',
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1
)
str(KMA2014_Annual_Regional_Stratified)
dput(x = KMA2014_Annual_Regional_Stratified, file = "Estimates objects/KMA2014_Annual_Regional_Stratified.txt")
KMA2014_Annual_Regional_Stratified_EstimatesStats <- KMA2014_Annual_Regional_Stratified$Stats
dput(x = KMA2014_Annual_Regional_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_Regional_Stratified_EstimatesStats.txt")

# Create a list object with all Stratified Annual_Regional Harvest
#~~~~~~~~~~~~~~~~~~
KMA2014_Annual_Regional_Stratified_HarvestEstimatesStats <-
  cbind(KMA2014_Annual_Regional_Stratified_EstimatesStats[, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2014_Final, na.rm = TRUE),
        KMA2014_Annual_Regional_Stratified_EstimatesStats[, c("P=0", "GR")])
dput(x = KMA2014_Annual_Regional_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014_Annual_Regional_Stratified_HarvestEstimatesStats.txt")
str(KMA2014_Annual_Regional_Stratified_HarvestEstimatesStats)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create KMA-wide temporal rollup
# 14RG
#~~~~~~~~~~~~~~~~~~
sapply(c("Early", "Middle", "Late"), function(tempmix) {
  assign(x = paste("KMA2014_", tempmix, "_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = grep(pattern = tempmix, x = KMA2014Strata, value = TRUE), 
           catchvec = t(HarvestByStrata2014)[!is.na(t(HarvestByStrata2014))][grep(pattern = tempmix, x = KMA2014Strata)], 
           newname = paste("KMA2014_", tempmix, "_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste("KMA2014_", tempmix, "_Stratified", sep = '')), file = paste("Estimates objects/KMA2014_", tempmix, "_Stratified.txt", sep = ''))
} ); beep(5)
str(KMA2014_Early_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2014_Temporal_Annual_EstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  Stats <- get(paste("KMA2014_", tempmix, "_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2014_Temporal_Annual_EstimatesStats)
dput(x = KMA2014_Temporal_Annual_EstimatesStats, file = "Estimates objects/Final/KMA2014_Temporal_Annual_EstimatesStats.txt")


# Create a list object with all Stratified Annual Harvest
KMA2014_Temporal_Annual_HarvestEstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  cbind(KMA2014_Temporal_Annual_EstimatesStats[[tempmix]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2014_Final[, grep(x = colnames(HarvestByStrata2014_Final), pattern = tempmix)], na.rm = TRUE),
        KMA2014_Temporal_Annual_EstimatesStats[[tempmix]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2014_Temporal_Annual_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014_Temporal_Annual_HarvestEstimatesStats.txt")
str(KMA2014_Temporal_Annual_HarvestEstimatesStats)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6RG
#~~~~~~~~~~~~~~~~~~
sapply(c("Early", "Middle", "Late"), function(tempmix) {
  assign(x = paste("KMA2014_", tempmix, "_Regional_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = grep(pattern = tempmix, x = KMA2014Strata, value = TRUE), 
           catchvec = t(HarvestByStrata2014)[!is.na(t(HarvestByStrata2014))][grep(pattern = tempmix, x = KMA2014Strata)], 
           newname = paste("KMA2014_", tempmix, "_Regional_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste("KMA2014_", tempmix, "_Regional_Stratified", sep = '')), file = paste("Estimates objects/KMA2014_", tempmix, "_Regional_Stratified.txt", sep = ''))
} ); beep(5)
str(KMA2014_Early_Regional_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2014_Temporal_Annual_Regional_EstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  Stats <- get(paste("KMA2014_", tempmix, "_Regional_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2014_Temporal_Annual_Regional_EstimatesStats)
dput(x = KMA2014_Temporal_Annual_Regional_EstimatesStats, file = "Estimates objects/Final/KMA2014_Temporal_Annual_Regional_EstimatesStats.txt")


# Create a list object with all Stratified Annual Harvest
KMA2014_Temporal_Annual_Regional_HarvestEstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  cbind(KMA2014_Temporal_Annual_Regional_EstimatesStats[[tempmix]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2014_Final[, grep(x = colnames(HarvestByStrata2014_Final), pattern = tempmix)], na.rm = TRUE),
        KMA2014_Temporal_Annual_Regional_EstimatesStats[[tempmix]][, c("P=0", "GR")])
}, simplify = FALSE )
str(KMA2014_Temporal_Annual_Regional_HarvestEstimatesStats)
dput(x = KMA2014_Temporal_Annual_Regional_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2014_Temporal_Annual_Regional_HarvestEstimatesStats.txt")






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Regional Roll-Ups 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HarvestByStrata2015_Final

KMA2015Strata_Regional_Estimates <- dget(file = "Estimates objects/KMA2015Strata_Regional_Estimates.txt")
which(sapply(KMA2015Strata_Regional_Estimates$Output, function(mix) {dim(mix)[1]}) == 200000)


# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(KMA2015[-3], function(geomix) {
  assign(x = paste(geomix, "_Annual_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = paste(geomix, colnames(HarvestByStrata2015)[!is.na(HarvestByStrata2015[geomix,])], sep = "_"), 
           catchvec = HarvestByStrata2015[geomix, !is.na(HarvestByStrata2015[geomix,])], 
           newname = paste(geomix, "_Annual_Stratified_Regional", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_Annual_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_Annual_Stratified_Regional.txt", sep = ''))
} ); beep(5)

str(SALITC15_Annual_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2015_Annual_Regional_EstimatesStats <- sapply(KMA2015[-3], function(geomix) {
  Stats <- get(paste(geomix, "_Annual_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2015_Annual_Regional_EstimatesStats)

# Add Igvak 2_Middle 2015
SIGVAC15_2_Middle_EstimatesStats <- dget(file = "Estimates objects/KMA2015Strata_Regional_EstimatesStats.txt")["SIGVAC15_2_Middle"]
str(SIGVAC15_2_Middle_EstimatesStats)

KMA2015_Annual_Regional_EstimatesStats <- c(KMA2015_Annual_Regional_EstimatesStats[1:2], list("SIGVAC15" = SIGVAC15_2_Middle_EstimatesStats[[1]]), KMA2015_Annual_Regional_EstimatesStats[3:5])
str(KMA2015_Annual_Regional_EstimatesStats)
dput(x = KMA2015_Annual_Regional_EstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_Regional_EstimatesStats.txt")




# Create a list object with all Stratified Annual Harvest
KMA2015_Annual_Regional_HarvestEstimatesStats <- sapply(KMA2015, function(strata) {
  cbind(KMA2015_Annual_Regional_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2015_Final[strata, ], na.rm = TRUE),
        KMA2015_Annual_Regional_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2015_Annual_Regional_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_Regional_HarvestEstimatesStats.txt")
str(KMA2015_Annual_Regional_HarvestEstimatesStats)




# Create Strata_Regional_HarvestEstimatesStats
KMA2015Strata_Regional_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_Regional_EstimatesStats.txt")
KMA2015Strata_Regional_HarvestEstimatesStats <- sapply(names(KMA2015Strata_Regional_EstimatesStats), function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(KMA2015Strata_Regional_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2015_Final[strata.split[1], strata.split[2]],
        KMA2015Strata_Regional_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )
str(KMA2015Strata_Regional_HarvestEstimatesStats)

dput(x = KMA2015Strata_Regional_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015Strata_Regional_HarvestEstimatesStats.txt")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create Temporal Stratified rollup
# 14RG
#~~~~~~~~~~~~~~~~~~
# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(c("1_Early", "2_Middle", "3_Late"), function(tempmix) {
  assign(x = paste0("KMA2015_", tempmix, "_Temporal_Stratified"), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC2, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = paste(rownames(HarvestByStrata2015)[!is.na(HarvestByStrata2015[,tempmix])], tempmix, sep = "_"), 
           catchvec = HarvestByStrata2015[!is.na(HarvestByStrata2015[, tempmix]), tempmix], 
           newname = paste0("KMA2015_", tempmix, "_Temporal_Stratified"), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste0("KMA2015_", tempmix, "_Temporal_Stratified")), file = paste("Estimates objects/KMA2015_", tempmix, "_Temporal_Stratified.txt", sep = ''))
} ); beep(5)


# Create a list object with all Stratified Temporal Rollups for the "Estimates objects/Final" folder
KMA2015_Temporal_Stratified_EstimatesStats <- sapply(c("1_Early", "2_Middle", "3_Late"), function(tempmix) {
  Stats <- get(paste0("KMA2015_", tempmix, "_Temporal_Stratified"))$Stats
  Stats
}, simplify = FALSE)
str(KMA2015_Temporal_Stratified_EstimatesStats)
dput(x = KMA2015_Temporal_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2015_Temporal_Stratified_EstimatesStats.txt")


# Create a list object with all Stratified Annual Harvest
KMA2015_Temporal_Stratified_HarvestEstimatesStats <- sapply(c("1_Early", "2_Middle", "3_Late"), function(tempmix) {
  cbind(KMA2015_Temporal_Stratified_EstimatesStats[[tempmix]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2015_Final[, tempmix], na.rm = TRUE),
        KMA2015_Temporal_Stratified_EstimatesStats[[tempmix]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2015_Temporal_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015_Temporal_Stratified_HarvestEstimatesStats.txt")
str(KMA2015_Temporal_Stratified_HarvestEstimatesStats)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create KMA-wide annual rollup
# 14RG
#~~~~~~~~~~~~~~~~~~
KMA2015_Annual_Stratified <- StratifiedEstimator.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
  mixvec = KMA2015Strata, 
  catchvec = as.vector(na.omit(as.vector(t(HarvestByStrata2015_Final[sort(KMA2015), ])))), 
  newname = "KMA2015_Annual_Stratified", priorname = '',
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1
)
str(KMA2015_Annual_Stratified)
dput(x = KMA2015_Annual_Stratified, file = "Estimates objects/KMA2015_Annual_Stratified.txt")
KMA2015_Annual_Stratified_EstimatesStats <- KMA2015_Annual_Stratified$Stats
dput(x = KMA2015_Annual_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_Stratified_EstimatesStats.txt")

# Create a list object with all Stratified Annual Harvest
#~~~~~~~~~~~~~~~~~~
KMA2015_Annual_Stratified_HarvestEstimatesStats <-
  cbind(KMA2015_Annual_Stratified_EstimatesStats[, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2015_Final, na.rm = TRUE),
        KMA2015_Annual_Stratified_EstimatesStats[, c("P=0", "GR")])
dput(x = KMA2015_Annual_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_Stratified_HarvestEstimatesStats.txt")
str(KMA2015_Annual_Stratified_HarvestEstimatesStats)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6RG
#~~~~~~~~~~~~~~~~~~
KMA2015_Annual_Regional_Stratified <- StratifiedEstimator.GCL(
  groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
  mixvec = KMA2015Strata, 
  catchvec = as.vector(na.omit(as.vector(t(HarvestByStrata2015_Final[sort(KMA2015), ])))), 
  newname = "KMA2015_Annual_Regional_Stratified", priorname = '',
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1
)
str(KMA2015_Annual_Regional_Stratified)
dput(x = KMA2015_Annual_Regional_Stratified, file = "Estimates objects/KMA2015_Annual_Regional_Stratified.txt")
KMA2015_Annual_Regional_Stratified_EstimatesStats <- KMA2015_Annual_Regional_Stratified$Stats
dput(x = KMA2015_Annual_Regional_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_Regional_Stratified_EstimatesStats.txt")

# Create a list object with all Stratified Annual_Regional Harvest
#~~~~~~~~~~~~~~~~~~
KMA2015_Annual_Regional_Stratified_HarvestEstimatesStats <-
  cbind(KMA2015_Annual_Regional_Stratified_EstimatesStats[, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2015_Final, na.rm = TRUE),
        KMA2015_Annual_Regional_Stratified_EstimatesStats[, c("P=0", "GR")])
dput(x = KMA2015_Annual_Regional_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_Regional_Stratified_HarvestEstimatesStats.txt")
str(KMA2015_Annual_Regional_Stratified_HarvestEstimatesStats)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create KMA-wide temporal rollup
# 14RG
#~~~~~~~~~~~~~~~~~~
sapply(c("Early", "Middle", "Late"), function(tempmix) {
  assign(x = paste("KMA2015_", tempmix, "_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = grep(pattern = tempmix, x = KMA2015Strata, value = TRUE), 
           catchvec = t(HarvestByStrata2015)[!is.na(t(HarvestByStrata2015))][grep(pattern = tempmix, x = KMA2015Strata)], 
           newname = paste("KMA2015_", tempmix, "_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste("KMA2015_", tempmix, "_Stratified", sep = '')), file = paste("Estimates objects/KMA2015_", tempmix, "_Stratified.txt", sep = ''))
} ); beep(5)
str(KMA2015_Early_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2015_Temporal_Annual_EstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  Stats <- get(paste("KMA2015_", tempmix, "_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2015_Temporal_Annual_EstimatesStats)
dput(x = KMA2015_Temporal_Annual_EstimatesStats, file = "Estimates objects/Final/KMA2015_Temporal_Annual_EstimatesStats.txt")


# Create a list object with all Stratified Annual Harvest
KMA2015_Temporal_Annual_HarvestEstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  cbind(KMA2015_Temporal_Annual_EstimatesStats[[tempmix]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2015_Final[, grep(x = colnames(HarvestByStrata2015_Final), pattern = tempmix)], na.rm = TRUE),
        KMA2015_Temporal_Annual_EstimatesStats[[tempmix]][, c("P=0", "GR")])
}, simplify = FALSE )
str(KMA2015_Temporal_Annual_HarvestEstimatesStats)
dput(x = KMA2015_Temporal_Annual_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015_Temporal_Annual_HarvestEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6RG
#~~~~~~~~~~~~~~~~~~
sapply(c("Early", "Middle", "Late"), function(tempmix) {
  assign(x = paste("KMA2015_", tempmix, "_Regional_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = grep(pattern = tempmix, x = KMA2015Strata, value = TRUE), 
           catchvec = t(HarvestByStrata2015)[!is.na(t(HarvestByStrata2015))][grep(pattern = tempmix, x = KMA2015Strata)], 
           newname = paste("KMA2015_", tempmix, "_Regional_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste("KMA2015_", tempmix, "_Regional_Stratified", sep = '')), file = paste("Estimates objects/KMA2015_", tempmix, "_Regional_Stratified.txt", sep = ''))
} ); beep(5)
str(KMA2015_Early_Regional_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2015_Temporal_Annual_Regional_EstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  Stats <- get(paste("KMA2015_", tempmix, "_Regional_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2015_Temporal_Annual_Regional_EstimatesStats)
dput(x = KMA2015_Temporal_Annual_Regional_EstimatesStats, file = "Estimates objects/Final/KMA2015_Temporal_Annual_Regional_EstimatesStats.txt")


# Create a list object with all Stratified Annual Harvest
KMA2015_Temporal_Annual_Regional_HarvestEstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  cbind(KMA2015_Temporal_Annual_Regional_EstimatesStats[[tempmix]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2015_Final[, grep(x = colnames(HarvestByStrata2015_Final), pattern = tempmix)], na.rm = TRUE),
        KMA2015_Temporal_Annual_Regional_EstimatesStats[[tempmix]][, c("P=0", "GR")])
}, simplify = FALSE )
str(KMA2015_Temporal_Annual_Regional_HarvestEstimatesStats)
dput(x = KMA2015_Temporal_Annual_Regional_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015_Temporal_Annual_Regional_HarvestEstimatesStats.txt")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stratified Regional Roll-Ups 2016 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HarvestByStrata2016_Final

KMA2016Strata_Regional_Estimates <- dget(file = "Estimates objects/KMA2016Strata_Regional_Estimates.txt")
which(sapply(KMA2016Strata_Regional_Estimates$Output, function(mix) {dim(mix)[1]}) == 200000)


# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(KMA2016, function(geomix) {
  assign(x = paste(geomix, "_Annual_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = paste(geomix, colnames(HarvestByStrata2016)[!is.na(HarvestByStrata2016[geomix,])], sep = "_"), 
           catchvec = HarvestByStrata2016[geomix, !is.na(HarvestByStrata2016[geomix,])], 
           newname = paste(geomix, "_Annual_Stratified_Regional", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_Annual_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_Annual_Stratified_Regional.txt", sep = ''))
} ); beep(5)

str(SALITC16_Annual_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2016_Annual_Regional_EstimatesStats <- sapply(KMA2016, function(geomix) {
  Stats <- get(paste(geomix, "_Annual_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2016_Annual_Regional_EstimatesStats)
dput(x = KMA2016_Annual_Regional_EstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_Regional_EstimatesStats.txt")




# Create a list object with all Stratified Annual Harvest
KMA2016_Annual_Regional_HarvestEstimatesStats <- sapply(KMA2016, function(strata) {
  cbind(KMA2016_Annual_Regional_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2016_Final[strata, ], na.rm = TRUE),
        KMA2016_Annual_Regional_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2016_Annual_Regional_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_Regional_HarvestEstimatesStats.txt")
str(KMA2016_Annual_Regional_HarvestEstimatesStats)




# Create Strata_Regional_HarvestEstimatesStats
KMA2016Strata_Regional_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_Regional_EstimatesStats.txt")
KMA2016Strata_Regional_HarvestEstimatesStats <- sapply(names(KMA2016Strata_Regional_EstimatesStats), function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(KMA2016Strata_Regional_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2016_Final[strata.split[1], strata.split[2]],
        KMA2016Strata_Regional_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )
str(KMA2016Strata_Regional_HarvestEstimatesStats)

dput(x = KMA2016Strata_Regional_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016Strata_Regional_HarvestEstimatesStats.txt")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create Temporal Stratified rollup
# 14RG
#~~~~~~~~~~~~~~~~~~
# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(c("1_Early", "2_Middle", "3_Late"), function(tempmix) {
  assign(x = paste0("KMA2016_", tempmix, "_Temporal_Stratified"), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC2, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = paste(rownames(HarvestByStrata2016)[!is.na(HarvestByStrata2016[,tempmix])], tempmix, sep = "_"), 
           catchvec = HarvestByStrata2016[!is.na(HarvestByStrata2016[, tempmix]), tempmix], 
           newname = paste0("KMA2016_", tempmix, "_Temporal_Stratified"), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste0("KMA2016_", tempmix, "_Temporal_Stratified")), file = paste("Estimates objects/KMA2016_", tempmix, "_Temporal_Stratified.txt", sep = ''))
} ); beep(5)


# Create a list object with all Stratified Temporal Rollups for the "Estimates objects/Final" folder
KMA2016_Temporal_Stratified_EstimatesStats <- sapply(c("1_Early", "2_Middle", "3_Late"), function(tempmix) {
  Stats <- get(paste0("KMA2016_", tempmix, "_Temporal_Stratified"))$Stats
  Stats
}, simplify = FALSE)
str(KMA2016_Temporal_Stratified_EstimatesStats)
dput(x = KMA2016_Temporal_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2016_Temporal_Stratified_EstimatesStats.txt")


# Create a list object with all Stratified Annual Harvest
KMA2016_Temporal_Stratified_HarvestEstimatesStats <- sapply(c("1_Early", "2_Middle", "3_Late"), function(tempmix) {
  cbind(KMA2016_Temporal_Stratified_EstimatesStats[[tempmix]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2016_Final[, tempmix], na.rm = TRUE),
        KMA2016_Temporal_Stratified_EstimatesStats[[tempmix]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2016_Temporal_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016_Temporal_Stratified_HarvestEstimatesStats.txt")
str(KMA2016_Temporal_Stratified_HarvestEstimatesStats)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create KMA-wide annual rollup
# 14RG
#~~~~~~~~~~~~~~~~~~
KMA2016_Annual_Stratified <- StratifiedEstimator.GCL(
  groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
  mixvec = KMA2016Strata, 
  catchvec = as.vector(na.omit(as.vector(t(HarvestByStrata2016_Final[sort(KMA2016), ])))), 
  newname = "KMA2016_Annual_Stratified", priorname = '',
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1
)
str(KMA2016_Annual_Stratified)
dput(x = KMA2016_Annual_Stratified, file = "Estimates objects/KMA2016_Annual_Stratified.txt")
KMA2016_Annual_Stratified_EstimatesStats <- KMA2016_Annual_Stratified$Stats
dput(x = KMA2016_Annual_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_Stratified_EstimatesStats.txt")

# Create a list object with all Stratified Annual Harvest
#~~~~~~~~~~~~~~~~~~
KMA2016_Annual_Stratified_HarvestEstimatesStats <-
  cbind(KMA2016_Annual_Stratified_EstimatesStats[, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2016_Final, na.rm = TRUE),
        KMA2016_Annual_Stratified_EstimatesStats[, c("P=0", "GR")])
dput(x = KMA2016_Annual_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_Stratified_HarvestEstimatesStats.txt")
str(KMA2016_Annual_Stratified_HarvestEstimatesStats)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6RG
#~~~~~~~~~~~~~~~~~~
KMA2016_Annual_Regional_Stratified <- StratifiedEstimator.GCL(
  groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
  mixvec = KMA2016Strata, 
  catchvec = as.vector(na.omit(as.vector(t(HarvestByStrata2016_Final[sort(KMA2016), ])))), 
  newname = "KMA2016_Annual_Regional_Stratified", priorname = '',
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1
)
str(KMA2016_Annual_Regional_Stratified)
dput(x = KMA2016_Annual_Regional_Stratified, file = "Estimates objects/KMA2016_Annual_Regional_Stratified.txt")
KMA2016_Annual_Regional_Stratified_EstimatesStats <- KMA2016_Annual_Regional_Stratified$Stats
dput(x = KMA2016_Annual_Regional_Stratified_EstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_Regional_Stratified_EstimatesStats.txt")

# Create a list object with all Stratified Annual_Regional Harvest
#~~~~~~~~~~~~~~~~~~
KMA2016_Annual_Regional_Stratified_HarvestEstimatesStats <-
  cbind(KMA2016_Annual_Regional_Stratified_EstimatesStats[, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2016_Final, na.rm = TRUE),
        KMA2016_Annual_Regional_Stratified_EstimatesStats[, c("P=0", "GR")])
dput(x = KMA2016_Annual_Regional_Stratified_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016_Annual_Regional_Stratified_HarvestEstimatesStats.txt")
str(KMA2016_Annual_Regional_Stratified_HarvestEstimatesStats)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create KMA-wide temporal rollup
# 14RG
#~~~~~~~~~~~~~~~~~~
sapply(c("Early", "Middle", "Late"), function(tempmix) {
  assign(x = paste("KMA2016_", tempmix, "_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA14GroupsPC), groupnames = KMA14GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = grep(pattern = tempmix, x = KMA2016Strata, value = TRUE), 
           catchvec = t(HarvestByStrata2016)[!is.na(t(HarvestByStrata2016))][grep(pattern = tempmix, x = KMA2016Strata)], 
           newname = paste("KMA2016_", tempmix, "_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste("KMA2016_", tempmix, "_Stratified", sep = '')), file = paste("Estimates objects/KMA2016_", tempmix, "_Stratified.txt", sep = ''))
} ); beep(5)
str(KMA2016_Early_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2016_Temporal_Annual_EstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  Stats <- get(paste("KMA2016_", tempmix, "_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2016_Temporal_Annual_EstimatesStats)
dput(x = KMA2016_Temporal_Annual_EstimatesStats, file = "Estimates objects/Final/KMA2016_Temporal_Annual_EstimatesStats.txt")


# Create a list object with all Stratified Annual Harvest
KMA2016_Temporal_Annual_HarvestEstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  cbind(KMA2016_Temporal_Annual_EstimatesStats[[tempmix]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2016_Final[, grep(x = colnames(HarvestByStrata2016_Final), pattern = tempmix)], na.rm = TRUE),
        KMA2016_Temporal_Annual_EstimatesStats[[tempmix]][, c("P=0", "GR")])
}, simplify = FALSE )
str(KMA2016_Temporal_Annual_HarvestEstimatesStats)
dput(x = KMA2016_Temporal_Annual_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016_Temporal_Annual_HarvestEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6RG
#~~~~~~~~~~~~~~~~~~
sapply(c("Early", "Middle", "Late"), function(tempmix) {
  assign(x = paste("KMA2016_", tempmix, "_Regional_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = KMA473PopsGroupVec6_14RG, groupnames = KMA6GroupsPC, maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
           mixvec = grep(pattern = tempmix, x = KMA2016Strata, value = TRUE), 
           catchvec = t(HarvestByStrata2016)[!is.na(t(HarvestByStrata2016))][grep(pattern = tempmix, x = KMA2016Strata)], 
           newname = paste("KMA2016_", tempmix, "_Regional_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste("KMA2016_", tempmix, "_Regional_Stratified", sep = '')), file = paste("Estimates objects/KMA2016_", tempmix, "_Regional_Stratified.txt", sep = ''))
} ); beep(5)
str(KMA2016_Early_Regional_Stratified)


# Create a list object with all Stratified Annual Rollups for the "Estimates objects/Final" folder
KMA2016_Temporal_Annual_Regional_EstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  Stats <- get(paste("KMA2016_", tempmix, "_Regional_Stratified", sep = ''))$Stats
  Stats
}, simplify = FALSE)
str(KMA2016_Temporal_Annual_Regional_EstimatesStats)
dput(x = KMA2016_Temporal_Annual_Regional_EstimatesStats, file = "Estimates objects/Final/KMA2016_Temporal_Annual_Regional_EstimatesStats.txt")


# Create a list object with all Stratified Annual Harvest
KMA2016_Temporal_Annual_Regional_HarvestEstimatesStats <- sapply(c("Early", "Middle", "Late"), function(tempmix) {
  cbind(KMA2016_Temporal_Annual_Regional_EstimatesStats[[tempmix]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2016_Final[, grep(x = colnames(HarvestByStrata2016_Final), pattern = tempmix)], na.rm = TRUE),
        KMA2016_Temporal_Annual_Regional_EstimatesStats[[tempmix]][, c("P=0", "GR")])
}, simplify = FALSE )
str(KMA2016_Temporal_Annual_Regional_HarvestEstimatesStats)
dput(x = KMA2016_Temporal_Annual_Regional_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2016_Temporal_Annual_Regional_HarvestEstimatesStats.txt")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Percentages for KMA Strata Mixtures 14RG/6RG ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")
KMA2015Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")
KMA2016Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_EstimatesStats.txt")
str(KMA2014Strata_EstimatesStats)

KMA2014Strata_Regional_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_Regional_EstimatesStats.txt")
KMA2015Strata_Regional_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_Regional_EstimatesStats.txt")
KMA2016Strata_Regional_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_Regional_EstimatesStats.txt")
str(KMA2014Strata_Regional_EstimatesStats)


sapply(c(KMA2014Strata_EstimatesStats, KMA2015Strata_EstimatesStats, KMA2016Strata_EstimatesStats), function(mix) {table(mix[, "GR"] > 1.2)})

round(KMA2014Strata_EstimatesStats$SUYAKC14_2_Middle[, c("mean", "GR")], 3)
round(KMA2015Strata_EstimatesStats$SUGANC15_3_Late[, c("mean", "GR")], 3)
round(KMA2016Strata_EstimatesStats$SKARLC16_1_Early[, c("mean", "GR")], 3)

sapply(c(KMA2014Strata_Regional_EstimatesStats, KMA2015Strata_Regional_EstimatesStats, KMA2016Strata_Regional_EstimatesStats), function(mix) {table(mix[, "GR"] > 1.2)})





# Three barplot layout
layoutmat <- matrix(data=c(  1, 2, 3,
                             1, 4, 5,
                             1, 6, 7,
                             8, 9, 10), nrow = 4, ncol = 3, byrow = TRUE)



GeoMix <- c("SUGANC", "SUYAKC", "SKARLC", "SAYAKC", "SALITC", "SIGVAC")

filenames <- setNames(object = c("Uganik Proportions 2014-2016", 
                                 "Uyak Proportions 2014-2016",
                                 "Karluk Proportions 2014-2016",
                                 "Ayakulik Proportions 2014-2016",
                                 "Alitak Proportions 2014-2016",
                                 "Igvak Proportions 2014-2016"), nm = GeoMix)

# If showing proportions (percetages) use blue, otherwise green as "low"
ProportionColors <- colorpanel(n = 3, low = "blue", high = "white")


#~~~~~~~~~~~~~~~~~~
# 2014
TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend14 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend14 <- sapply(KMA2014, function(geo) {
  Legend14[sapply(TempMix14[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempProportionColors14 <- sapply(KMA2014, function(geo) {
  ProportionColors[sapply(TempMix14[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates14 <- KMA2014Strata_EstimatesStats

RegionalEstimates14 <- KMA2014Strata_Regional_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2015
TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend15 <- setNames(object = c("June 1-July 3", "July 4-August 1", "August 2-29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend15 <- sapply(KMA2015, function(geo) {
  Legend15[sapply(TempMix15[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempProportionColors15 <- sapply(KMA2015, function(geo) {
  ProportionColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates15 <- KMA2015Strata_EstimatesStats

RegionalEstimates15 <- KMA2015Strata_Regional_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2016
TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

Legend16 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"),
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend16 <- sapply(KMA2016, function(geo) {
  Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempProportionColors16 <- sapply(KMA2016, function(geo) {
  ProportionColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates16 <- KMA2016Strata_EstimatesStats

RegionalEstimates16 <- KMA2016Strata_Regional_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
SubRegGroups <- KMA14GroupsPC[2:11]
cex.lab <- 1.5
cex.xaxis <- 0.5  # 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)

sapply(GeoMix, function(geomix) {
  
  emf(file = paste("Figures/All Years/", filenames[geomix], ".emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")
  
  
  layout(mat = layoutmat, widths = c(0.075, 0.375, 0.625), heights = c(0.9, 0.9, 0.9, 0.15))
  par(mar = rep(0, 4))
  par(family = "times")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Y-axis label
  plot.new()
  text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2014 Barplot Regional
  if(length(grep(pattern = geomix, x = names(TempMix14), value = TRUE)) == 0) {
    par(mar = c(1, 1, 1, 0))
    Barplot14 <- barplot2(height = rep(0, 6), beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd, ci.l = rep(0, 6),  ci.u = rep(0, 6), 
                          ylim = c(0, 100), col = 1, yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
    abline(h = 0, xpd = FALSE)
    text(x = mean(Barplot14), y = 50, labels = "No estimates avaialable for 2014", cex = cex.leg)
  } else {
    geomix14 <- grep(pattern = geomix, x = names(TempMix14), value = TRUE)
    par(mar = c(1, 1, 1, 0))
    Barplot14 <- barplot2(height = t(sapply(TempMix14[[geomix14]], function(tempmix) {RegionalEstimates14[[tempmix]][, "median"]})) * 100, 
                          beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                          ci.l = t(sapply(TempMix14[[geomix14]], function(tempmix) {RegionalEstimates14[[tempmix]][, "5%"]})) * 100, 
                          ci.u = t(sapply(TempMix14[[geomix14]], function(tempmix) {RegionalEstimates14[[tempmix]][, "95%"]})) * 100, 
                          ylim = c(0, 100), col = TempProportionColors14[[geomix14]], yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
    # legend(legend = TempLegend14[[geomix14]], x = "topleft", fill = TempProportionColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2014")
    abline(h = 0, xpd = FALSE)
  }
  #~~~~~~~~~~~~~~~~
  ## 2014 Barplot SubRegional
  if(length(grep(pattern = geomix, x = names(TempMix14), value = TRUE)) == 0) {
    par(mar = c(1, 1, 1, 0))
    Barplot14 <- barplot2(height = rep(0, 10), beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd, ci.l = rep(0, 10),  ci.u = rep(0, 10), 
                          ylim = c(0, 100), col = 1, yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
    abline(h = 0, xpd = FALSE)
    text(x = mean(Barplot14), y = 50, labels = "No estimates avaialable for 2014", cex = cex.leg)
  } else {
    geomix14 <- grep(pattern = geomix, x = names(TempMix14), value = TRUE)
    par(mar = c(1, 1, 1, 0))
    Barplot14 <- barplot2(height = t(sapply(TempMix14[[geomix14]], function(tempmix) {Estimates14[[tempmix]][SubRegGroups, "median"]})) * 100, 
                          beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                          ci.l = t(sapply(TempMix14[[geomix14]], function(tempmix) {Estimates14[[tempmix]][SubRegGroups, "5%"]})) * 100, 
                          ci.u = t(sapply(TempMix14[[geomix14]], function(tempmix) {Estimates14[[tempmix]][SubRegGroups, "95%"]})) * 100, 
                          ylim = c(0, 100), col = TempProportionColors14[[geomix14]], yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
    legend(legend = TempLegend14[[geomix14]], x = min(Barplot14[, 6]), y = 100, fill = TempProportionColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2014")
    abline(h = 0, xpd = FALSE)
    abline(v = mean(Barplot14[, 2:3]), lty = 2)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2015 Barplot Regional
  geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
  par(mar = c(1, 1, 1, 0))
  Barplot15 <- barplot2(height = t(sapply(TempMix15[[geomix15]], function(tempmix) {RegionalEstimates15[[tempmix]][, "median"]})) * 100, 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix15[[geomix15]], function(tempmix) {RegionalEstimates15[[tempmix]][, "5%"]})) * 100, 
                        ci.u = t(sapply(TempMix15[[geomix15]], function(tempmix) {RegionalEstimates15[[tempmix]][, "95%"]})) * 100, 
                        ylim = c(0, 100), col = TempProportionColors15[[geomix15]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  # legend(legend = TempLegend15[[geomix15]], x = "topleft", fill = TempProportionColors15[[geomix15]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  #~~~~~~~~~~~~~~~~
  ## 2015 Barplot SubRegional
  geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
  par(mar = c(1, 1, 1, 0))
  Barplot15 <- barplot2(height = t(sapply(TempMix15[[geomix15]], function(tempmix) {Estimates15[[tempmix]][SubRegGroups, "median"]})) * 100, 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix15[[geomix15]], function(tempmix) {Estimates15[[tempmix]][SubRegGroups, "5%"]})) * 100, 
                        ci.u = t(sapply(TempMix15[[geomix15]], function(tempmix) {Estimates15[[tempmix]][SubRegGroups, "95%"]})) * 100, 
                        ylim = c(0, 100), col = TempProportionColors15[[geomix15]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix15]], x = min(Barplot15[, 6]), y = 100, fill = TempProportionColors15[[geomix15]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  abline(v = mean(Barplot15[, 2:3]), lty = 2)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2016 Barplot Regional
  geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
  par(mar = c(1, 1, 1, 0))
  Barplot16 <- barplot2(height = t(sapply(TempMix16[[geomix16]], function(tempmix) {RegionalEstimates16[[tempmix]][, "median"]})) * 100, 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix16[[geomix16]], function(tempmix) {RegionalEstimates16[[tempmix]][, "5%"]})) * 100, 
                        ci.u = t(sapply(TempMix16[[geomix16]], function(tempmix) {RegionalEstimates16[[tempmix]][, "95%"]})) * 100, 
                        ylim = c(0, 100), col = TempProportionColors16[[geomix16]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  # legend(legend = TempLegend16[[geomix16]], x = "topleft", fill = TempProportionColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  mtext(text = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), side = 1, line = 1, at = colMeans(Barplot16), adj = 0.5, cex = cex.xaxis)
  # text(x = colMeans(Barplot16), y = -6, labels = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), adj = 0.5, cex = cex.xaxis, srt = 45, xpd = TRUE)
  #~~~~~~~~~~~~~~~~
  ## 2016 Barplot SubRegional
  geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
  par(mar = c(1, 1, 1, 0))
  Barplot16 <- barplot2(height = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][SubRegGroups, "median"]})) * 100, 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][SubRegGroups, "5%"]})) * 100, 
                        ci.u = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][SubRegGroups, "95%"]})) * 100, 
                        ylim = c(0, 100), col = TempProportionColors16[[geomix16]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
  legend(legend = TempLegend16[[geomix16]], x = min(Barplot16[, 6]), y = 100, fill = TempProportionColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  abline(v = mean(Barplot16[, 2:3]), lty = 2)
  mtext(text = Groups2Rows[2:11], side = 1, line = 1, at = colMeans(Barplot16) + c(rep(0, 2), 0.25, 0.5, 0.25, -0.25, 0, 0.25, 0.25, 0), adj = 0.5, cex = cex.xaxis)
  # text(x = colMeans(Barplot16), y = -6, labels = Groups2Rows[2:11], adj = 0.5, cex = cex.xaxis, srt = 45, xpd = TRUE)
  # - c(rep(0.5, 3), rep(0, 2), rep(0.5, 2), rep(0, 3))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Blank Corner
  par(mar = rep(0, 4))
  plot.new()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## x-axis label
  par(mar = rep(0, 4))
  plot.new()
  text(x = 0.5, y = 0.25, labels = "Regional Reporting Group", cex = cex.lab)
  
  par(mar = rep(0, 4))
  plot.new()
  text(x = 0.5, y = 0.25, labels = "Subregional Reporting Group", cex = cex.lab)
  
  dev.off()
})






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

emf(file = "Figures/All Years/Blank.emf", width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 0.375, 0.625), heights = c(0.9, 0.9, 0.9, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2014 Barplot Regional
par(mar = c(1, 1, 1, 0))
Barplot14 <- barplot2(height = matrix(data = 0, nrow = 3, ncol = 6), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = matrix(data = 0, nrow = 3, ncol = 6), 
                      ci.u = matrix(data = 0, nrow = 3, ncol = 6), 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
abline(h = 0, xpd = FALSE)

#~~~~~~~~~~~~~~~~
## 2014 Barplot SubRegional
par(mar = c(1, 1, 1, 0))
Barplot14 <- barplot2(height = matrix(data = 0, nrow = 3, ncol = 10), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = matrix(data = 0, nrow = 3, ncol = 10), 
                      ci.u = matrix(data = 0, nrow = 3, ncol = 10), 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = c("Early", "Middle", "Late"), x = min(Barplot14[, 6]), y = 100, fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="2014")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot14[, 2:3]), lty = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015 Barplot Regional
par(mar = c(1, 1, 1, 0))
Barplot15 <- barplot2(height = matrix(data = 0, nrow = 3, ncol = 6), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = matrix(data = 0, nrow = 3, ncol = 6), 
                      ci.u = matrix(data = 0, nrow = 3, ncol = 6), 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
abline(h = 0, xpd = FALSE)

#~~~~~~~~~~~~~~~~
## 2015 Barplot SubRegional
par(mar = c(1, 1, 1, 0))
Barplot15 <- barplot2(height = matrix(data = 0, nrow = 3, ncol = 10), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = matrix(data = 0, nrow = 3, ncol = 10), 
                      ci.u = matrix(data = 0, nrow = 3, ncol = 10), 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = c("Early", "Middle", "Late"), x = min(Barplot15[, 6]), y = 100, fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="2015")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot15[, 2:3]), lty = 2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2016 Barplot Regional
par(mar = c(1, 1, 1, 0))
Barplot16 <- barplot2(height = matrix(data = 0, nrow = 3, ncol = 6), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = matrix(data = 0, nrow = 3, ncol = 6), 
                      ci.u = matrix(data = 0, nrow = 3, ncol = 6), 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
abline(h = 0, xpd = FALSE)
mtext(text = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), side = 1, line = 1, at = colMeans(Barplot16), adj = 0.5, cex = cex.xaxis)
#~~~~~~~~~~~~~~~~
## 2016 Barplot SubRegional
par(mar = c(1, 1, 1, 0))
Barplot16 <- barplot2(height = matrix(data = 0, nrow = 3, ncol = 10), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = matrix(data = 0, nrow = 3, ncol = 10), 
                      ci.u = matrix(data = 0, nrow = 3, ncol = 10), 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = c("Early", "Middle", "Late"), x = min(Barplot16[, 6]), y = 100, fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="2016")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot16[, 2:3]), lty = 2)
mtext(text = Groups2Rows[2:11], side = 1, line = 1, at = colMeans(Barplot16) + c(rep(0, 2), 0.25, 0.5, 0.25, -0.25, 0, 0.25, 0.25, 0), adj = 0.5, cex = cex.xaxis)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Regional Reporting Group", cex = cex.lab)

par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Subregional Reporting Group", cex = cex.lab)

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Where were stocks caught? ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")
KMA2015Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")
KMA2016Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_EstimatesStats.txt")

KMA2014Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
KMA2015Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_HarvestEstimatesStats.txt")
KMA2016Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_HarvestEstimatesStats.txt")

KMA2014_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_Stratified_EstimatesStats.txt")
KMA2015_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_Stratified_EstimatesStats.txt")
KMA2016_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_Stratified_EstimatesStats.txt")

KMA2014_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_Stratified_HarvestEstimatesStats.txt")
KMA2015_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_Stratified_HarvestEstimatesStats.txt")
KMA2016_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_Stratified_HarvestEstimatesStats.txt")


EstimateStats <- c(KMA2014Strata_EstimatesStats, KMA2015Strata_EstimatesStats, KMA2016Strata_EstimatesStats)
str(EstimateStats)
length(EstimateStats)


HarvestEstimatesStats <- c(KMA2014Strata_HarvestEstimatesStats, KMA2015Strata_HarvestEstimatesStats, KMA2016Strata_HarvestEstimatesStats)

StrataMedians <- sapply(KMA14GroupsPC, function(RG) {
  sapply(EstimateStats, function(strata) {
    round(strata[RG, "median"], 3)
  })
})

HarvestMedians <- sapply(KMA14GroupsPC, function(RG) {
  sapply(HarvestEstimatesStats, function(strata) {
    round(strata[RG, "median"])
  })
})

# If >= 5% of mixture, where?
sapply(KMA14GroupsPC, function(RG) {which(StrataMedians[, RG] >= 0.05)})

# Plot histogram
invisible(sapply(KMA14GroupsPC, function(RG) {
  hist(StrataMedians[, RG], col = 8, main = RG, ylab = "Proportion", xlim = c(0, 1), breaks = seq(0, 1, 0.01))
  abline(v = 0.05, lwd = 3, col = 2)}))

# If >= 5% of mixtures, which ones and how much?
sapply(KMA14GroupsPC, function(RG) {
  list("n" = sum(StrataMedians[, RG] >= 0.05),
       "Strata" = rbind("Percent" = StrataMedians[StrataMedians[, RG] >= 0.05, RG] * 100,
                        "Harvest" = HarvestMedians[StrataMedians[, RG] >= 0.05, RG]),
       "Annual" = round(rbind("2014" = c(KMA2014_Annual_Stratified_EstimatesStats[RG, "median"] * 100, KMA2014_Annual_Stratified_HarvestEstimatesStats[RG, "median"]),
                              "2015" = c(KMA2015_Annual_Stratified_EstimatesStats[RG, "median"] * 100, KMA2015_Annual_Stratified_HarvestEstimatesStats[RG, "median"]),
                              "2016" = c(KMA2016_Annual_Stratified_EstimatesStats[RG, "median"] * 100, KMA2016_Annual_Stratified_HarvestEstimatesStats[RG, "median"])), 1))}, simplify = FALSE )

apply(StrataMedians, 2, function(x) {x>=0.05})




z <- cbind(Annual2014_Stratified_HarvestEstimates,
           Annual2015_Stratified_HarvestEstimates,
           Annual2016_Stratified_HarvestEstimates)
str(z)
sum(z["Prince William Sound", c(5:6, 11:12, 17:18)]) / sum(z["Prince William Sound", ])




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Harvest for KMA Strata Mixtures 14RG/6RG ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
KMA2015Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_HarvestEstimatesStats.txt")
KMA2016Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_HarvestEstimatesStats.txt")
str(KMA2014Strata_HarvestEstimatesStats)

KMA2014Strata_Regional_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_Regional_HarvestEstimatesStats.txt")
KMA2015Strata_Regional_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_Regional_HarvestEstimatesStats.txt")
KMA2016Strata_Regional_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016Strata_Regional_HarvestEstimatesStats.txt")
str(KMA2014Strata_Regional_HarvestEstimatesStats)


sapply(c(KMA2014Strata_HarvestEstimatesStats, KMA2015Strata_HarvestEstimatesStats, KMA2016Strata_HarvestEstimatesStats), function(mix) {table(mix[, "GR"] > 1.2)})

round(KMA2014Strata_HarvestEstimatesStats$SUYAKC14_2_Middle[, c("mean", "GR")], 3)
round(KMA2015Strata_HarvestEstimatesStats$SUGANC15_3_Late[, c("mean", "GR")], 3)
round(KMA2016Strata_HarvestEstimatesStats$SKARLC16_1_Early[, c("mean", "GR")], 3)

sapply(c(KMA2014Strata_Regional_HarvestEstimatesStats, KMA2015Strata_Regional_HarvestEstimatesStats, KMA2016Strata_Regional_HarvestEstimatesStats), function(mix) {table(mix[, "GR"] > 1.2)})





# Three barplot layout
layoutmat <- matrix(data=c(  1, 2, 3,
                             1, 4, 5,
                             1, 6, 7,
                             8, 9, 10), nrow = 4, ncol = 3, byrow = TRUE)



GeoMix <- c("SUGANC", "SUYAKC", "SKARLC", "SAYAKC", "SALITC", "SIGVAC")

filenames <- setNames(object = c("Uganik Harvest 2014-2016", 
                                 "Uyak Harvest 2014-2016",
                                 "Karluk Harvest 2014-2016",
                                 "Ayakulik Harvest 2014-2016",
                                 "Alitak Harvest 2014-2016",
                                 "Igvak Harvest 2014-2016"), nm = GeoMix)

# If showing proportions (percetages) use blue, otherwise green as "low"
HarvestColors <- colorpanel(n = 3, low = "green", high = "white")

HarvestColors <- c("darkgreen", "green", "white")

#~~~~~~~~~~~~~~~~~~
# 2014
TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_HarvestEstimatesStats), value = TRUE)}, simplify = FALSE)

Legend14 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend14 <- sapply(KMA2014, function(geo) {
  Legend14[sapply(TempMix14[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempHarvestColors14 <- sapply(KMA2014, function(geo) {
  HarvestColors[sapply(TempMix14[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

HarvestEstimates14 <- KMA2014Strata_HarvestEstimatesStats

RegionalHarvestEstimates14 <- KMA2014Strata_Regional_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2015
TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_HarvestEstimatesStats), value = TRUE)}, simplify = FALSE)

Legend15 <- setNames(object = c("June 1-July 3", "July 4-August 1", "August 2-29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend15 <- sapply(KMA2015, function(geo) {
  Legend15[sapply(TempMix15[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempHarvestColors15 <- sapply(KMA2015, function(geo) {
  HarvestColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

HarvestEstimates15 <- KMA2015Strata_HarvestEstimatesStats

RegionalHarvestEstimates15 <- KMA2015Strata_Regional_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2016
TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_HarvestEstimatesStats), value = TRUE)}, simplify = FALSE)

Legend16 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"),
                     nm = c("1_Early", "2_Middle", "3_Late"))
TempLegend16 <- sapply(KMA2016, function(geo) {
  Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
}, simplify = FALSE)


TempHarvestColors16 <- sapply(KMA2016, function(geo) {
  HarvestColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

HarvestEstimates16 <- KMA2016Strata_HarvestEstimatesStats

RegionalHarvestEstimates16 <- KMA2016Strata_Regional_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
SubRegGroups <- KMA14GroupsPC[2:11]
cex.lab <- 1.5
cex.xaxis <- 0.5  # 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5
ymax <- 200000  # max(sapply(c(RegionalHarvestEstimates14, RegionalHarvestEstimates15, RegionalHarvestEstimates16), function(strata) {strata[, "95%"]}))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)

sapply(GeoMix, function(geomix) {
  
  emf(file = paste("Figures/All Years/", filenames[geomix], "DarkGreen.emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")
  
  
  layout(mat = layoutmat, widths = c(0.075, 0.375, 0.625), heights = c(0.9, 0.9, 0.9, 0.15))
  par(mar = rep(0, 4))
  par(family = "times")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Y-axis label
  plot.new()
  text(x = 0.25, y = 0.5, labels = "Number of Fish Harvested (Thousands)", srt = 90, cex = cex.lab)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2014 Barplot Regional
  if(length(grep(pattern = geomix, x = names(TempMix14), value = TRUE)) == 0) {
    par(mar = c(1, 1, 1, 0))
    Barplot14 <- barplot2(height = rep(0, 6), beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd, ci.l = rep(0, 6),  ci.u = rep(0, 6), 
                          ylim = c(0, ymax), col = 1, yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, ymax, 50000), labels = formatC(x = seq(0, ymax, 50000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
    abline(h = 0, xpd = FALSE)
    text(x = mean(Barplot14), y = 100000, labels = "No estimates avaialable for 2014", cex = cex.leg)
  } else {
    geomix14 <- grep(pattern = geomix, x = names(TempMix14), value = TRUE)
    par(mar = c(1, 1, 1, 0))
    Barplot14 <- barplot2(height = t(sapply(TempMix14[[geomix14]], function(tempmix) {RegionalHarvestEstimates14[[tempmix]][, "median"]})), 
                          beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                          ci.l = t(sapply(TempMix14[[geomix14]], function(tempmix) {RegionalHarvestEstimates14[[tempmix]][, "5%"]})), 
                          ci.u = t(sapply(TempMix14[[geomix14]], function(tempmix) {RegionalHarvestEstimates14[[tempmix]][, "95%"]})), 
                          ylim = c(0, ymax), col = TempHarvestColors14[[geomix14]], yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, ymax, 50000), labels = formatC(x = seq(0, ymax, 50000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
    # legend(legend = TempLegend14[[geomix14]], x = "topleft", fill = TempHarvestColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2014")
    abline(h = 0, xpd = FALSE)
  }
  #~~~~~~~~~~~~~~~~
  ## 2014 Barplot SubRegional
  if(length(grep(pattern = geomix, x = names(TempMix14), value = TRUE)) == 0) {
    par(mar = c(1, 1, 1, 0))
    Barplot14 <- barplot2(height = rep(0, 10), beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd, ci.l = rep(0, 10),  ci.u = rep(0, 10), 
                          ylim = c(0, ymax), col = 1, yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, ymax, 50000), labels = FALSE, cex.axis = cex.yaxis)
    abline(h = 0, xpd = FALSE)
    text(x = mean(Barplot14), y = 100000, labels = "No estimates avaialable for 2014", cex = cex.leg)
  } else {
    geomix14 <- grep(pattern = geomix, x = names(TempMix14), value = TRUE)
    par(mar = c(1, 1, 1, 0))
    Barplot14 <- barplot2(height = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][SubRegGroups, "median"]})), 
                          beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                          ci.l = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][SubRegGroups, "5%"]})), 
                          ci.u = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][SubRegGroups, "95%"]})), 
                          ylim = c(0, ymax), col = TempHarvestColors14[[geomix14]], yaxt = "n", xaxt = 'n')
    axis(side = 2, at = seq(0, ymax, 50000), labels = FALSE, cex.axis = cex.yaxis)
    legend(legend = TempLegend14[[geomix14]], x = min(Barplot14[, 6]), y = ymax, fill = TempHarvestColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2014")
    abline(h = 0, xpd = FALSE)
    abline(v = mean(Barplot14[, 2:3]), lty = 2)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2015 Barplot Regional
  geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
  par(mar = c(1, 1, 1, 0))
  Barplot15 <- barplot2(height = t(sapply(TempMix15[[geomix15]], function(tempmix) {RegionalHarvestEstimates15[[tempmix]][, "median"]})), 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix15[[geomix15]], function(tempmix) {RegionalHarvestEstimates15[[tempmix]][, "5%"]})), 
                        ci.u = t(sapply(TempMix15[[geomix15]], function(tempmix) {RegionalHarvestEstimates15[[tempmix]][, "95%"]})), 
                        ylim = c(0, ymax), col = TempHarvestColors15[[geomix15]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 50000), labels = formatC(x = seq(0, ymax, 50000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  # legend(legend = TempLegend15[[geomix15]], x = "topleft", fill = TempHarvestColors15[[geomix15]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  #~~~~~~~~~~~~~~~~
  ## 2015 Barplot SubRegional
  geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
  par(mar = c(1, 1, 1, 0))
  Barplot15 <- barplot2(height = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][SubRegGroups, "median"]})), 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][SubRegGroups, "5%"]})), 
                        ci.u = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][SubRegGroups, "95%"]})), 
                        ylim = c(0, ymax), col = TempHarvestColors15[[geomix15]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 50000), labels = FALSE, cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix15]], x = min(Barplot15[, 6]), y = ymax, fill = TempHarvestColors15[[geomix15]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  abline(v = mean(Barplot15[, 2:3]), lty = 2)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2016 Barplot Regional
  geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
  par(mar = c(1, 1, 1, 0))
  Barplot16 <- barplot2(height = t(sapply(TempMix16[[geomix16]], function(tempmix) {RegionalHarvestEstimates16[[tempmix]][, "median"]})), 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix16[[geomix16]], function(tempmix) {RegionalHarvestEstimates16[[tempmix]][, "5%"]})), 
                        ci.u = t(sapply(TempMix16[[geomix16]], function(tempmix) {RegionalHarvestEstimates16[[tempmix]][, "95%"]})), 
                        ylim = c(0, ymax), col = TempHarvestColors16[[geomix16]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 50000), labels = formatC(x = seq(0, ymax, 50000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  # legend(legend = TempLegend16[[geomix16]], x = "topleft", fill = TempHarvestColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  mtext(text = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), side = 1, line = 1, at = colMeans(Barplot16), adj = 0.5, cex = cex.xaxis)
  # text(x = colMeans(Barplot16), y = -6, labels = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), adj = 0.5, cex = cex.xaxis, srt = 45, xpd = TRUE)
  #~~~~~~~~~~~~~~~~
  ## 2016 Barplot SubRegional
  geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
  par(mar = c(1, 1, 1, 0))
  Barplot16 <- barplot2(height = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][SubRegGroups, "median"]})), 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][SubRegGroups, "5%"]})), 
                        ci.u = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][SubRegGroups, "95%"]})), 
                        ylim = c(0, ymax), col = TempHarvestColors16[[geomix16]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 50000), labels = FALSE, cex.axis = cex.yaxis)
  legend(legend = TempLegend16[[geomix16]], x = min(Barplot16[, 6]), y = ymax, fill = TempHarvestColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  abline(v = mean(Barplot16[, 2:3]), lty = 2)
  mtext(text = Groups2Rows[2:11], side = 1, line = 1, at = colMeans(Barplot16) + c(rep(0, 2), 0.25, 0.5, 0.25, -0.25, 0, 0.25, 0.25, 0), adj = 0.5, cex = cex.xaxis)
  # text(x = colMeans(Barplot16), y = -6, labels = Groups2Rows[2:11], adj = 0.5, cex = cex.xaxis, srt = 45, xpd = TRUE)
  # - c(rep(0.5, 3), rep(0, 2), rep(0.5, 2), rep(0, 3))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Blank Corner
  par(mar = rep(0, 4))
  plot.new()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## x-axis label
  par(mar = rep(0, 4))
  plot.new()
  text(x = 0.5, y = 0.25, labels = "Regional Reporting Group", cex = cex.lab)
  
  par(mar = rep(0, 4))
  plot.new()
  text(x = 0.5, y = 0.25, labels = "Subregional Reporting Group", cex = cex.lab)
  
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Annual KMA Percentages ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_Stratified_EstimatesStats.txt")
KMA2015_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_Stratified_EstimatesStats.txt")
KMA2016_Annual_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_Stratified_EstimatesStats.txt")

KMA2014_Annual_Regional_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_Regional_Stratified_EstimatesStats.txt")
KMA2015_Annual_Regional_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_Regional_Stratified_EstimatesStats.txt")
KMA2016_Annual_Regional_Stratified_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_Regional_Stratified_EstimatesStats.txt")


# Three barplot layout
layoutmat <- matrix(data=c(  1, 2, 3,
                             4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)

ProportionColors <- colorpanel(n = 3, low = "blue", high = "white")

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
SubRegGroups <- KMA14GroupsPC[2:11]
cex.lab <- 1.5
cex.xaxis <- 0.5  # 0.5
cex.yaxis <- 1.3
cex.leg <- 1.8  # 1.1
ci.lwd <- 2.5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)



emf(file = "Figures/All Years/KMA Proportions 2014-2016.emf", width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 0.375, 0.6725), heights = c(2.7, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of KMA Harvest", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplot Regional
par(mar = c(1, 1, 1.5, 0))
Barplot <- barplot2(height = t(cbind(KMA2014_Annual_Regional_Stratified_EstimatesStats[, "median"],
                                     KMA2015_Annual_Regional_Stratified_EstimatesStats[, "median"],
                                     KMA2016_Annual_Regional_Stratified_EstimatesStats[, "median"])) * 100, 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = t(cbind(KMA2014_Annual_Regional_Stratified_EstimatesStats[, "5%"],
                                   KMA2015_Annual_Regional_Stratified_EstimatesStats[, "5%"],
                                   KMA2016_Annual_Regional_Stratified_EstimatesStats[, "5%"])) * 100, 
                    ci.u = t(cbind(KMA2014_Annual_Regional_Stratified_EstimatesStats[, "95%"],
                                   KMA2015_Annual_Regional_Stratified_EstimatesStats[, "95%"],
                                   KMA2016_Annual_Regional_Stratified_EstimatesStats[, "95%"])) * 100, 
                    ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
# legend(legend = 2014:2016, x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="")
abline(h = 0, xpd = FALSE)
mtext(text = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), side = 1, line = 1, at = colMeans(Barplot), adj = 0.5, cex = cex.xaxis)

## Barplot Subregional
par(mar = c(1, 1, 1.5, 0))
Barplot <- barplot2(height = t(cbind(KMA2014_Annual_Stratified_EstimatesStats[SubRegGroups, "median"],
                                     KMA2015_Annual_Stratified_EstimatesStats[SubRegGroups, "median"],
                                     KMA2016_Annual_Stratified_EstimatesStats[SubRegGroups, "median"])) * 100, 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = t(cbind(KMA2014_Annual_Stratified_EstimatesStats[SubRegGroups, "5%"],
                                   KMA2015_Annual_Stratified_EstimatesStats[SubRegGroups, "5%"],
                                   KMA2016_Annual_Stratified_EstimatesStats[SubRegGroups, "5%"])) * 100, 
                    ci.u = t(cbind(KMA2014_Annual_Stratified_EstimatesStats[SubRegGroups, "95%"],
                                   KMA2015_Annual_Stratified_EstimatesStats[SubRegGroups, "95%"],
                                   KMA2016_Annual_Stratified_EstimatesStats[SubRegGroups, "95%"])) * 100, 
                    ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = 2014:2016, x = min(Barplot[, 5]), y = 100, fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="June 1-August 29")
#legend(legend = rep("", 3), x = min(Barplot[, 6]), y = 100, border = "black", bty = "n", cex = cex.title, title="June 1-August 29")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot[, 2:3]), lty = 2)
mtext(text = Groups2Rows[2:11], side = 1, line = 1, at = colMeans(Barplot) + c(rep(0, 2), 0.25, 0.5, 0.25, -0.25, 0, 0.25, 0.25, 0), adj = 0.5, cex = cex.xaxis)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Regional Reporting Group", cex = cex.lab)

par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Subregional Reporting Group", cex = cex.lab)


dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


emf(file = "Figures/All Years/KMA Proportions 2014-2016 Blank.emf", width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 0.375, 0.6725), heights = c(2.7, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of KMA Harvest", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplot Regional
par(mar = c(1, 1, 1.5, 0))
Barplot <- barplot2(height = t(cbind(rep(0, 6),
                                     rep(0, 6),
                                     rep(0, 6))) * 100, 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = t(cbind(rep(0, 6),
                                   rep(0, 6),
                                   rep(0, 6))) * 100, 
                    ci.u = t(cbind(rep(0, 6),
                                   rep(0, 6),
                                   rep(0, 6))) * 100, 
                    ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
# legend(legend = 2014:2016, x = "topleft", fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="")
abline(h = 0, xpd = FALSE)
mtext(text = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), side = 1, line = 1, at = colMeans(Barplot), adj = 0.5, cex = cex.xaxis)

## Barplot Subregional
par(mar = c(1, 1, 1.5, 0))
Barplot <- barplot2(height = t(cbind(rep(0, 10),
                                     rep(0, 10),
                                     rep(0, 10))) * 100, 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = t(cbind(rep(0, 10),
                                   rep(0, 10),
                                   rep(0, 10))) * 100, 
                    ci.u = t(cbind(rep(0, 10),
                                   rep(0, 10),
                                   rep(0, 10))) * 100, 
                    ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = 2014:2016, x = min(Barplot[, 5]), y = 100, fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="June 1-August 29", adj = 0)
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot[, 2:3]), lty = 2)
mtext(text = Groups2Rows[2:11], side = 1, line = 1, at = colMeans(Barplot) + c(rep(0, 2), 0.25, 0.5, 0.25, -0.25, 0, 0.25, 0.25, 0), adj = 0.5, cex = cex.xaxis)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Regional Reporting Group", cex = cex.lab)

par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Subregional Reporting Group", cex = cex.lab)


dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Annual KMA Harvest ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_Stratified_HarvestEstimatesStats.txt")
KMA2015_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_Stratified_HarvestEstimatesStats.txt")
KMA2016_Annual_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_Stratified_HarvestEstimatesStats.txt")

KMA2014_Annual_Regional_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Annual_Regional_Stratified_HarvestEstimatesStats.txt")
KMA2015_Annual_Regional_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Annual_Regional_Stratified_HarvestEstimatesStats.txt")
KMA2016_Annual_Regional_Stratified_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Annual_Regional_Stratified_HarvestEstimatesStats.txt")


# Three barplot layout
layoutmat <- matrix(data=c(  1, 2, 3,
                             4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE)

HarvestColors <- colorpanel(n = 3, low = "green", high = "white")

HarvestColors <- c("darkgreen", "green", "white")

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
SubRegGroups <- KMA14GroupsPC[2:11]
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.yaxis <- 1.3
cex.leg <- 1.8  # 1.1
ci.lwd <- 2.5
ymax <- 1400000  # max(sapply(list(KMA2014_Annual_Regional_Stratified_HarvestEstimatesStats, KMA2015_Annual_Regional_Stratified_HarvestEstimatesStats, KMA2016_Annual_Regional_Stratified_HarvestEstimatesStats), function(strata) {strata[, "95%"]}))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)



emf(file = "Figures/All Years/KMA Harvest 2014-2016DarkGreen.emf", width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 0.375, 0.6725), heights = c(2.7, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Number of Fish Harvested (Thousands)", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Barplot Regional
par(mar = c(1, 1, 1.5, 0))
Barplot <- barplot2(height = t(cbind(KMA2014_Annual_Regional_Stratified_HarvestEstimatesStats[, "median"],
                                     KMA2015_Annual_Regional_Stratified_HarvestEstimatesStats[, "median"],
                                     KMA2016_Annual_Regional_Stratified_HarvestEstimatesStats[, "median"])), 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = t(cbind(KMA2014_Annual_Regional_Stratified_HarvestEstimatesStats[, "5%"],
                                   KMA2015_Annual_Regional_Stratified_HarvestEstimatesStats[, "5%"],
                                   KMA2016_Annual_Regional_Stratified_HarvestEstimatesStats[, "5%"])), 
                    ci.u = t(cbind(KMA2014_Annual_Regional_Stratified_HarvestEstimatesStats[, "95%"],
                                   KMA2015_Annual_Regional_Stratified_HarvestEstimatesStats[, "95%"],
                                   KMA2016_Annual_Regional_Stratified_HarvestEstimatesStats[, "95%"])), 
                    ylim = c(0, ymax), col = HarvestColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 200000), labels = formatC(x = seq(0, ymax, 200000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
# legend(legend = 2014:2016, x = "topleft", fill = HarvestColors, border = "black", bty = "n", cex = cex.leg, title="")
abline(h = 0, xpd = FALSE)

mtext(text = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), side = 1, line = 1, at = colMeans(Barplot), adj = 0.5, cex = cex.xaxis)

## Barplot Subregional
par(mar = c(1, 1, 1.5, 0))
Barplot <- barplot2(height = t(cbind(KMA2014_Annual_Stratified_HarvestEstimatesStats[SubRegGroups, "median"],
                                     KMA2015_Annual_Stratified_HarvestEstimatesStats[SubRegGroups, "median"],
                                     KMA2016_Annual_Stratified_HarvestEstimatesStats[SubRegGroups, "median"])), 
                    beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                    ci.l = t(cbind(KMA2014_Annual_Stratified_HarvestEstimatesStats[SubRegGroups, "5%"],
                                   KMA2015_Annual_Stratified_HarvestEstimatesStats[SubRegGroups, "5%"],
                                   KMA2016_Annual_Stratified_HarvestEstimatesStats[SubRegGroups, "5%"])), 
                    ci.u = t(cbind(KMA2014_Annual_Stratified_HarvestEstimatesStats[SubRegGroups, "95%"],
                                   KMA2015_Annual_Stratified_HarvestEstimatesStats[SubRegGroups, "95%"],
                                   KMA2016_Annual_Stratified_HarvestEstimatesStats[SubRegGroups, "95%"])), 
                    ylim = c(0, ymax), col = HarvestColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 200000), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = 2014:2016, x = min(Barplot[, 5]), y = ymax, fill = HarvestColors, border = "black", bty = "n", cex = cex.leg, title="June 1-August 29")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot[, 2:3]), lty = 2)
mtext(text = Groups2Rows[2:11], side = 1, line = 1, at = colMeans(Barplot) + c(rep(0, 2), 0.25, 0.5, 0.25, -0.25, 0, 0.25, 0.25, 0), adj = 0.5, cex = cex.xaxis)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Regional Reporting Group", cex = cex.lab)

par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Subregional Reporting Group", cex = cex.lab)

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Annual Temporal KMA Percentages ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014_Temporal_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Temporal_Annual_EstimatesStats.txt")
KMA2015_Temporal_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Temporal_Annual_EstimatesStats.txt")
KMA2016_Temporal_Annual_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Temporal_Annual_EstimatesStats.txt")
str(KMA2014_Temporal_Annual_EstimatesStats)

KMA2014_Temporal_Annual_Regional_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Temporal_Annual_Regional_EstimatesStats.txt")
KMA2015_Temporal_Annual_Regional_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Temporal_Annual_Regional_EstimatesStats.txt")
KMA2016_Temporal_Annual_Regional_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Temporal_Annual_Regional_EstimatesStats.txt")
str(KMA2016_Temporal_Annual_Regional_EstimatesStats)


# Three barplot layout
layoutmat <- matrix(data=c(  1, 2, 3,
                             1, 4, 5,
                             1, 6, 7,
                             8, 9, 10), nrow = 4, ncol = 3, byrow = TRUE)


# If showing proportions (percetages) use blue, otherwise green as "low"
ProportionColors <- colorpanel(n = 3, low = "blue", high = "white")


#~~~~~~~~~~~~~~~~~~
# 2014
Legend14 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))

Estimates14 <- KMA2014_Temporal_Annual_EstimatesStats

RegionalEstimates14 <- KMA2014_Temporal_Annual_Regional_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2015
Legend15 <- setNames(object = c("June 1-July 3", "July 4-August 1", "August 2-29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))

Estimates15 <- KMA2015_Temporal_Annual_EstimatesStats

RegionalEstimates15 <- KMA2015_Temporal_Annual_Regional_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2016
Legend16 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"),
                     nm = c("1_Early", "2_Middle", "3_Late"))

Estimates16 <- KMA2016_Temporal_Annual_EstimatesStats

RegionalEstimates16 <- KMA2016_Temporal_Annual_Regional_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
SubRegGroups <- KMA14GroupsPC[2:11]
cex.lab <- 1.5
cex.xaxis <- 0.5  # 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)


emf(file = paste("Figures/All Years/KMA Temporal Proportions 2014-2016.emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 0.375, 0.625), heights = c(0.9, 0.9, 0.9, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2014 Barplot Regional
par(mar = c(1, 1, 1, 0))
Barplot14 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalEstimates14[[tempmix]][, "median"]})) * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalEstimates14[[tempmix]][, "5%"]})) * 100, 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalEstimates14[[tempmix]][, "95%"]})) * 100, 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
abline(h = 0, xpd = FALSE)

#~~~~~~~~~~~~~~~~
## 2014 Barplot SubRegional

par(mar = c(1, 1, 1, 0))
Barplot14 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {Estimates14[[tempmix]][SubRegGroups, "median"]})) * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {Estimates14[[tempmix]][SubRegGroups, "5%"]})) * 100, 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {Estimates14[[tempmix]][SubRegGroups, "95%"]})) * 100, 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = Legend14, x = min(Barplot14[, 6]), y = 100, fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="2014")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot14[, 2:3]), lty = 2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015 Barplot Regional
par(mar = c(1, 1, 1, 0))
Barplot15 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalEstimates15[[tempmix]][, "median"]})) * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalEstimates15[[tempmix]][, "5%"]})) * 100, 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalEstimates15[[tempmix]][, "95%"]})) * 100, 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
abline(h = 0, xpd = FALSE)
#~~~~~~~~~~~~~~~~
## 2015 Barplot SubRegional
geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
par(mar = c(1, 1, 1, 0))
Barplot15 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {Estimates15[[tempmix]][SubRegGroups, "median"]})) * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {Estimates15[[tempmix]][SubRegGroups, "5%"]})) * 100, 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {Estimates15[[tempmix]][SubRegGroups, "95%"]})) * 100, 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = Legend15, x = min(Barplot15[, 6]), y = 100, fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="2015")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot15[, 2:3]), lty = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2016 Barplot Regional
geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
par(mar = c(1, 1, 1, 0))
Barplot16 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalEstimates16[[tempmix]][, "median"]})) * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalEstimates16[[tempmix]][, "5%"]})) * 100, 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalEstimates16[[tempmix]][, "95%"]})) * 100, 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
abline(h = 0, xpd = FALSE)
mtext(text = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), side = 1, line = 1, at = colMeans(Barplot16), adj = 0.5, cex = cex.xaxis)
# text(x = colMeans(Barplot16), y = -6, labels = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), adj = 0.5, cex = cex.xaxis, srt = 45, xpd = TRUE)
#~~~~~~~~~~~~~~~~
## 2016 Barplot SubRegional
geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
par(mar = c(1, 1, 1, 0))
Barplot16 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {Estimates16[[tempmix]][SubRegGroups, "median"]})) * 100, 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {Estimates16[[tempmix]][SubRegGroups, "5%"]})) * 100, 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {Estimates16[[tempmix]][SubRegGroups, "95%"]})) * 100, 
                      ylim = c(0, 100), col = ProportionColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 100, 25), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = Legend16, x = min(Barplot16[, 6]), y = 100, fill = ProportionColors, border = "black", bty = "n", cex = cex.leg, title="2016")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot16[, 2:3]), lty = 2)
mtext(text = Groups2Rows[2:11], side = 1, line = 1, at = colMeans(Barplot16) + c(rep(0, 2), 0.25, 0.5, 0.25, -0.25, 0, 0.25, 0.25, 0), adj = 0.5, cex = cex.xaxis)
# text(x = colMeans(Barplot16), y = -6, labels = Groups2Rows[2:11], adj = 0.5, cex = cex.xaxis, srt = 45, xpd = TRUE)
# - c(rep(0.5, 3), rep(0, 2), rep(0.5, 2), rep(0, 3))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Regional Reporting Group", cex = cex.lab)

par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Subregional Reporting Group", cex = cex.lab)

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Harvest for KMA Strata Mixtures 14RG/6RG ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014_Temporal_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Temporal_Annual_HarvestEstimatesStats.txt")
KMA2015_Temporal_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Temporal_Annual_HarvestEstimatesStats.txt")
KMA2016_Temporal_Annual_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Temporal_Annual_HarvestEstimatesStats.txt")
str(KMA2014_Temporal_Annual_HarvestEstimatesStats)

KMA2014_Temporal_Annual_Regional_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014_Temporal_Annual_Regional_HarvestEstimatesStats.txt")
KMA2015_Temporal_Annual_Regional_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015_Temporal_Annual_Regional_HarvestEstimatesStats.txt")
KMA2016_Temporal_Annual_Regional_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2016_Temporal_Annual_Regional_HarvestEstimatesStats.txt")
str(KMA2014_Temporal_Annual_Regional_HarvestEstimatesStats)


# Three barplot layout
layoutmat <- matrix(data=c(  1, 2, 3,
                             1, 4, 5,
                             1, 6, 7,
                             8, 9, 10), nrow = 4, ncol = 3, byrow = TRUE)

# If showing proportions (percetages) use blue, otherwise green as "low"
HarvestColors <- colorpanel(n = 3, low = "green", high = "white")


#~~~~~~~~~~~~~~~~~~
# 2014
Legend14 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))

HarvestEstimates14 <- KMA2014_Temporal_Annual_HarvestEstimatesStats

RegionalHarvestEstimates14 <- KMA2014_Temporal_Annual_Regional_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2015
Legend15 <- setNames(object = c("June 1-July 3", "July 4-August 1", "August 2-29"), 
                     nm = c("1_Early", "2_Middle", "3_Late"))

HarvestEstimates15 <- KMA2015_Temporal_Annual_HarvestEstimatesStats

RegionalHarvestEstimates15 <- KMA2015_Temporal_Annual_Regional_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# 2016
Legend16 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"),
                     nm = c("1_Early", "2_Middle", "3_Late"))

HarvestEstimates16 <- KMA2016_Temporal_Annual_HarvestEstimatesStats

RegionalHarvestEstimates16 <- KMA2016_Temporal_Annual_Regional_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA14GroupsPC
Groups2Rows <- KMA14GroupsPC2Rows
SubRegGroups <- KMA14GroupsPC[2:11]
cex.lab <- 1.5
cex.xaxis <- 0.5  # 0.5
cex.yaxis <- 1.3
cex.leg <- 1.1
ci.lwd <- 2.5
ymax <- 500000  # max(sapply(c(RegionalHarvestEstimates14, RegionalHarvestEstimates15, RegionalHarvestEstimates16), function(strata) {strata[, "95%"]}))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)


emf(file = paste("Figures/All Years/KMA Temporal Harvest 2014-2016.emf", sep = ''), width = 6, height = 5.75, family = "serif", bg = "white")


layout(mat = layoutmat, widths = c(0.075, 0.375, 0.625), heights = c(0.9, 0.9, 0.9, 0.15))
par(mar = rep(0, 4))
par(family = "times")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Y-axis label
plot.new()
text(x = 0.25, y = 0.5, labels = "Number of Fish Harvested (Thousands)", srt = 90, cex = cex.lab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2014 Barplot Regional
par(mar = c(1, 1, 1, 0))
Barplot14 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalHarvestEstimates14[[tempmix]][, "median"]})), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalHarvestEstimates14[[tempmix]][, "5%"]})), 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalHarvestEstimates14[[tempmix]][, "95%"]})), 
                      ylim = c(0, ymax), col = HarvestColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 125000), labels = formatC(x = seq(0, ymax, 125000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
abline(h = 0, xpd = FALSE)

#~~~~~~~~~~~~~~~~
## 2014 Barplot SubRegional
par(mar = c(1, 1, 1, 0))
Barplot14 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {HarvestEstimates14[[tempmix]][SubRegGroups, "median"]})), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {HarvestEstimates14[[tempmix]][SubRegGroups, "5%"]})), 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {HarvestEstimates14[[tempmix]][SubRegGroups, "95%"]})), 
                      ylim = c(0, ymax), col = HarvestColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 125000), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = Legend14, x = min(Barplot14[, 6]), y = ymax, fill = HarvestColors, border = "black", bty = "n", cex = cex.leg, title="2014")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot14[, 2:3]), lty = 2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015 Barplot Regional
geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
par(mar = c(1, 1, 1, 0))
Barplot15 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalHarvestEstimates15[[tempmix]][, "median"]})), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalHarvestEstimates15[[tempmix]][, "5%"]})), 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalHarvestEstimates15[[tempmix]][, "95%"]})), 
                      ylim = c(0, ymax), col = HarvestColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 125000), labels = formatC(x = seq(0, ymax, 125000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
abline(h = 0, xpd = FALSE)
#~~~~~~~~~~~~~~~~
## 2015 Barplot SubRegional
geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
par(mar = c(1, 1, 1, 0))
Barplot15 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {HarvestEstimates15[[tempmix]][SubRegGroups, "median"]})), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {HarvestEstimates15[[tempmix]][SubRegGroups, "5%"]})), 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {HarvestEstimates15[[tempmix]][SubRegGroups, "95%"]})), 
                      ylim = c(0, ymax), col = HarvestColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 125000), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = Legend15, x = min(Barplot15[, 6]), y = ymax, fill = HarvestColors, border = "black", bty = "n", cex = cex.leg, title="2015")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot15[, 2:3]), lty = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2016 Barplot Regional
geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
par(mar = c(1, 1, 1, 0))
Barplot16 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalHarvestEstimates16[[tempmix]][, "median"]})), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalHarvestEstimates16[[tempmix]][, "5%"]})), 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {RegionalHarvestEstimates16[[tempmix]][, "95%"]})), 
                      ylim = c(0, ymax), col = HarvestColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 125000), labels = formatC(x = seq(0, ymax, 125000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
abline(h = 0, xpd = FALSE)
mtext(text = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), side = 1, line = 1, at = colMeans(Barplot16), adj = 0.5, cex = cex.xaxis)
# text(x = colMeans(Barplot16), y = -6, labels = c(Groups2Rows[1], "Chignik\n", "Kodiak\n", Groups2Rows[12:14]), adj = 0.5, cex = cex.xaxis, srt = 45, xpd = TRUE)
#~~~~~~~~~~~~~~~~
## 2016 Barplot SubRegional
geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
par(mar = c(1, 1, 1, 0))
Barplot16 <- barplot2(height = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {HarvestEstimates16[[tempmix]][SubRegGroups, "median"]})), 
                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                      ci.l = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {HarvestEstimates16[[tempmix]][SubRegGroups, "5%"]})), 
                      ci.u = t(sapply(c("Early", "Middle", "Late"), function(tempmix) {HarvestEstimates16[[tempmix]][SubRegGroups, "95%"]})), 
                      ylim = c(0, ymax), col = HarvestColors, yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, ymax, 125000), labels = FALSE, cex.axis = cex.yaxis)
legend(legend = Legend16, x = min(Barplot16[, 6]), y = ymax, fill = HarvestColors, border = "black", bty = "n", cex = cex.leg, title="2016")
abline(h = 0, xpd = FALSE)
abline(v = mean(Barplot16[, 2:3]), lty = 2)
mtext(text = Groups2Rows[2:11], side = 1, line = 1, at = colMeans(Barplot16) + c(rep(0, 2), 0.25, 0.5, 0.25, -0.25, 0, 0.25, 0.25, 0), adj = 0.5, cex = cex.xaxis)
# text(x = colMeans(Barplot16), y = -6, labels = Groups2Rows[2:11], adj = 0.5, cex = cex.xaxis, srt = 45, xpd = TRUE)
# - c(rep(0.5, 3), rep(0, 2), rep(0.5, 2), rep(0, 3))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank Corner
par(mar = rep(0, 4))
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Regional Reporting Group", cex = cex.lab)

par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.25, labels = "Subregional Reporting Group", cex = cex.lab)

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table Regional Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Dates
# DatesStrata2014 <- read.table(file = "Harvest/2014DatesByStrata.txt", header = TRUE, sep = "\t", as.is = TRUE)
# DatesStrata2014.mat <- as.matrix(DatesStrata2014[-1])
# dimnames(DatesStrata2014.mat) <- list(DatesStrata2014$location, c("1_Early", "2_Middle", "3_Late"))
# dput(x = DatesStrata2014.mat, file = "Objects/DatesStrata2014_Final.txt"); rm(DatesStrata2014.mat)
DatesStrata2014_Final <- dget(file = "Objects/DatesStrata2014_Final.txt")


# DatesStrata2015 <- read.table(file = "Harvest/2015DatesByStrata.txt", header = TRUE, sep = "\t", as.is = TRUE)
# DatesStrata2015.mat <- as.matrix(DatesStrata2015[-1])
# dimnames(DatesStrata2015.mat) <- list(DatesStrata2015$location, c("1_Early", "2_Middle", "3_Late"))
# dput(x = DatesStrata2015.mat, file = "Objects/DatesStrata2015_Final.txt"); rm(DatesStrata2015.mat)
DatesStrata2015_Final <- dget(file = "Objects/DatesStrata2015_Final.txt")


# DatesStrata2016 <- read.table(file = "Harvest/2016DatesByStrata.txt", header = TRUE, sep = "\t", as.is = TRUE)
# DatesStrata2016.mat <- as.matrix(DatesStrata2016[-1])
# dimnames(DatesStrata2016.mat) <- list(DatesStrata2016$location, c("1_Early", "2_Middle", "3_Late"))
# dput(x = DatesStrata2016.mat, file = "Objects/DatesStrata2016_Final.txt"); rm(DatesStrata2016.mat)
DatesStrata2016_Final <- dget(file = "Objects/DatesStrata2016_Final.txt")


## Sample sizes
# KMA2014_2016Strata_SampleSizes_Final <- KMA2014_2016Strata_SampleSizes[, "Final"]
# LateLateStrata <- grep(pattern = "LateLate", x = names(KMA2014_2016Strata_SampleSizes_Final))
# KMA2014_2016Strata_SampleSizes_Final[LateLateStrata-1] <- KMA2014_2016Strata_SampleSizes_Final[LateLateStrata-1] + KMA2014_2016Strata_SampleSizes_Final[LateLateStrata]
# KMA2014_2016Strata_SampleSizes_Final_Condense <- KMA2014_2016Strata_SampleSizes_Final[-LateLateStrata]
# dput(x = KMA2014_2016Strata_SampleSizes_Final_Condense, file = "Objects/KMA2014_2016Strata_SampleSizes_Final_Condense.txt")
KMA2014_2016Strata_SampleSizes_Final_Condense <- dget(file = "Objects/KMA2014_2016Strata_SampleSizes_Final_Condense.txt")


## Geographic headers
# GeoHeader <- setNames(object = c(paste0("Alitak (statistical areas 257-10, 20, 50, 60, 70)"),
#                                  paste0("Ayakulik-Halibut Bay (statistical areas 256-10", "\u2013", "256-30)"),
#                                  paste0("Igvak (statistical areas 262-75, 80, 90, 95)"),
#                                  paste0("Karluk-Sturgeon (statistical areas 255-10, 20; 256-40)"),
#                                  paste0("Uganik-Kupreanof (statistical areas 253-10", "\u2013", "253-35)"),
#                                  paste0("Uyak (statistical areas 254-10)", "\u2013", "254-41)")),
#                       nm = unlist(strsplit(x = KMA2016, split = "16")))
GeoHeader <- setNames(object = c(paste0("Alitak (statistical areas 257-10, 20, 50, 60, 70)"),
                                 paste0("Ayakulik-Halibut Bay (statistical areas 256-10, 15, 20, 25, 30)"),
                                 paste0("Igvak (statistical areas 262-75, 80, 90, 95)"),
                                 paste0("Karluk-Sturgeon (statistical areas 255-10, 20; 256-40)"),
                                 paste0("Uganik-Kupreanof (statistical areas 253-11, 12, 13, 14, 31, 32, 33, 34, 35)"),
                                 paste0("Uyak (statistical areas 254-10, 20, 21, 30, 31, 40, 41)")),
                      nm = unlist(strsplit(x = KMA2016, split = "16")))
dput(x = GeoHeader, file = "Objects/GeoHeader.txt")
GeoHeader <- dget(file = "Objects/GeoHeader.txt")



# Get Final Estimates Objects
KMAfinalestimatesobjects <- list.files(path = "Estimates objects/Final", recursive = FALSE)
invisible(sapply(KMAfinalestimatesobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste("Estimates objects/Final", objct, sep = "/")), pos = 1) })); beep(2)
KMAfinalestimatesobjects; rm(KMAfinalestimatesobjects)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Defining caption variables

EstimatesStats <- c(KMA2014Strata_EstimatesStats, KMA2014_Annual_EstimatesStats,
                    KMA2015Strata_EstimatesStats, KMA2015_Annual_EstimatesStats,
                    KMA2016Strata_EstimatesStats, KMA2016_Annual_EstimatesStats)

HarvestEstimatesStats <- c(KMA2014Strata_HarvestEstimatesStats, KMA2014_Annual_HarvestEstimatesStats,
                           KMA2015Strata_HarvestEstimatesStats, KMA2015_Annual_HarvestEstimatesStats,
                           KMA2016Strata_HarvestEstimatesStats, KMA2016_Annual_HarvestEstimatesStats)

Regional_EstimatesStats <- c(KMA2014Strata_Regional_EstimatesStats, KMA2014_Annual_Regional_EstimatesStats,
                             KMA2015Strata_Regional_EstimatesStats, KMA2015_Annual_Regional_EstimatesStats,
                             KMA2016Strata_Regional_EstimatesStats, KMA2016_Annual_Regional_EstimatesStats)

Regional_HarvestEstimatesStats <- c(KMA2014Strata_Regional_HarvestEstimatesStats, KMA2014_Annual_Regional_HarvestEstimatesStats,
                                    KMA2015Strata_Regional_HarvestEstimatesStats, KMA2015_Annual_Regional_HarvestEstimatesStats,
                                    KMA2016Strata_Regional_HarvestEstimatesStats, KMA2016_Annual_Regional_HarvestEstimatesStats)


SheetNames <- sort(names(EstimatesStats))
names(SheetNames) <- SheetNames
mixvec <- SheetNames


# mix <- SheetNames[2]
harvest <- rbind(HarvestByStrata2014_Final, HarvestByStrata2015_Final, HarvestByStrata2016_Final)
dates <- rbind(DatesStrata2014_Final, DatesStrata2015_Final, DatesStrata2016_Final)
sampsize <- KMA2014_2016Strata_SampleSizes_Final_Condense

SubRegGroups <- KMA14GroupsPC[2:11]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

for(mix in SheetNames) {
  
  # String split by "_
  SheetNames.split <- unlist(strsplit(x = mix, split = "_"))
  
  
  # The first element is the geographic area + year
  geomix <- SheetNames.split[1]
  yr <- unlist(strsplit(x = geomix, split = "C"))[2]
  geo <- unlist(strsplit(x = geomix, split = "1"))[1]
  
  
  # If it is not an annual roll-up, then get the strata number + strata name
  if(length(SheetNames.split) > 1) {
    tempmix <- paste(c(SheetNames.split[2], SheetNames.split[3]), collapse = "_")
    
    Caption <- paste("Table X.-Regional and subregional (within Chignik and Kodiak) estimates of stock composition (%) and stock-specific harvest for temporal stratum ",
                     SheetNames.split[2], " (", dates[geomix, tempmix],
                     "; Harvest=", formatC(x = harvest[geomix, tempmix], format = "f", digits = 0, big.mark = ","),
                     "; n=", sampsize[mix], ")", " of ", GeoHeader[geo], ", 20", yr,
                     ". Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and standard deviation (SD).",
                     sep = '')
  } else {
    Caption <- paste("Table X.-Annual regional and subregional (within Chignik and Kodiak) estimates of stock composition (%) and stock-specific harvest for ", GeoHeader[geo], ", 20", yr,
                     ". Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and standard deviation (SD).",
                     sep = '')
  }
  
  Disclaimer <- "Note: Stock composition estimates may not sum to 100% and stock-specific harvest estimates may not sum to the total harvest due to rounding error."
  
  
  TableX <- matrix(data = "", nrow = 24, ncol = 14)
  
  TableX[1, 1] <- Caption
  TableX[2, c(3, 10)] <- c("Stock Composition", "Stock-specific Harvest")
  TableX[3, c(1, 4, 11)] <- c("Reporting Group", rep("90% CI", 2))
  TableX[4, c(1, 2, 3:5, 7:8, 10:14, 6)] <- c("Regional", "Subregional", rep(c("Median", "5%", "95%", "Mean", "SD"), 2), "P=0")
  TableX[5:10, 1] <- c(KMA14GroupsPC[1], "Chignik", "Kodiak", KMA14GroupsPC[12:14])
  TableX[5:10, c(3:5, 7:8)] <- formatC(x = Regional_EstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
  TableX[5:10, 6] <- formatC(x = Regional_EstimatesStats[[mix]][, "P=0"], digits = 2, format = "f")
  TableX[5:10, 10:14] <- formatC(x = Regional_HarvestEstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")
  TableX[11, 12:13] <- c("Total", formatC(x = sum(Regional_HarvestEstimatesStats[[mix]][, "mean"]), digits = 0, format = "f", big.mark = ","))
  TableX[c(13, 16), 1] <- c("Chignik", "Kodiak")
  TableX[c(13:14, 16:23), 2] <- c(SubRegGroups[1:3], "Ayakulik / Frazer", SubRegGroups[5:10])
  TableX[c(13:14, 16:23), c(3:5, 7:8)] <- formatC(x = EstimatesStats[[mix]][SubRegGroups, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
  TableX[c(13:14, 16:23), 6] <- formatC(x = EstimatesStats[[mix]][SubRegGroups, "P=0"], digits = 2, format = "f")
  TableX[c(13:14, 16:23), 10:14] <- formatC(x = HarvestEstimatesStats[[mix]][SubRegGroups, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")
  TableX[24, 1] <- Disclaimer
  
  write.xlsx(x = as.data.frame(TableX), 
             file = "Estimates tables/KMA Sockeye Estimates Tables Regional.xlsx",
             col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = mix)
}; beep(5)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Table Annual KMA Results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMApercent <- setNames(object = c(46.7, 55.2, 62.4), nm = 14:16)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("Estimates tables")
require(xlsx)

for(yr in 14:16){
  
  EstimatesStats <- dget(file = paste0("Estimates objects/Final/KMA20", yr, "_Annual_Stratified_EstimatesStats.txt"))
  HarvestEstimatesStats <- dget(file = paste0("Estimates objects/Final/KMA20", yr, "_Annual_Stratified_HarvestEstimatesStats.txt"))
  
  Regional_EstimatesStats <- dget(file = paste0("Estimates objects/Final/KMA20", yr, "_Annual_Regional_Stratified_EstimatesStats.txt"))
  Regional_HarvestEstimatesStats <- dget(file = paste0("Estimates objects/Final/KMA20", yr, "_Annual_Regional_Stratified_HarvestEstimatesStats.txt"))
  
  
  Caption <- paste("Table X.-Annual regional and subregional (within Chignik and Kodiak) estimates of stock composition (%) and stock-specific harvest for KMA, 20", yr,
                   ". Note that these annual summaries only include strata sampled for this project, which account for ", KMApercent[as.character(yr)],"% of the KMA commercial sockeye salmon harvest. Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and standard deviation (SD).",
                   sep = '')
  
  
  Disclaimer <- "Note: Stock composition estimates may not sum to 100% and stock-specific harvest estimates may not sum to the total harvest due to rounding error."
  
  
  TableX <- matrix(data = "", nrow = 24, ncol = 14)
  
  TableX[1, 1] <- Caption
  TableX[2, c(3, 10)] <- c("Stock Composition", "Stock-specific Harvest")
  TableX[3, c(1, 4, 11)] <- c("Reporting Group", rep("90% CI", 2))
  TableX[4, c(1, 2, 3:5, 7:8, 10:14, 6)] <- c("Regional", "Subregional", rep(c("Median", "5%", "95%", "Mean", "SD"), 2), "P=0")
  TableX[5:10, 1] <- c(KMA14GroupsPC[1], "Chignik", "Kodiak", KMA14GroupsPC[12:14])
  TableX[5:10, c(3:5, 7:8)] <- formatC(x = Regional_EstimatesStats[, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
  TableX[5:10, 6] <- formatC(x = Regional_EstimatesStats[, "P=0"], digits = 2, format = "f")
  TableX[5:10, 10:14] <- formatC(x = Regional_HarvestEstimatesStats[, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")
  TableX[11, 12:13] <- c("Total", formatC(x = sum(Regional_HarvestEstimatesStats[, "mean"]), digits = 0, format = "f", big.mark = ","))
  TableX[c(13, 16), 1] <- c("Chignik", "Kodiak")
  TableX[c(13:14, 16:23), 2] <- c(SubRegGroups[1:3], "Ayakulik / Frazer", SubRegGroups[5:10])
  TableX[c(13:14, 16:23), c(3:5, 7:8)] <- formatC(x = EstimatesStats[SubRegGroups, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
  TableX[c(13:14, 16:23), 6] <- formatC(x = EstimatesStats[SubRegGroups, "P=0"], digits = 2, format = "f")
  TableX[c(13:14, 16:23), 10:14] <- formatC(x = HarvestEstimatesStats[SubRegGroups, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")
  TableX[24, 1] <- Disclaimer
  
  write.xlsx(x = as.data.frame(TableX), 
             file = "Estimates tables/KMA Sockeye Estimates Tables Regional.xlsx",
             col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = paste0("KMA20", yr))
}; beep(5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2014 Late Late Harvest for Appendix
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Defining caption variables

EstimatesStats <- dget(file = "Estimates objects/Round4Mixtures_2014_EstimatesStats.txt")
str(EstimatesStats)

HarvestEstimatesStats <- sapply(names(EstimatesStats), function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2014[strata.split[1], strata.split[2]],
        EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )
str(HarvestEstimatesStats)

Regional_EstimatesStats <- dget(file = "Estimates objects/KMA2014Strata_Regional_EstimatesStats.txt")
str(Regional_EstimatesStats)

Regional_HarvestEstimatesStats <- sapply(names(EstimatesStats), function(strata) {
  strata.split <- unlist(strsplit(x = strata, split = "_"))
  strata.split <- c(strata.split[1], paste(c(strata.split[2], strata.split[3]), collapse = "_"))
  
  cbind(Regional_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * HarvestByStrata2014[strata.split[1], strata.split[2]],
        Regional_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )
str(Regional_HarvestEstimatesStats)


SheetNames <- sort(names(EstimatesStats))
names(SheetNames) <- SheetNames
mixvec <- SheetNames


# mix <- SheetNames[2]
harvest <- HarvestByStrata2014
dates <- matrix(data = rep("August 25-29", 3), nrow = 3, ncol = 1, 
                dimnames = list(c("SKARLC14", "SUGANC14", "SUYAKC14"),
                                "4_LateLate"))
sampsize <- KMA2014_2015Strata_SampleSizes[, "Final"]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

for(mix in SheetNames) {
  
  # String split by "_
  SheetNames.split <- unlist(strsplit(x = mix, split = "_"))
  
  
  # The first element is the geographic area + year
  geomix <- SheetNames.split[1]
  yr <- unlist(strsplit(x = geomix, split = "C"))[2]
  geo <- unlist(strsplit(x = geomix, split = "1"))[1]
  
  
  # If it is not an annual roll-up, then get the strata number + strata name
  tempmix <- paste(c(SheetNames.split[2], SheetNames.split[3]), collapse = "_")
  
  Caption <- paste("Table X.-Regional and subregional (within Chignik and Kodiak) estimates of stock composition (%) and stock-specific harvest for temporal stratum ",
                   SheetNames.split[2], " (", dates[geomix, tempmix],
                   "; Harvest=", formatC(x = harvest[geomix, tempmix], format = "f", digits = 0, big.mark = ","),
                   "; n=", sampsize[mix], ")", " of ", GeoHeader[geo], ", 20", yr,
                   ". Estimates include median, 90% credibility interval (CI), the probability that the group estimate is equal to zero (P=0), mean, and standard deviation (SD).",
                   sep = '')
  
  Disclaimer <- "Note: Stock composition estimates may not sum to 100% and stock-specific harvest estimates may not sum to the total harvest due to rounding error."
  
  
  TableX <- matrix(data = "", nrow = 24, ncol = 14)
  
  TableX[1, 1] <- Caption
  TableX[2, c(3, 10)] <- c("Stock Composition", "Stock-specific Harvest")
  TableX[3, c(1, 4, 11)] <- c("Reporting Group", rep("90% CI", 2))
  TableX[4, c(1, 2, 3:5, 7:8, 10:14, 6)] <- c("Regional", "Subregional", rep(c("Median", "5%", "95%", "Mean", "SD"), 2), "P=0")
  TableX[5:10, 1] <- c(KMA14GroupsPC[1], "Chignik", "Kodiak", KMA14GroupsPC[12:14])
  TableX[5:10, c(3:5, 7:8)] <- formatC(x = Regional_EstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
  TableX[5:10, 6] <- formatC(x = Regional_EstimatesStats[[mix]][, "P=0"], digits = 2, format = "f")
  TableX[5:10, 10:14] <- formatC(x = Regional_HarvestEstimatesStats[[mix]][, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")
  TableX[11, 12:13] <- c("Total", formatC(x = sum(Regional_HarvestEstimatesStats[[mix]][, "mean"]), digits = 0, format = "f", big.mark = ","))
  TableX[c(13, 16), 1] <- c("Chignik", "Kodiak")
  TableX[c(13:14, 16:23), 2] <- c(SubRegGroups[1:3], "Ayakulik / Frazer", SubRegGroups[5:10])
  TableX[c(13:14, 16:23), c(3:5, 7:8)] <- formatC(x = EstimatesStats[[mix]][SubRegGroups, c("median", "5%", "95%", "mean", "sd")] * 100, digits = 1, format = "f")
  TableX[c(13:14, 16:23), 6] <- formatC(x = EstimatesStats[[mix]][SubRegGroups, "P=0"], digits = 2, format = "f")
  TableX[c(13:14, 16:23), 10:14] <- formatC(x = HarvestEstimatesStats[[mix]][SubRegGroups, c("median", "5%", "95%", "mean", "sd")], digits = 0, format = "f", big.mark = ",")
  TableX[24, 1] <- Disclaimer
  
  write.xlsx(x = as.data.frame(TableX), 
             file = "Estimates tables/KMA Sockeye Estimates Tables Regional.xlsx",
             col.names = FALSE, row.names = FALSE, append = TRUE, sheetName = mix)
}; beep(5)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Double Check Final Numbers ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compare Cook Inlet mean values between 14RG and 6RG estimates

names(KMA2014Strata_EstimatesStats)
names(KMA2014Strata_Regional_EstimatesStats)

sapply(KMA2014Strata[-grep(pattern = "LateLate", x = KMA2014Strata)], function(mix) {
  KMA2014Strata_EstimatesStats[[mix]]["Cook Inlet", "mean"] - KMA2014Strata_Regional_EstimatesStats[[mix]]["Cook Inlet", "mean"]
})


sapply(KMA2015Strata, function(mix) {
  KMA2015Strata_EstimatesStats[[mix]]["Cook Inlet", "mean"] - KMA2015Strata_Regional_EstimatesStats[[mix]]["Cook Inlet", "mean"]
})

sapply(KMA2016Strata, function(mix) {
  KMA2016Strata_EstimatesStats[[mix]]["Cook Inlet", "mean"] - KMA2016Strata_Regional_EstimatesStats[[mix]]["Cook Inlet", "mean"]
})



sapply(KMA2014Strata[-grep(pattern = "LateLate", x = KMA2014Strata)], function(mix) {
  KMA2014Strata_HarvestEstimatesStats[[mix]]["Cook Inlet", "mean"] - KMA2014Strata_Regional_HarvestEstimatesStats[[mix]]["Cook Inlet", "mean"]
})


sapply(KMA2015Strata, function(mix) {
  KMA2015Strata_HarvestEstimatesStats[[mix]]["Cook Inlet", "mean"] - KMA2015Strata_Regional_HarvestEstimatesStats[[mix]]["Cook Inlet", "mean"]
})

sapply(KMA2016Strata, function(mix) {
  KMA2016Strata_HarvestEstimatesStats[[mix]]["Cook Inlet", "mean"] - KMA2016Strata_Regional_HarvestEstimatesStats[[mix]]["Cook Inlet", "mean"]
})



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize to Frazer Reporting Group ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ayakulik_Strata <- grep(pattern = "AYAKC", x = KMA2014_2016Strata, value = TRUE)


Ayakulik_Strata_31RG_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
  maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
  mixvec = Ayakulik_Strata, prior = "",  
  ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE); beep(4)

str(Ayakulik_Strata_31RG_Estimates)


# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Ayakulik_Strata_31RG_Estimates, file = "Estimates objects/Ayakulik_Strata_31RG_Estimates.txt")
dput(Ayakulik_Strata_31RG_Estimates$Stats, file = "Estimates objects/Ayakulik_Strata_31RG_EstimatesStats.txt")

Ayakulik_Strata_31RG_Estimates <- dget(file = "Estimates objects/Ayakulik_Strata_31RG_Estimates.txt")
Ayakulik_Strata_31RG_EstimatesStats <- dget(file = "Estimates objects/Ayakulik_Strata_31RG_EstimatesStats.txt")


sapply(Ayakulik_Strata_31RG_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(Ayakulik_Strata_31RG_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})


QuickBarplot(mixvec = Ayakulik_Strata, estimatesstats = Ayakulik_Strata_31RG_EstimatesStats, groups = KMA31GroupsPC, header = setNames(object = paste("Ayakulik", rep(c(2014:2016), each = 3), c("Early", "Middle", "Late")), nm = Ayakulik_Strata))

# Extrapolate to harvest
HarvestByStrata2014_2016_Final <- rbind(HarvestByStrata2014_Final, HarvestByStrata2015_Final, HarvestByStrata2016_Final)

Ayakulik_Strata_31RG_Harvest_Medians <- sapply(Ayakulik_Strata, function(silly) {
  geo <- unlist(strsplit(x = Ayakulik_Strata[1], split = "_"))[1]
  temp <- paste(unlist(strsplit(x = Ayakulik_Strata[1], split = "_"))[2:3], collapse = "_")
  round(Ayakulik_Strata_31RG_EstimatesStats[[silly]][, "median"] * HarvestByStrata2014_2016_Final[geo, temp])
})
Ayakulik_Strata_31RG_Harvest_Medians["Susitna Yetna", ]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Alitak_Strata <- grep(pattern = "ALITC", x = KMA2014_2016Strata, value = TRUE)


Alitak_Strata_31RG_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
  maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
  mixvec = Alitak_Strata, prior = "",  
  ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE); beep(4)

str(Alitak_Strata_31RG_Estimates)


# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Alitak_Strata_31RG_Estimates, file = "Estimates objects/Alitak_Strata_31RG_Estimates.txt")
dput(Alitak_Strata_31RG_Estimates$Stats, file = "Estimates objects/Alitak_Strata_31RG_EstimatesStats.txt")

Alitak_Strata_31RG_Estimates <- dget(file = "Estimates objects/Alitak_Strata_31RG_Estimates.txt")
Alitak_Strata_31RG_EstimatesStats <- dget(file = "Estimates objects/Alitak_Strata_31RG_EstimatesStats.txt")


sapply(Alitak_Strata_31RG_EstimatesStats, function(Mix) {Mix[, "GR"]})
sapply(Alitak_Strata_31RG_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})


QuickBarplot(mixvec = Alitak_Strata, estimatesstats = Alitak_Strata_31RG_EstimatesStats, groups = KMA31GroupsPC, header = setNames(object = paste("Alitak", rep(c(2014:2016), each = 3), c("Early", "Middle", "Late"))[-1], nm = Alitak_Strata))



PlotPosterior(mixvec = Alitak_Strata, output = Alitak_Strata_31RG_Estimates$Output, 
              groups = KMA31GroupsPC, colors = rep("black", 31), 
              header = setNames(
                object = paste("Alitak", rep(c(2014:2016), each = 3), c("Early", "Middle", "Late"))[-1], 
                nm = Alitak_Strata),
              set.mfrow = c(7, 5), thin = 10)


# Extrapolate to harvest
HarvestByStrata2014_2016_Final <- rbind(HarvestByStrata2014_Final, HarvestByStrata2015_Final, HarvestByStrata2016_Final)

Alitak_Strata_31RG_Harvest_Medians <- sapply(Alitak_Strata, function(silly) {
  geo <- unlist(strsplit(x = Alitak_Strata[1], split = "_"))[1]
  temp <- paste(unlist(strsplit(x = Alitak_Strata[1], split = "_"))[2:3], collapse = "_")
  round(Alitak_Strata_31RG_EstimatesStats[[silly]][, "median"] * HarvestByStrata2014_2016_Final[geo, temp])
})
Alitak_Strata_31RG_Harvest_Medians["Susitna Yetna", ]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA_Strata_31RG_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
  maindir = "BAYES/2014-2016 Mixtures 46loci 14RG/Output", 
  mixvec = KMA2014_2016Strata, prior = "",  
  ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE); beep(4)


dput(KMA_Strata_31RG_Estimates, file = "Estimates objects/KMA_Strata_31RG_Estimates.txt")
dput(KMA_Strata_31RG_Estimates$Stats, file = "Estimates objects/KMA_Strata_31RG_EstimatesStats.txt")

KMA_Strata_31RG_Estimates <- dget(file = "Estimates objects/KMA_Strata_31RG_Estimates.txt")
KMA_Strata_31RG_EstimatesStats <- dget(file = "Estimates objects/KMA_Strata_31RG_EstimatesStats.txt")

# Extrapolate to harvest
HarvestByStrata2014_2016 <- rbind(HarvestByStrata2014, 
                                  cbind(HarvestByStrata2015, "4_LateLate" = rep(NA, 6)), 
                                  cbind(HarvestByStrata2016, "4_LateLate" = rep(NA, 6)))

KMA_Strata_31RG_Harvest_Medians <- sapply(KMA2014_2016Strata, function(silly) {
  geo <- unlist(strsplit(x = silly, split = "_"))[1]
  temp <- paste(unlist(strsplit(x = silly, split = "_"))[2:3], collapse = "_")
  round(KMA_Strata_31RG_EstimatesStats[[silly]][, "median"] * HarvestByStrata2014_2016[geo, temp])
})

KMA_Strata_31RG_Harvest_Means <- sapply(KMA2014_2016Strata, function(silly) {
  geo <- unlist(strsplit(x = silly, split = "_"))[1]
  temp <- paste(unlist(strsplit(x = silly, split = "_"))[2:3], collapse = "_")
  round(KMA_Strata_31RG_EstimatesStats[[silly]][, "mean"] * HarvestByStrata2014_2016[geo, temp])
})


sum(KMA_Strata_31RG_Harvest_Medians["Susitna Yetna", KMA2014Strata])
sum(KMA_Strata_31RG_Harvest_Medians["Susitna Yetna", KMA2015Strata])
sum(KMA_Strata_31RG_Harvest_Medians["Susitna Yetna", KMA2016Strata])


sum(KMA_Strata_31RG_Harvest_Means["Chignik Lake", KMA2014Strata])
sum(KMA_Strata_31RG_Harvest_Means["Chignik Lake", KMA2015Strata])
sum(KMA_Strata_31RG_Harvest_Means["Chignik Lake", KMA2016Strata])


t(sapply(KMA31GroupsPC, function(RG) {
  c("2014" = sum(KMA_Strata_31RG_Harvest_Medians[RG, KMA2014Strata]),
    "2015" = sum(KMA_Strata_31RG_Harvest_Medians[RG, KMA2015Strata]),
    "2016" = sum(KMA_Strata_31RG_Harvest_Medians[RG, KMA2016Strata]))
} ))


t(sapply(KMA31GroupsPC, function(RG) {
  c("2014" = sum(KMA_Strata_31RG_Harvest_Means[RG, KMA2014Strata]),
    "2015" = sum(KMA_Strata_31RG_Harvest_Means[RG, KMA2015Strata]),
    "2016" = sum(KMA_Strata_31RG_Harvest_Means[RG, KMA2016Strata]))
} ))



require(ggplot2)
par(mar = c(6.6, 4.1, 1.1, 1.1))
KMA31Groups_Annual_Barplot <- barplot2(height = sapply(KMA31GroupsPC, function(RG) {
  c("2014" = sum(KMA_Strata_31RG_Harvest_Means[RG, KMA2014Strata]),
    "2015" = sum(KMA_Strata_31RG_Harvest_Means[RG, KMA2015Strata]),
    "2016" = sum(KMA_Strata_31RG_Harvest_Means[RG, KMA2016Strata]))
} ), 
         beside = TRUE, plot.ci = FALSE,
         ylim = c(0, 400000), col = c("darkgreen", "green", "white"), yaxt = "n", xaxt = 'n')
axis(side = 2, at = seq(0, 400000, 50000), labels = formatC(x = seq(0, 400000, 50000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = 1)
legend(legend = 2014:2016, x = "topleft", fill = c("darkgreen", "green", "white"), border = "black", bty = "n", cex = 1.5, title="")
abline(h = 0, xpd = FALSE)
mtext(side = 2, line = 2.5, text = "Number of Fish Harvested (Thousands)", cex = 1.5)
mtext(side = 1, line = 5.5, text = "Reporting Group", cex = 1.5)
text(x = colMeans(KMA31Groups_Annual_Barplot), y = -5000, labels = KMA31GroupsPC, srt = 90, cex = 0.7, xpd = TRUE, adj = 1)

















#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2014 Harvest Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
harvest14 <- read.csv(file = "Harvest/Kodiak Sockeye Salmon Catch by Day and Stat Area 2014.csv", as.is = TRUE)
str(harvest14)

# Convert Date
harvest14$Date.Landed <- as.Date(harvest14$Date.Landed, format = "%Y-%m-%d")
harvest14$Date.Fishing.Began <- as.Date(harvest14$Date.Fishing.Began, format = "%Y-%m-%d")
harvest14$Date.Fishing.Ended <- as.Date(harvest14$Date.Fishing.Ended, format = "%Y-%m-%d")
range(harvest14$Date.Landed)

# Create Temporal Strata
harvest14$Strata <- ifelse(harvest14$Date.Fishing.Began >= as.Date("2014-06-01") & harvest14$Date.Fishing.Began <= as.Date("2014-06-27"),
                           "Early",
                           ifelse(harvest14$Date.Fishing.Began >= as.Date("2014-06-28") & harvest14$Date.Fishing.Began <= as.Date("2014-07-25"),
                                  "Middle",
                                  ifelse(harvest14$Date.Fishing.Began >= as.Date("2014-07-26") & harvest14$Date.Fishing.Began <= as.Date("2014-08-29"),
                                         "Late",
                                         ifelse(harvest14$Date.Fishing.Began > as.Date("2014-08-29"),
                                                "Post",NA))))

harvest14$Strata <- factor(harvest14$Strata, levels = c("Early", "Middle", "Late", "Post"))

# Create Geographic Strata
stat.areas <- unique(harvest14$Stat.Area)

Stat.Area.Uganik <- stat.areas[which(stat.areas >= 25300 & stat.areas <= 25399)]
Stat.Area.Uyak <- stat.areas[which(stat.areas >= 25400 & stat.areas <= 25441)]
Stat.Area.Spiridon <- stat.areas[which(stat.areas == 25450)]
Stat.Area.NorthCape <- stat.areas[which(stat.areas >= 25930 & stat.areas <= 25939)]
Stat.Area.SWNWAfognak <- stat.areas[which(stat.areas >= 25110 & stat.areas <= 25190)]
Stat.Area.NESEAfognak <- stat.areas[which(stat.areas >= 25210 & stat.areas <= 25239)]
Stat.Area.Sitkalidak <- stat.areas[which(stat.areas >= 25800 & stat.areas <= 25899)]
Stat.Area.NEKodiak <- stat.areas[c(which(stat.areas >= 25910 & stat.areas <= 25927), which(stat.areas >= 25940 & stat.areas <= 25946))]
Stat.Area.Alitak <- stat.areas[c(which(stat.areas >= 25710 & stat.areas <= 25720), which(stat.areas >= 25750 & stat.areas <= 25770))]
Stat.Area.MoserOlga <- stat.areas[which(stat.areas >= 25740 & stat.areas <= 25743)]
Stat.Area.Karluk <- stat.areas[c(which(stat.areas >= 25510 & stat.areas <= 25520), which(stat.areas ==25640))]
Stat.Area.Ayakulik <- stat.areas[which(stat.areas >= 25610 & stat.areas <= 25630)]
Stat.Area.NorthShelikof <- stat.areas[which(stat.areas >= 26210 & stat.areas <= 26255)]
Stat.Area.Katmai <- stat.areas[which(stat.areas >= 26260 & stat.areas <= 26270)]
Stat.Area.CapeIgvak <- stat.areas[which(stat.areas >= 26275 & stat.areas <= 26295)]

length(unlist(sapply(objects(pattern = "Stat.Area"), function(x) {get(x)})))
length(stat.areas)

unique(harvest14$Stat.Area)[!stat.areas %in% unlist(sapply(objects(pattern = "Stat.Area"), function(x) {get(x)}))]


harvest14$Geo <- ifelse(harvest14$Stat.Area %in% Stat.Area.Spiridon, "Spiridon",
                        ifelse(harvest14$Stat.Area %in% Stat.Area.Uganik, "Uganik",
                        ifelse(harvest14$Stat.Area %in% Stat.Area.Uyak, "Uyak",
                               ifelse(harvest14$Stat.Area %in% Stat.Area.NorthCape, "NorthCape",
                                      ifelse(harvest14$Stat.Area %in% Stat.Area.SWNWAfognak, "SWNWAfognak",
                                             ifelse(harvest14$Stat.Area %in% Stat.Area.NESEAfognak, "NESEAfognak",
                                                    ifelse(harvest14$Stat.Area %in% Stat.Area.Sitkalidak, "Sitkalidak",
                                                           ifelse(harvest14$Stat.Area %in% Stat.Area.NEKodiak, "NEKodiak",
                                                                  ifelse(harvest14$Stat.Area %in% Stat.Area.Alitak, "Alitak",
                                                                         ifelse(harvest14$Stat.Area %in% Stat.Area.MoserOlga, "MoserOlga",
                                                                                ifelse(harvest14$Stat.Area %in% Stat.Area.Karluk, "Karluk",
                                                                                       ifelse(harvest14$Stat.Area %in% Stat.Area.Ayakulik, "Ayakulik",
                                                                                              ifelse(harvest14$Stat.Area %in% Stat.Area.NorthShelikof, "NorthShelikof",
                                                                                                     ifelse(harvest14$Stat.Area %in% Stat.Area.Katmai, "Katmai",
                                                                                                            ifelse(harvest14$Stat.Area %in% Stat.Area.CapeIgvak, "CapeIgvak", NA
                                                                                                            )))))))))))))))

harvest14$Geo <- factor(harvest14$Geo, levels = c("Spiridon", "Uganik", "Uyak", "NorthCape", "SWNWAfognak", "NESEAfognak", "Sitkalidak", "NEKodiak", "Alitak", "MoserOlga", "Karluk", "Ayakulik", "NorthShelikof", "Katmai", "CapeIgvak"))

table(harvest14$Strata, harvest14$Geo)

require(plyr)
require(reshape)
daply(.data = harvest14, ~Geo+Strata, summarise, harvest = sum(Number))

cast(aggregate(Number ~ Geo + Strata, data = harvest14, sum), Geo ~ Strata, value = "Number")

md <- melt(data = harvest14, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE)
cast(md, Geo~Strata, sum)

aggregate(Number ~ Date.Fishing.Began, data = subset(harvest14, subset = harvest14$Geo == "Uyak" & harvest14$Strata == "Late"), sum)




t(data.matrix(cast(aggregate(Number ~ Geo + Strata, data = harvest14, sum), Geo ~ Strata, value = "Number"))[c(1,2,8,10,11,14), ])

t(data.frame(cast(aggregate(Number ~ Geo + Strata, data = harvest14, sum), Geo ~ Strata, value = "Number"))[c(1,2,8,10,11,14), ])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2015 Harvest Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
harvest15 <- read.csv(file = "Harvest/Kodiak Sockeye Salmon Catch by Day and Stat Area 2015.csv", as.is = TRUE)
str(harvest15)

# Convert Date
harvest15$Date.Landed <- as.Date(harvest15$Date.Landed, format = "%Y-%m-%d")
harvest15$Date.Fishing.Began <- as.Date(harvest15$Date.Fishing.Began, format = "%Y-%m-%d")
harvest15$Date.Fishing.Ended <- as.Date(harvest15$Date.Fishing.Ended, format = "%Y-%m-%d")
range(harvest15$Date.Landed)

# Create Temporal Strata
harvest15$Strata <- ifelse(harvest15$Date.Fishing.Began >= as.Date("2015-06-01") & harvest15$Date.Fishing.Began <= as.Date("2015-07-03"),
                           "Early",
                           ifelse(harvest15$Date.Fishing.Began >= as.Date("2015-07-04") & harvest15$Date.Fishing.Began <= as.Date("2015-08-01"),
                                  "Middle",
                                  ifelse(harvest15$Date.Fishing.Began >= as.Date("2015-08-02") & harvest15$Date.Fishing.Began <= as.Date("2015-08-29"),
                                         "Late",
                                         NA)))

harvest15$Strata <- factor(harvest15$Strata, levels = c("Early", "Middle", "Late"))

# Create Geographic Strata
stat.areas <- unique(harvest15$Stat.Area)

Stat.Area.Uganik <- stat.areas[which(stat.areas >= 25300 & stat.areas <= 25399)]
Stat.Area.Uyak <- stat.areas[which(stat.areas >= 25400 & stat.areas <= 25441)]
Stat.Area.NorthCape <- stat.areas[which(stat.areas >= 25930 & stat.areas <= 25939)]
Stat.Area.SWNWAfognak <- stat.areas[which(stat.areas >= 25110 & stat.areas <= 25190)]
Stat.Area.NESEAfognak <- stat.areas[which(stat.areas >= 25210 & stat.areas <= 25239)]
Stat.Area.Sitkalidak <- stat.areas[which(stat.areas >= 25800 & stat.areas <= 25899)]
Stat.Area.NEKodiak <- stat.areas[c(which(stat.areas >= 25910 & stat.areas <= 25927), which(stat.areas >= 25940 & stat.areas <= 25946))]
Stat.Area.Alitak <- stat.areas[c(which(stat.areas >= 25710 & stat.areas <= 25720), which(stat.areas >= 25750 & stat.areas <= 25770))]
Stat.Area.MoserOlga <- stat.areas[which(stat.areas >= 25740 & stat.areas <= 25743)]
Stat.Area.Karluk <- stat.areas[c(which(stat.areas >= 25510 & stat.areas <= 25520), which(stat.areas ==25640))]
Stat.Area.Ayakulik <- stat.areas[which(stat.areas >= 25610 & stat.areas <= 25630)]
Stat.Area.NorthShelikof <- stat.areas[which(stat.areas >= 26210 & stat.areas <= 26255)]
Stat.Area.Katmai <- stat.areas[which(stat.areas >= 26260 & stat.areas <= 26270)]
Stat.Area.CapeIgvak <- stat.areas[which(stat.areas >= 26275 & stat.areas <= 26295)]

length(unlist(sapply(objects(pattern = "Stat.Area"), function(x) {get(x)})))
length(stat.areas)

unique(harvest15$Stat.Area)[!stat.areas %in% unlist(sapply(objects(pattern = "Stat.Area"), function(x) {get(x)}))]


harvest15$Geo <- ifelse(harvest15$Stat.Area %in% Stat.Area.Uganik, "Uganik",
                        ifelse(harvest15$Stat.Area %in% Stat.Area.Uyak, "Uyak",
                               ifelse(harvest15$Stat.Area %in% Stat.Area.NorthCape, "NorthCape",
                                      ifelse(harvest15$Stat.Area %in% Stat.Area.SWNWAfognak, "SWNWAfognak",
                                             ifelse(harvest15$Stat.Area %in% Stat.Area.NESEAfognak, "NESEAfognak",
                                                    ifelse(harvest15$Stat.Area %in% Stat.Area.Sitkalidak, "Sitkalidak",
                                                           ifelse(harvest15$Stat.Area %in% Stat.Area.NEKodiak, "NEKodiak",
                                                                  ifelse(harvest15$Stat.Area %in% Stat.Area.Alitak, "Alitak",
                                                                         ifelse(harvest15$Stat.Area %in% Stat.Area.MoserOlga, "MoserOlga",
                                                                                ifelse(harvest15$Stat.Area %in% Stat.Area.Karluk, "Karluk",
                                                                                       ifelse(harvest15$Stat.Area %in% Stat.Area.Ayakulik, "Ayakulik",
                                                                                              ifelse(harvest15$Stat.Area %in% Stat.Area.NorthShelikof, "NorthShelikof",
                                                                                                     ifelse(harvest15$Stat.Area %in% Stat.Area.Katmai, "Katmai",
                                                                                                            ifelse(harvest15$Stat.Area %in% Stat.Area.CapeIgvak, "CapeIgvak", NA
                                                                                                            ))))))))))))))

harvest15$Geo <- factor(harvest15$Geo, levels = c("Uganik", "Uyak", "NorthCape", "SWNWAfognak", "NESEAfognak", "Sitkalidak", "NEKodiak", "Alitak", "MoserOlga", "Karluk", "Ayakulik", "NorthShelikof", "Katmai", "CapeIgvak"))

table(harvest15$Strata, harvest15$Geo)

require(plyr)
require(reshape)
daply(.data = harvest15, ~Geo+Strata, summarise, harvest = sum(Number))

cast(aggregate(Number ~ Geo + Strata, data = harvest15, sum), Geo ~ Strata, value = "Number")

md <- melt(data = harvest15, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE)
cast(md, Geo~Strata, sum)

aggregate(Number ~ Date.Fishing.Began, data = subset(harvest15, subset = harvest15$Geo == "Uyak" & harvest15$Strata == "Middle"), sum)
aggregate(Number ~ Fishery.Name, data = subset(harvest15, subset = harvest15$Geo == "Uyak" & harvest15$Strata == "Middle"), sum)
cast(aggregate(Number ~ Fishery.Name + Date.Fishing.Began, data = subset(harvest15, subset = harvest15$Geo == "Uyak" & harvest15$Strata == "Middle"), sum), Date.Fishing.Began ~ Fishery.Name, value = "Number")



aggregate(Number ~ Date.Fishing.Began, data = subset(harvest15, subset = harvest15$Geo == "Uganik" & harvest15$Strata == "Early"), sum)
aggregate(Number ~ Date.Fishing.Began, data = subset(harvest15, subset = harvest15$Geo == "Uganik" & harvest15$Strata == "Late"), sum)
aggregate(Number ~ Date.Fishing.Began, data = subset(harvest15, subset = harvest15$Geo == "Ayakulik" & harvest15$Strata == "Early"), sum)
aggregate(Number ~ Date.Fishing.Began, data = subset(harvest15, subset = harvest15$Geo == "Ayakulik" & harvest15$Strata == "Middle"), sum)
aggregate(Number ~ Date.Fishing.Began, data = subset(harvest15, subset = harvest15$Geo == "Alitak" & harvest15$Strata == "Middle"), sum)
aggregate(Number ~ Date.Fishing.Began, data = subset(harvest15, subset = harvest15$Geo == "Alitak" & harvest15$Strata == "Late"), sum)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2016 Harvest Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
harvest16 <- read.csv(file = "Harvest/Kodiak Sockeye Salmon Catch by Day and Stat Area 2016.csv", as.is = TRUE)
str(harvest16)

# Convert Date
harvest16$Date.Landed <- as.Date(harvest16$Date.Landed, format = "%Y-%m-%d")
harvest16$Date.Fishing.Began <- as.Date(harvest16$Date.Fishing.Began, format = "%Y-%m-%d")
harvest16$Date.Fishing.Ended <- as.Date(harvest16$Date.Fishing.Ended, format = "%Y-%m-%d")
range(harvest16$Date.Landed)

# Create Temporal Strata
harvest16$Strata <- ifelse(harvest16$Date.Fishing.Began >= as.Date("2016-06-01") & harvest16$Date.Fishing.Began <= as.Date("2016-06-27"),
                           "Early",
                           ifelse(harvest16$Date.Fishing.Began >= as.Date("2016-06-28") & harvest16$Date.Fishing.Began <= as.Date("2016-07-25"),
                                  "Middle",
                                  ifelse(harvest16$Date.Fishing.Began >= as.Date("2016-07-26") & harvest16$Date.Fishing.Began <= as.Date("2016-08-28"),
                                         "Late",
                                         NA)))

harvest16$Strata <- factor(harvest16$Strata, levels = c("Early", "Middle", "Late"))

# Create Geographic Strata
stat.areas <- unique(harvest16$Stat.Area)

Stat.Area.Uganik <- stat.areas[which(stat.areas >= 25300 & stat.areas <= 25399)]
Stat.Area.Uyak <- stat.areas[which(stat.areas >= 25400 & stat.areas < 25450)]
Stat.Area.NorthCape <- stat.areas[which(stat.areas >= 25930 & stat.areas <= 25939)]
Stat.Area.SWNWAfognak <- stat.areas[which(stat.areas >= 25110 & stat.areas <= 25190)]
Stat.Area.NESEAfognak <- stat.areas[which(stat.areas >= 25210 & stat.areas <= 25239)]
Stat.Area.Sitkalidak <- stat.areas[which(stat.areas >= 25800 & stat.areas <= 25899)]
Stat.Area.NEKodiak <- stat.areas[c(which(stat.areas >= 25910 & stat.areas <= 25927), which(stat.areas >= 25940 & stat.areas <= 25946))]
Stat.Area.Alitak <- stat.areas[c(which(stat.areas >= 25710 & stat.areas <= 25720), which(stat.areas >= 25750 & stat.areas <= 25770))]
Stat.Area.MoserOlga <- stat.areas[which(stat.areas >= 25740 & stat.areas <= 25743)]
Stat.Area.Karluk <- stat.areas[c(which(stat.areas >= 25510 & stat.areas <= 25520), which(stat.areas ==25640))]
Stat.Area.Ayakulik <- stat.areas[which(stat.areas >= 25610 & stat.areas <= 25630)]
Stat.Area.NorthShelikof <- stat.areas[which(stat.areas >= 26210 & stat.areas <= 26255)]
Stat.Area.Katmai <- stat.areas[which(stat.areas >= 26260 & stat.areas <= 26270)]
Stat.Area.CapeIgvak <- stat.areas[which(stat.areas >= 26275 & stat.areas <= 26295)]

length(unlist(sapply(objects(pattern = "Stat.Area"), function(x) {get(x)})))
length(stat.areas)

unique(harvest16$Stat.Area)[!stat.areas %in% unlist(sapply(objects(pattern = "Stat.Area"), function(x) {get(x)}))]


harvest16$Geo <- ifelse(harvest16$Stat.Area %in% Stat.Area.Uganik, "Uganik",
                        ifelse(harvest16$Stat.Area %in% Stat.Area.Uyak, "Uyak",
                               ifelse(harvest16$Stat.Area %in% Stat.Area.NorthCape, "NorthCape",
                                      ifelse(harvest16$Stat.Area %in% Stat.Area.SWNWAfognak, "SWNWAfognak",
                                             ifelse(harvest16$Stat.Area %in% Stat.Area.NESEAfognak, "NESEAfognak",
                                                    ifelse(harvest16$Stat.Area %in% Stat.Area.Sitkalidak, "Sitkalidak",
                                                           ifelse(harvest16$Stat.Area %in% Stat.Area.NEKodiak, "NEKodiak",
                                                                  ifelse(harvest16$Stat.Area %in% Stat.Area.Alitak, "Alitak",
                                                                         ifelse(harvest16$Stat.Area %in% Stat.Area.MoserOlga, "MoserOlga",
                                                                                ifelse(harvest16$Stat.Area %in% Stat.Area.Karluk, "Karluk",
                                                                                       ifelse(harvest16$Stat.Area %in% Stat.Area.Ayakulik, "Ayakulik",
                                                                                              ifelse(harvest16$Stat.Area %in% Stat.Area.NorthShelikof, "NorthShelikof",
                                                                                                     ifelse(harvest16$Stat.Area %in% Stat.Area.Katmai, "Katmai",
                                                                                                            ifelse(harvest16$Stat.Area %in% Stat.Area.CapeIgvak, "CapeIgvak", NA
                                                                                                            ))))))))))))))

harvest16$Geo <- factor(harvest16$Geo, levels = c("Uganik", "Uyak", "NorthCape", "SWNWAfognak", "NESEAfognak", "Sitkalidak", "NEKodiak", "Alitak", "MoserOlga", "Karluk", "Ayakulik", "NorthShelikof", "Katmai", "CapeIgvak"))

table(harvest16$Strata, harvest16$Geo, useNA = "always")

# Harvest by Spatio-Temporal Strata
require(plyr)
require(reshape)
daply(.data = harvest16, ~Geo+Strata, summarise, harvest = sum(Number))

cast(aggregate(Number ~ Geo + Strata, data = harvest16, sum), Geo ~ Strata, value = "Number")

md <- melt(data = harvest16, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE)
cast(md, Geo~Strata, sum)


sum(aggregate(Number ~ Geo + Strata, data = harvest14, sum)[, "Number"])
sum(aggregate(Number ~ Geo + Strata, data = harvest15, sum)[, "Number"])
sum(aggregate(Number ~ Geo + Strata, data = harvest16, sum)[, "Number"])



# Harvest by Date, By Strata
harvest16.Geo.Day <- cast(aggregate(Number ~ Geo + Date.Fishing.Began, data = harvest16, sum), Date.Fishing.Began ~ Geo, value = "Number")
harvest16.Geo.Day[is.na(harvest16.Geo.Day)] <- 0

write.csv(x = harvest16.Geo.Day, file = "Harvest/Kodiak Sockeye Salmon Catch by Day and Geographic Area 2016.csv", row.names = FALSE)

cast(aggregate(Number ~ Geo + Date.Fishing.Began, data = subset(x = harvest16, subset = Strata == "Early"), sum), Date.Fishing.Began ~ Geo, value = "Number")
cast(aggregate(Number ~ Geo + Date.Fishing.Began, data = subset(x = harvest16, subset = Strata == "Middle"), sum), Date.Fishing.Began ~ Geo, value = "Number")
cast(aggregate(Number ~ Geo + Date.Fishing.Began, data = subset(x = harvest16, subset = Strata == "Late"), sum), Date.Fishing.Began ~ Geo, value = "Number")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Harvest Data Formula ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA_Harvest.f <- function(file, name, yr) {
  
  harvest <- read.csv(file = file, as.is = TRUE)

  # Convert Date
  harvest$Date.Landed <- as.Date(harvest$Date.Landed, format = "%Y-%m-%d")
  harvest$Date.Fishing.Began <- as.Date(harvest$Date.Fishing.Began, format = "%Y-%m-%d")
  harvest$Date.Fishing.Ended <- as.Date(harvest$Date.Fishing.Ended, format = "%Y-%m-%d")
  range(harvest$Date.Landed)
  
  # Create Temporal Strata
  if(yr == 2014) {
    harvest$Strata <- ifelse(harvest$Date.Fishing.Began < as.Date("2014-06-01"),
                             "Pre",
                             ifelse(harvest$Date.Fishing.Began >= as.Date("2014-06-01") & harvest$Date.Fishing.Began <= as.Date("2014-06-27"),
                                    "Early",
                                    ifelse(harvest$Date.Fishing.Began >= as.Date("2014-06-28") & harvest$Date.Fishing.Began <= as.Date("2014-07-25"),
                                           "Middle",
                                           ifelse(harvest$Date.Fishing.Began >= as.Date("2014-07-26") & harvest$Date.Fishing.Began <= as.Date("2014-08-29"),
                                                  "Late",
                                                  ifelse(harvest$Date.Fishing.Began > as.Date("2014-08-29"),
                                                         "Post",
                                                         NA)))))
  }
  
  if(yr == 2015) {
    harvest$Strata <- ifelse(harvest$Date.Fishing.Began < as.Date("2015-06-01"),
                             "Pre",
                             ifelse(harvest$Date.Fishing.Began >= as.Date("2015-06-01") & harvest$Date.Fishing.Began <= as.Date("2015-07-03"),
                                    "Early",
                                    ifelse(harvest$Date.Fishing.Began >= as.Date("2015-07-04") & harvest$Date.Fishing.Began <= as.Date("2015-08-01"),
                                           "Middle",
                                           ifelse(harvest$Date.Fishing.Began >= as.Date("2015-08-02") & harvest$Date.Fishing.Began <= as.Date("2015-08-29"),
                                                  "Late",
                                                  ifelse(harvest$Date.Fishing.Began > as.Date("2015-08-29"),
                                                         "Post",
                                                         NA)))))
  }
  
  if(yr == 2016) {
    harvest$Strata <- ifelse(harvest$Date.Fishing.Began < as.Date("2016-06-01"),
                             "Pre",
                             ifelse(harvest$Date.Fishing.Began >= as.Date("2016-06-01") & harvest$Date.Fishing.Began <= as.Date("2016-06-27"),
                                    "Early",
                                    ifelse(harvest$Date.Fishing.Began >= as.Date("2016-06-28") & harvest$Date.Fishing.Began <= as.Date("2016-07-25"),
                                           "Middle",
                                           ifelse(harvest$Date.Fishing.Began >= as.Date("2016-07-26") & harvest$Date.Fishing.Began <= as.Date("2016-08-29"),
                                                  "Late",
                                                  ifelse(harvest$Date.Fishing.Began > as.Date("2016-08-29"),
                                                         "Post",
                                                         NA)))))
  }
  
  
  harvest$Strata <- factor(harvest$Strata, levels = c("Pre", "Early", "Middle", "Late", "Post"))
  
  # Create Geographic Strata
  stat.areas <- unique(harvest$Stat.Area)
  
  Stat.Area.Uganik <- stat.areas[which(stat.areas >= 25300 & stat.areas <= 25399)]
  Stat.Area.Uyak <- stat.areas[which(stat.areas >= 25400 & stat.areas < 25450)]
  Stat.Area.Spiridon <- stat.areas[which(stat.areas == 25450)]
  Stat.Area.NorthCape <- stat.areas[which(stat.areas >= 25930 & stat.areas <= 25939)]
  Stat.Area.SWNWAfognak <- stat.areas[which(stat.areas >= 25110 & stat.areas <= 25190)]
  Stat.Area.NESEAfognak <- stat.areas[which(stat.areas >= 25210 & stat.areas <= 25239)]
  Stat.Area.Sitkalidak <- stat.areas[which(stat.areas >= 25800 & stat.areas <= 25899)]
  Stat.Area.NEKodiak <- stat.areas[c(which(stat.areas >= 25910 & stat.areas <= 25927), which(stat.areas >= 25940 & stat.areas <= 25946))]
  Stat.Area.Alitak <- stat.areas[c(which(stat.areas >= 25710 & stat.areas <= 25720), which(stat.areas >= 25750 & stat.areas <= 25770))]
  Stat.Area.MoserOlga <- stat.areas[which(stat.areas >= 25740 & stat.areas <= 25743)]
  Stat.Area.Karluk <- stat.areas[c(which(stat.areas >= 25510 & stat.areas <= 25520), which(stat.areas ==25640))]
  Stat.Area.Ayakulik <- stat.areas[which(stat.areas >= 25610 & stat.areas <= 25630)]
  Stat.Area.NorthShelikof <- stat.areas[which(stat.areas >= 26210 & stat.areas <= 26255)]
  Stat.Area.Katmai <- stat.areas[which(stat.areas >= 26260 & stat.areas <= 26270)]
  Stat.Area.CapeIgvak <- stat.areas[which(stat.areas >= 26275 & stat.areas <= 26295)]
  
  length(unlist(sapply(objects(pattern = "Stat.Area"), function(x) {get(x)})))
  length(stat.areas)
  
  unique(harvest$Stat.Area)[!stat.areas %in% unlist(sapply(objects(pattern = "Stat.Area"), function(x) {get(x)}))]
  
  
  harvest$Geo <- ifelse(harvest$Stat.Area %in% Stat.Area.Uganik, "Uganik",
                        ifelse(harvest$Stat.Area %in% Stat.Area.Uyak, "Uyak",
                               ifelse(harvest$Stat.Area %in% Stat.Area.Spiridon, "Spiridon",
                                      ifelse(harvest$Stat.Area %in% Stat.Area.Karluk, "Karluk",
                                             ifelse(harvest$Stat.Area %in% Stat.Area.Ayakulik, "Ayakulik",
                                                    ifelse(harvest$Stat.Area %in% Stat.Area.Alitak, "Alitak",
                                                           ifelse(harvest$Stat.Area %in% Stat.Area.MoserOlga, "MoserOlga",
                                                                  ifelse(harvest$Stat.Area %in% Stat.Area.Sitkalidak, "Sitkalidak",
                                                                         ifelse(harvest$Stat.Area %in% Stat.Area.NEKodiak, "NEKodiak",
                                                                                ifelse(harvest$Stat.Area %in% Stat.Area.NorthCape, "NorthCape",
                                                                                       ifelse(harvest$Stat.Area %in% Stat.Area.NESEAfognak, "NESEAfognak",
                                                                                              ifelse(harvest$Stat.Area %in% Stat.Area.SWNWAfognak, "SWNWAfognak",
                                                                                                     ifelse(harvest$Stat.Area %in% Stat.Area.NorthShelikof, "NorthShelikof",
                                                                                                            ifelse(harvest$Stat.Area %in% Stat.Area.Katmai, "Katmai",
                                                                                                                   ifelse(harvest$Stat.Area %in% Stat.Area.CapeIgvak, "CapeIgvak", NA
                                                                                                                   )))))))))))))))
                        
  harvest$Geo <- factor(harvest$Geo, levels = c("Uganik", "Uyak", "Spiridon", "Karluk", "Ayakulik", "Alitak", "MoserOlga", "Sitkalidak", "NEKodiak", "NorthCape", "NESEAfognak", "SWNWAfognak", "NorthShelikof", "Katmai", "CapeIgvak"))
  
  assign(x = name, value = harvest, pos = 1)
}





KMA_Harvest.f(file = "Harvest/Kodiak Sockeye Salmon Catch by Day and Stat Area 2014.csv", name = "KMASockeyeHarvest_2014", yr = 2014)
KMA_Harvest.f(file = "Harvest/Kodiak Sockeye Salmon Catch by Day and Stat Area 2015.csv", name = "KMASockeyeHarvest_2015", yr = 2015)
KMA_Harvest.f(file = "Harvest/Kodiak Sockeye Salmon Catch by Day and Stat Area 2016.csv", name = "KMASockeyeHarvest_2016", yr = 2016)

KMA_Harvest.f(file = "Harvest/Kodiak Pink Salmon Catch by Day and Stat Area 2014.csv", name = "KMAPinkHarvest_2014", yr = 2014)
KMA_Harvest.f(file = "Harvest/Kodiak Pink Salmon Catch by Day and Stat Area 2015.csv", name = "KMAPinkHarvest_2015", yr = 2015)
KMA_Harvest.f(file = "Harvest/Kodiak Pink Salmon Catch by Day and Stat Area 2016.csv", name = "KMAPinkHarvest_2016", yr = 2016)


require(reshape)
s14 <- cast(melt(data = KMASockeyeHarvest_2014, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE), Geo~Strata, sum)
s15 <- cast(melt(data = KMASockeyeHarvest_2015, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE), Geo~Strata, sum)
s16 <- cast(melt(data = KMASockeyeHarvest_2016, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE), Geo~Strata, sum)

p14 <- cast(melt(data = KMAPinkHarvest_2014, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE), Geo~Strata, sum)
p15 <- cast(melt(data = KMAPinkHarvest_2015, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE), Geo~Strata, sum)
p16 <- cast(melt(data = KMAPinkHarvest_2016, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE), Geo~Strata, sum)

max(p14, p15, p16)
max(s14, s15, s16)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA_Harvest_Simple.f <- function(file, name, yr) {
  
  harvest <- read.csv(file = file, as.is = TRUE)
  
  # Convert Date
  harvest$Date.Landed <- as.Date(harvest$Date.Landed, format = "%Y-%m-%d")
  harvest$Date.Fishing.Began <- as.Date(harvest$Date.Fishing.Began, format = "%Y-%m-%d")
  harvest$Date.Fishing.Ended <- as.Date(harvest$Date.Fishing.Ended, format = "%Y-%m-%d")
  range(harvest$Date.Landed)
  
  # Create Temporal Strata
  if(yr == 2014) {
    harvest$Strata <- ifelse(harvest$Date.Fishing.Began < as.Date("2014-06-01"),
                             "Pre",
                             ifelse(harvest$Date.Fishing.Began >= as.Date("2014-06-01") & harvest$Date.Fishing.Began <= as.Date("2014-06-27"),
                                    "Early",
                                    ifelse(harvest$Date.Fishing.Began >= as.Date("2014-06-28") & harvest$Date.Fishing.Began <= as.Date("2014-07-25"),
                                           "Middle",
                                           ifelse(harvest$Date.Fishing.Began >= as.Date("2014-07-26") & harvest$Date.Fishing.Began <= as.Date("2014-08-29"),
                                                  "Late",
                                                  ifelse(harvest$Date.Fishing.Began > as.Date("2014-08-29"),
                                                         "Post",
                                                         NA)))))
  }
  
  if(yr == 2015) {
    harvest$Strata <- ifelse(harvest$Date.Fishing.Began < as.Date("2015-06-01"),
                             "Pre",
                             ifelse(harvest$Date.Fishing.Began >= as.Date("2015-06-01") & harvest$Date.Fishing.Began <= as.Date("2015-07-03"),
                                    "Early",
                                    ifelse(harvest$Date.Fishing.Began >= as.Date("2015-07-04") & harvest$Date.Fishing.Began <= as.Date("2015-08-01"),
                                           "Middle",
                                           ifelse(harvest$Date.Fishing.Began >= as.Date("2015-08-02") & harvest$Date.Fishing.Began <= as.Date("2015-08-29"),
                                                  "Late",
                                                  ifelse(harvest$Date.Fishing.Began > as.Date("2015-08-29"),
                                                         "Post",
                                                         NA)))))
  }
  
  if(yr == 2016) {
    harvest$Strata <- ifelse(harvest$Date.Fishing.Began < as.Date("2016-06-01"),
                             "Pre",
                             ifelse(harvest$Date.Fishing.Began >= as.Date("2016-06-01") & harvest$Date.Fishing.Began <= as.Date("2016-06-27"),
                                    "Early",
                                    ifelse(harvest$Date.Fishing.Began >= as.Date("2016-06-28") & harvest$Date.Fishing.Began <= as.Date("2016-07-25"),
                                           "Middle",
                                           ifelse(harvest$Date.Fishing.Began >= as.Date("2016-07-26") & harvest$Date.Fishing.Began <= as.Date("2016-08-29"),
                                                  "Late",
                                                  ifelse(harvest$Date.Fishing.Began > as.Date("2016-08-29"),
                                                         "Post",
                                                         NA)))))
  }
  
  
  harvest$Strata <- factor(harvest$Strata, levels = c("Pre", "Early", "Middle", "Late", "Post"))
  
  # Create Geographic Strata
  stat.areas <- unique(harvest$Stat.Area)
  
  Stat.Area.Uganik <- stat.areas[which(stat.areas >= 25300 & stat.areas <= 25399)]
  Stat.Area.Uyak <- stat.areas[which(stat.areas >= 25400 & stat.areas < 25450)]
  Stat.Area.Alitak <- stat.areas[c(which(stat.areas >= 25710 & stat.areas <= 25720), which(stat.areas >= 25750 & stat.areas <= 25770))]
  Stat.Area.Karluk <- stat.areas[c(which(stat.areas >= 25510 & stat.areas <= 25520), which(stat.areas ==25640))]
  Stat.Area.Ayakulik <- stat.areas[which(stat.areas >= 25610 & stat.areas <= 25630)]
  Stat.Area.CapeIgvak <- stat.areas[which(stat.areas >= 26275 & stat.areas <= 26295)]
  
  length(unlist(sapply(objects(pattern = "Stat.Area"), function(x) {get(x)})))
  length(stat.areas)
  
  unique(harvest$Stat.Area)[!stat.areas %in% unlist(sapply(objects(pattern = "Stat.Area"), function(x) {get(x)}))]
  
  
  harvest$Geo <- ifelse(harvest$Stat.Area %in% Stat.Area.Uganik, "Uganik",
                        ifelse(harvest$Stat.Area %in% Stat.Area.Uyak, "Uyak",
                               ifelse(harvest$Stat.Area %in% Stat.Area.Karluk, "Karluk",
                                      ifelse(harvest$Stat.Area %in% Stat.Area.Ayakulik, "Ayakulik",
                                             ifelse(harvest$Stat.Area %in% Stat.Area.Alitak, "Alitak",
                                                    ifelse(harvest$Stat.Area %in% Stat.Area.CapeIgvak, "Igvak", "Unsampled"
                                                    ))))))
  
  harvest$Geo <- factor(harvest$Geo, levels = rev(c("Uganik", "Uyak", "Karluk", "Ayakulik", "Alitak", "Igvak", "Unsampled")))
  
  assign(x = name, value = harvest, pos = 1)
}

KMA_Harvest_Simple.f(file = "Harvest/Kodiak Sockeye Salmon Catch by Day and Stat Area 2014.csv", name = "KMASockeyeHarvestSimple_2014", yr = 2014)
KMA_Harvest_Simple.f(file = "Harvest/Kodiak Sockeye Salmon Catch by Day and Stat Area 2015.csv", name = "KMASockeyeHarvestSimple_2015", yr = 2015)
KMA_Harvest_Simple.f(file = "Harvest/Kodiak Sockeye Salmon Catch by Day and Stat Area 2016.csv", name = "KMASockeyeHarvestSimple_2016", yr = 2016)

require(reshape)
s14 <- cast(melt(data = KMASockeyeHarvestSimple_2014, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE), Geo~Strata, sum)
s15 <- cast(melt(data = KMASockeyeHarvestSimple_2015, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE), Geo~Strata, sum)
s16 <- cast(melt(data = KMASockeyeHarvestSimple_2016, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE), Geo~Strata, sum)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Plot_KMA_Harvest.f <- function(dat, species, yr, maxdat = NULL) {
  if(is.null(maxdat)) {maxdat <- max(dat)}
  require(lattice)
  new.colors <- colorRampPalette(c("white", "black"))
  levelplot(t(apply(t(dat), 1, rev)), col.regions = new.colors, xlab = "Temporal Strata", 
            ylab = "Geographic Area", main = paste(species, yr), at = seq(0, maxdat, length.out = 100), 
            scales = list(x = list(rot = 45)),
            panel = function(...) {
              panel.levelplot(...)
            })
}

Plot_KMA_Harvest.f(dat = s14, species = "Sockeye", yr = 2014, maxdat = 5e5)
Plot_KMA_Harvest.f(dat = s15, species = "Sockeye", yr = 2015, maxdat = 5e5)
Plot_KMA_Harvest.f(dat = s16, species = "Sockeye", yr = 2016, maxdat = 5e5)

Plot_KMA_Harvest.f(dat = p14, species = "Pink", yr = 2014, maxdat = 6e6)
Plot_KMA_Harvest.f(dat = p15, species = "Pink", yr = 2015, maxdat = 6e6)
Plot_KMA_Harvest.f(dat = p16, species = "Pink", yr = 2016, maxdat = 6e6)


Plot_KMA_Harvest.f(dat = s14, species = "Sockeye", yr = 2014, maxdat = 5.2e5)
Plot_KMA_Harvest.f(dat = p14, species = "Pink", yr = 2014, maxdat = 5.2e5)


Plot_KMA_StatArea_Harvest.f <- function(species, yr, geos, maxdat = NULL){
  dat <- cast(melt(data = subset(x = get(paste0("KMA", species, "Harvest_", yr)), subset = Geo %in% geos), id.vars = c("Strata", "Stat.Area"), measure.vars = "Number", na.rm = TRUE), Stat.Area~Strata, sum)
  if(is.null(maxdat)) {maxdat <- max(dat)}
  require(lattice)
  new.colors <- colorRampPalette(c("white", "black"))
  levelplot(t(apply(t(dat), 1, rev)), col.regions = new.colors, xlab = "Temporal Strata", 
            ylab = "Geographic Area", main = paste(species, yr), at = seq(0, maxdat, length.out = 100), 
            scales = list(x = list(rot = 45)),
            panel = function(...) {
              panel.levelplot(...)
            })
}



Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2014, geos = c("Uganik"), maxdat = 1.2e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2015, geos = c("Uganik"), maxdat = 1.2e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2016, geos = c("Uganik"), maxdat = 1.2e5)

Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2014, geos = c("Uyak", "Spiridon"), maxdat = 1.2e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2015, geos = c("Uyak", "Spiridon"), maxdat = 1.2e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2016, geos = c("Uyak", "Spiridon"), maxdat = 1.2e5)

Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2014, geos = c("Karluk"), maxdat = 3.5e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2015, geos = c("Karluk"), maxdat = 3.5e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2016, geos = c("Karluk"), maxdat = 3.5e5)

Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2014, geos = c("Ayakulik"), maxdat = 3.5e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2015, geos = c("Ayakulik"), maxdat = 3.5e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2016, geos = c("Ayakulik"), maxdat = 3.5e5)

Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2014, geos = c("Alitak", "MoserOlga"), maxdat = 1.2e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2015, geos = c("Alitak", "MoserOlga"), maxdat = 1.2e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2016, geos = c("Alitak", "MoserOlga"), maxdat = 1.2e5)

Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2014, geos = c("NorthShelikof", "Katmai", "CapeIgvak"), maxdat = 1.2e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2015, geos = c("NorthShelikof", "Katmai", "CapeIgvak"), maxdat = 1.2e5)
Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2016, geos = c("NorthShelikof", "Katmai", "CapeIgvak"), maxdat = 1.2e5)


Plot_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2015, geos = c("Alitak", "MoserOlga"))
Plot_KMA_StatArea_Harvest.f(species = "Pink", yr = 2015, geos = c("Alitak", "MoserOlga"))



Table_KMA_StatArea_Harvest.f <- function(species, yr){
  dat <- as.matrix(cast(melt(data = get(paste0("KMA", species, "Harvest_", yr)), id.vars = c("Strata", "Stat.Area"), measure.vars = "Number", na.rm = TRUE), Stat.Area~Strata, sum))
  dat
}


s14stat <- Table_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2014)
s15stat <- Table_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2015)
s16stat <- Table_KMA_StatArea_Harvest.f(species = "Sockeye", yr = 2016)
p14stat <- Table_KMA_StatArea_Harvest.f(species = "Pink", yr = 2014)
p15stat <- Table_KMA_StatArea_Harvest.f(species = "Pink", yr = 2015)
p16stat <- Table_KMA_StatArea_Harvest.f(species = "Pink", yr = 2016)

sapply(list(s14stat, s15stat, s16stat, p14stat, p15stat, p16stat), function(x) {dim(x)[2]})

KMA_StatAreas <- as.character(read.csv(file = "Figures/Maps/KMAStatAreas.csv")[,1])

KMA_Sockeye_Pink_Harvest_Strata.mat <- matrix(data = NA, 
                                              nrow = length(KMA_StatAreas), 
                                              ncol = sum(sapply(list(s14stat, s15stat, s16stat, p14stat, p15stat, p16stat), function(x) {dim(x)[2]})), 
                                              dimnames = list(KMA_StatAreas,
                                                              unlist(sapply(c(paste0("s", 14:16), paste0("p", 14:16)), function(x) {
                                                                y <- dimnames(get(paste0(x, "stat")))[2]
                                                                sapply(y, function(i) {paste(x, i, sep = "_")})
                                                              }))
                                              ))


KMA_Sockeye_Pink_Harvest_Strata.mat <- 
  Reduce(f = cbind, x = 
           sapply(list(s14stat, s15stat, s16stat, p14stat, p15stat, p16stat), function(y) {
             apply(y, 2, function(strata) {
               x <- rep(0, length(KMA_StatAreas))
               x[which(KMA_StatAreas %in% names(strata))] <- strata
               x
             })
           })
  )

dimnames(KMA_Sockeye_Pink_Harvest_Strata.mat) <- 
  list(KMA_StatAreas,
       unlist(sapply(c(paste0("s", 14:16), paste0("p", 14:16)), function(x) {
         y <- dimnames(get(paste0(x, "stat")))[2]
         sapply(y, function(i) {paste(x, i, sep = "_")})
       })))

str(KMA_Sockeye_Pink_Harvest_Strata.mat)

write.csv(x = KMA_Sockeye_Pink_Harvest_Strata.mat,
          file = "Figures/Maps/KMA_Sockeye_Pink_Harvest_Strata.csv")

max(KMA_Sockeye_Pink_Harvest_Strata.mat)

which.max(apply(KMA_Sockeye_Pink_Harvest_Strata.mat, 2, max))




Alitak_Pink_Harvest.f <- function(yr) {
  p <- cast(melt(data = subset(x = get(paste0("KMAPinkHarvest_", yr)), subset = Geo == "Alitak" | Geo == "MoserOlga"), id.vars = c("Strata", "Stat.Area"), measure.vars = "Number", na.rm = TRUE), Stat.Area~Strata, sum)
  
  require(lattice)
  new.colors <- colorRampPalette(c("white", "black"))
  levelplot(t(as.matrix(p)), col.regions = new.colors, xlab = "Temporal Strata", 
            ylab = "Stat Area", main = paste0("Pink ", yr), at = seq(0, max(p), length.out = 100), 
            scales = list(x = list(rot = 45)),
            panel = function(...) {
              panel.levelplot(...)
            })
}


Alitak_Sockeye_Harvest.f <- function(yr) {
  s <- cast(melt(data = subset(x = get(paste0("KMASockeyeHarvest_", yr)), subset = Geo == "Alitak" | Geo == "MoserOlga"), id.vars = c("Strata", "Stat.Area"), measure.vars = "Number", na.rm = TRUE), Stat.Area~Strata, sum)

  require(lattice)
  new.colors <- colorRampPalette(c("white", "black"))
  levelplot(t(as.matrix(s)), col.regions = new.colors, xlab = "Temporal Strata", 
            ylab = "Stat Area", main = paste0("Sockeye ", yr), at = seq(0, max(s), length.out = 100), 
            scales = list(x = list(rot = 45)),
            panel = function(...) {
              panel.levelplot(...)
            })
}



Alitak_Pink_Harvest.f(yr = 2014)
Alitak_Sockeye_Harvest.f(yr = 2014)

Alitak_Pink_Harvest.f(yr = 2015)
Alitak_Sockeye_Harvest.f(yr = 2015)

Alitak_Pink_Harvest.f(yr = 2016)
Alitak_Sockeye_Harvest.f(yr = 2016)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Maps ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
while(!require(maps)) {install.packages("maps")}
while(!require(mapdata)) {install.packages("mapdata")}
while(!require(maptools)) {install.packages("maptools")}
while(!require(GISTools)) {install.packages("GISTools")}
while(!require(rgeos)) {install.packages("rgeos")}
while(!require(sp)) {install.packages("sp")}
require(devEMF)

StatAreas.shp <- readShapePoly("Figures/Maps/pvs_stat_KMA.shp")
str(StatAreas.shp, max.level = 2)
str(StatAreas.shp@data)

str(StatAreas.shp[!is.na(StatAreas.shp@data$Statarea)], max.level = 2)  # fail

KMA_StatAreas <- as.character(read.csv(file = "Figures/Maps/KMAStatAreas.csv")[,1])
KMAStatAreas.shp <- subset(StatAreas.shp, StatAreas.shp@data$Statarea %in% KMA_StatAreas)
str(KMAStatAreas.shp)
max(KMAStatAreas.shp@data[, 9:21])  # Sockeye max stat area
max(KMAStatAreas.shp@data[, 22:33])  # Sockeye max stat area

which.max(apply(KMAStatAreas.shp@data[, 9:21], 2, max))


colorRampPalette(c("white", "black"))(100)
colorRampPalette(c("white", "black"))(100)[round(StatAreas.shp@data$s14_Early / 3.3e3) + 1]

map("worldHires", "usa", xlim = c(-156.6, -151.7), ylim = c(56.3, 59), col = "gray90", fill = TRUE)
plot(KMAStatAreas.shp, add = TRUE, col = "white", border = TRUE)
plot(KMAStatAreas.shp, add = TRUE, 
     col = colorRampPalette(c("white", "black"))(100)[round(KMAStatAreas.shp@data$s15_Middle / 3.3e3) + 1],
     border = TRUE)

plot(KMAStatAreas.shp, add = TRUE, 
     col = colorRampPalette(c("white", "black"))(100)[round(KMAStatAreas.shp@data[, "s14_Middle"] / 3.3e3) + 1],
     border = TRUE)


Plot_KMA_Harvest_Map.f <- function(area, max.col = NULL) {
  if(is.null(max.col)) {max.col <- max(KMAStatAreas.shp@data[, area])}
  
  map("worldHires", "usa", xlim = c(-156.6, -151.7), ylim = c(56.3, 59), col = "gray90", fill = TRUE)
  plot(KMAStatAreas.shp, add = TRUE, 
       col = colorRampPalette(c("white", "black"))(101)[round(KMAStatAreas.shp@data[, area] / (max.col/100)) + 1],
       border = TRUE)
}


setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Figures")
# dir.create("Harvest Maps")
setwd("Harvest Maps")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Uganik <- KMAStatAreas.shp[which(KMAStatAreas.shp@data$Statarea %in% seq(from = 25300, to = 25399, by = 1)), ]
Uganik@data$SampArea <- "Uganik"
Uganik.dis <- gUnaryUnion(Uganik, id = Uganik@data$SampArea)
Uyak <- KMAStatAreas.shp[which(KMAStatAreas.shp@data$Statarea %in% seq(from = 25400, to = 25499, by = 1)), ]
Uyak@data$SampArea <- "Uyak"
Uyak.dis <- gUnaryUnion(Uyak, id = Uyak@data$SampArea)
Karluk <- KMAStatAreas.shp[which(KMAStatAreas.shp@data$Statarea %in% c("25510", "25520", "25640")), ]
Karluk@data$SampArea <- "Karluk"
Karluk.dis <- gUnaryUnion(Karluk, id = Karluk@data$SampArea)
Ayakulik <- KMAStatAreas.shp[which(KMAStatAreas.shp@data$Statarea %in% seq(from = 25610, to = 25630, by = 1)), ]
Ayakulik@data$SampArea <- "Ayakulik"
Ayakulik.dis <- gUnaryUnion(Ayakulik, id = Ayakulik@data$SampArea)
Alitak <- KMAStatAreas.shp[which(KMAStatAreas.shp@data$Statarea %in% c("25710", "25720", "25750", "25760", "25770")), ]
Alitak@data$SampArea <- "Alitak"
Alitak.dis <- gUnaryUnion(Alitak, id = Alitak@data$SampArea)
Igvak <- KMAStatAreas.shp[which(KMAStatAreas.shp@data$Statarea %in% c("26275", "26280", "26290", "26295")), ]
Igvak@data$SampArea <- "Igvak"
Igvak.dis <- gUnaryUnion(Igvak, id = Igvak@data$SampArea)


Plot_KMA_Harvest_Map.f <- function(area, max.col = NULL, cex.lab = 0.9) {
  # emf(file = paste0("KMA Harvest by Stat Area ", area, ".emf"), width = 5.75, height = 5.75, family = "serif", bg = "white")
  if(is.null(max.col)) {max.col <- max(KMAStatAreas.shp@data[, area])}
  map("worldHires", "usa", xlim = c(-156.55, -151.71), ylim = c(56.36, 58.85), col = "gray90", fill = TRUE)
  color.ramp <- colorRampPalette(c("white", "green", "darkgreen", "black"))(101)
  plot(KMAStatAreas.shp, add = TRUE, 
       col = color.ramp[round(KMAStatAreas.shp@data[, area] / (max.col/100)) + 1],
       border = TRUE)
  legend("bottomleft", 
         legend = c(max.col, rep("", 99), 0),
         fill = rev(color.ramp), 
         border = NA,
         bty = 'n', x.intersp = 0.5, y.intersp = 0.07, lty = NULL)
  text(x = -156, y = 56.7, labels = "Harvest", cex = 1.3)
  
  yr <- paste0("20", paste0(unlist(strsplit(x = area, split = ""))[2:3], collapse = ""), collapse = "")
  strata <- unlist(strsplit(x = area, split = "_"))[2]
  
  text(x = -152.8, y = 56.70, labels = paste(yr, strata), cex = 1.3)
  maps::map.scale(x = -153.6, y = 56.46, ratio = FALSE, relwidth = 0.2)
  north.arrow(xb = -152, yb = 56.6, len = 0.05, lab = "N")
  
  plot(Uganik.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -153.7, y = 58.04, labels = "Uganik\nKupreanof", cex = cex.lab, adj = c(1, 0.5))
  plot(Uyak.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -154.15, y = 57.8, labels = "Uyak", cex = cex.lab, adj = c(1, 0.5))
  plot(Karluk.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -154.72, y = 57.62, labels = "Karluk\nSturgeon", cex = cex.lab, adj = c(1, 0.5))
  plot(Ayakulik.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -154.95, y = 57.3, labels = "Ayakulik\nHalibut Bay", cex = cex.lab, adj = c(1, 0.5))
  plot(Alitak.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -154.43, y = 56.84, labels = "Alitak", cex = cex.lab, adj = c(1, 0.5))
  plot(Igvak.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -155.85, y = 57.45, labels = "Igvak", cex = cex.lab, adj = c(0, 0.5))
    # dev.off()
  }

# What is sockeye max?
max(KMAStatAreas.shp@data[, 9:21])
which.max(apply(KMAStatAreas.shp@data[, 9:21], 2, max))

Plot_KMA_Harvest_Map.f(area = "s14_Early", max.col = 3.3e5)
Plot_KMA_Harvest_Map.f(area = "s15_Early", max.col = 3.3e5)
Plot_KMA_Harvest_Map.f(area = "s16_Early", max.col = 3.3e5)

Plot_KMA_Harvest_Map.f(area = "s14_Middle", max.col = 3.3e5)
Plot_KMA_Harvest_Map.f(area = "s15_Middle", max.col = 3.3e5)
Plot_KMA_Harvest_Map.f(area = "s16_Middle", max.col = 3.3e5)

Plot_KMA_Harvest_Map.f(area = "s14_Late", max.col = 3.3e5)
Plot_KMA_Harvest_Map.f(area = "s15_Late", max.col = 3.3e5)
Plot_KMA_Harvest_Map.f(area = "s16_Late", max.col = 3.3e5)

Plot_KMA_Harvest_Map.f(area = "s14_Post", max.col = 3.3e5)
Plot_KMA_Harvest_Map.f(area = "s15_Post", max.col = 3.3e5)
Plot_KMA_Harvest_Map.f(area = "s16_Post", max.col = 3.3e5)


# What is pink max?
Plot_KMA_Harvest_Map.f <- function(area, max.col = NULL, cex.lab = 0.9) {
  # emf(file = paste0("KMA Harvest by Stat Area ", area, ".emf"), width = 5.75, height = 5.75, family = "serif", bg = "white")
  if(is.null(max.col)) {max.col <- max(KMAStatAreas.shp@data[, area])}
  map("worldHires", "usa", xlim = c(-156.55, -151.71), ylim = c(56.36, 58.85), col = "gray90", fill = TRUE)
  color.ramp <- c(colorRampPalette(c("white", "green", "darkgreen", "black"))(11), colorRampPalette(c("lightpink", "red", "darkred", "purple", "purple4"))(90))
  plot(KMAStatAreas.shp, add = TRUE, 
       col = color.ramp[round(KMAStatAreas.shp@data[, area] / (max.col/100)) + 1],
       border = TRUE)
  legend("bottomleft", 
         legend = c(max.col, rep("", 89), 3.3e5, rep("", 9), 0),
         fill = rev(color.ramp), 
         border = NA,
         bty = 'n', x.intersp = 0.5, y.intersp = 0.07, lty = NULL)
  text(x = -156, y = 56.7, labels = "Harvest", cex = 1.3)
  
  yr <- paste0("20", paste0(unlist(strsplit(x = area, split = ""))[2:3], collapse = ""), collapse = "")
  strata <- unlist(strsplit(x = area, split = "_"))[2]
  
  text(x = -152.8, y = 56.70, labels = paste(yr, strata), cex = 1.3)
  maps::map.scale(x = -153.6, y = 56.46, ratio = FALSE, relwidth = 0.2)
  north.arrow(xb = -152, yb = 56.6, len = 0.05, lab = "N")
  
  plot(Uganik.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -153.7, y = 58.04, labels = "Uganik\nKupreanof", cex = cex.lab, adj = c(1, 0.5))
  plot(Uyak.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -154.15, y = 57.8, labels = "Uyak", cex = cex.lab, adj = c(1, 0.5))
  plot(Karluk.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -154.72, y = 57.62, labels = "Karluk\nSturgeon", cex = cex.lab, adj = c(1, 0.5))
  plot(Ayakulik.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -154.95, y = 57.3, labels = "Ayakulik\nHalibut Bay", cex = cex.lab, adj = c(1, 0.5))
  plot(Alitak.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -154.43, y = 56.84, labels = "Alitak", cex = cex.lab, adj = c(1, 0.5))
  plot(Igvak.dis, add = TRUE, border = "black", lwd = 3)
  text(x = -155.85, y = 57.45, labels = "Igvak", cex = cex.lab, adj = c(0, 0.5))
  # dev.off()
}

max(KMAStatAreas.shp@data[, 22:33])
which.max(apply(KMAStatAreas.shp@data[, 22:33], 2, max))

Plot_KMA_Harvest_Map.f(area = "p14_Early", max.col = 3.3e6)
Plot_KMA_Harvest_Map.f(area = "p15_Early", max.col = 3.3e6)
Plot_KMA_Harvest_Map.f(area = "p16_Early", max.col = 3.3e6)

Plot_KMA_Harvest_Map.f(area = "p14_Middle", max.col = 3.3e6)
Plot_KMA_Harvest_Map.f(area = "p15_Middle", max.col = 3.3e6)
Plot_KMA_Harvest_Map.f(area = "p16_Middle", max.col = 3.3e6)

Plot_KMA_Harvest_Map.f(area = "p14_Late", max.col = 3.3e6)
Plot_KMA_Harvest_Map.f(area = "p15_Late", max.col = 3.3e6)
Plot_KMA_Harvest_Map.f(area = "p16_Late", max.col = 3.3e6)

Plot_KMA_Harvest_Map.f(area = "p14_Post", max.col = 3.3e6)
Plot_KMA_Harvest_Map.f(area = "p15_Post", max.col = 3.3e6)
Plot_KMA_Harvest_Map.f(area = "p16_Post", max.col = 3.3e6)



# Lattice Plot
Plot_KMA_Harvest.f(dat = s14, species = "Sockeye", yr = 2014, maxdat = 5e5)
Plot_KMA_Harvest.f(dat = p14, species = "Pink", yr = 2014)
Plot_KMA_Harvest.f(dat = p14, species = "Pink", yr = 2014, maxdat = 1e6)
Plot_KMA_Harvest_Map.f(area = "s14_Early")
Plot_KMA_Harvest_Map.f(area = "p14_Early")
Plot_KMA_Harvest_Map.f(area = "s14_Middle")
Plot_KMA_Harvest_Map.f(area = "p14_Middle")
Plot_KMA_Harvest_Map.f(area = "s14_Late")
Plot_KMA_Harvest_Map.f(area = "p14_Late")

Plot_KMA_Harvest.f(dat = s15, species = "Sockeye", yr = 2015, maxdat = 5e5)
Plot_KMA_Harvest.f(dat = p15, species = "Pink", yr = 2015)
Plot_KMA_Harvest.f(dat = p15, species = "Pink", yr = 2015, maxdat = 1e6)
Plot_KMA_Harvest_Map.f(area = "s15_Early")
Plot_KMA_Harvest_Map.f(area = "p15_Early")
Plot_KMA_Harvest_Map.f(area = "s15_Middle")
Plot_KMA_Harvest_Map.f(area = "p15_Middle")
Plot_KMA_Harvest_Map.f(area = "s15_Late")
Plot_KMA_Harvest_Map.f(area = "p15_Late")


Plot_KMA_Harvest.f(dat = s16, species = "Sockeye", yr = 2016, maxdat = 5e5)
Plot_KMA_Harvest.f(dat = p16, species = "Pink", yr = 2016)
Plot_KMA_Harvest.f(dat = p16, species = "Pink", yr = 2016, maxdat = 1e6)
Plot_KMA_Harvest_Map.f(area = "s16_Early")
Plot_KMA_Harvest_Map.f(area = "p16_Early")
Plot_KMA_Harvest_Map.f(area = "s16_Middle")
Plot_KMA_Harvest_Map.f(area = "p16_Middle")
Plot_KMA_Harvest_Map.f(area = "s16_Late")
Plot_KMA_Harvest_Map.f(area = "p16_Late")



Plot_KMA_Harvest.f(dat = p14, species = "Pink", yr = 2014, maxdat = 6e6)
Plot_KMA_Harvest.f(dat = p15, species = "Pink", yr = 2015, maxdat = 6e6)
Plot_KMA_Harvest.f(dat = p16, species = "Pink", yr = 2016, maxdat = 6e6)





# What is pink max?
max(KMAStatAreas.shp@data[, 22:33])
which.max(apply(KMAStatAreas.shp@data[, 22:33], 2, max))

Plot_KMA_Harvest_Map.f(area = "p14_Early", max.col = 3.1e6)
Plot_KMA_Harvest_Map.f(area = "p15_Early", max.col = 3.1e6)
Plot_KMA_Harvest_Map.f(area = "p16_Early", max.col = 3.1e6)

Plot_KMA_Harvest_Map.f(area = "p14_Middle", max.col = 3.1e6)
Plot_KMA_Harvest_Map.f(area = "p15_Middle", max.col = 3.1e6)
Plot_KMA_Harvest_Map.f(area = "p16_Middle", max.col = 3.1e6)

Plot_KMA_Harvest_Map.f(area = "p14_Late", max.col = 3.1e6)
Plot_KMA_Harvest_Map.f(area = "p15_Late", max.col = 3.1e6)
Plot_KMA_Harvest_Map.f(area = "p16_Late", max.col = 3.1e6)


# Compare sockeye and pink
Plot_KMA_Harvest_Map.f(area = "s15_Middle")
Plot_KMA_Harvest_Map.f(area = "p15_Middle")
legend("bottomright", legend = c(0, 10000), fill = c("white", "black"), bty = 'n')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Harvest for Table 1 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

md <- melt(data = harvest14, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE)
x2014 <- t(data.frame(cast(md, Geo~Strata, sum))[c(1,2,10,11,8,14), ])
write.xlsx(x = x2014, file = "V:/Documents/4_Westward/Sockeye/KMA 2014-2016/KMA Mixtures FMS/Table 1.xlsx",
           col.names = FALSE, sheetName = "2014 Harvest", append = TRUE)

md <- melt(data = harvest15, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE)
x2015 <- t(data.frame(cast(md, Geo~Strata, sum))[c(1,2,10,11,8,14), ])
write.xlsx(x = x2015, file = "V:/Documents/4_Westward/Sockeye/KMA 2014-2016/KMA Mixtures FMS/Table 1.xlsx",
           col.names = FALSE, sheetName = "2015 Harvest", append = TRUE)

md <- melt(data = harvest16, id.vars = c("Strata", "Geo"), measure.vars = "Number", na.rm = TRUE)
x2016 <- t(data.frame(cast(md, Geo~Strata, sum))[c(1,2,10,11,8,14), ])
write.xlsx(x = x2016, file = "V:/Documents/4_Westward/Sockeye/KMA 2014-2016/KMA Mixtures FMS/Table 1.xlsx",
           col.names = FALSE, sheetName = "2016 Harvest", append = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samples14 <- read.csv(file = "Harvest/Bulk Tissue Inventory 2014.csv", as.is = TRUE)
str(samples14)

# Convert Date
samples14$Bulk.Date <- as.Date(samples14$Bulk.Date, format = "%Y-%m-%d")
samples14$Bulk.End.Date <- as.Date(samples14$Bulk.End.Date, format = "%Y-%m-%d")
range(samples14$Bulk.Date)

# Create Temporal Strata
samples14$Strata <- ifelse(samples14$Bulk.Date >= as.Date("2014-06-01") & samples14$Bulk.Date <= as.Date("2014-06-27"),
                           "Early",
                           ifelse(samples14$Bulk.Date >= as.Date("2014-06-28") & samples14$Bulk.Date <= as.Date("2014-07-25"),
                                  "Middle",
                                  ifelse(samples14$Bulk.Date >= as.Date("2014-07-26") & samples14$Bulk.Date <= as.Date("2014-08-29"),
                                         "Late",
                                         NA)))

samples14$Strata <- factor(samples14$Strata, levels = c("Early", "Middle", "Late"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samples15 <- read.csv(file = "Harvest/Bulk Tissue Inventory 2015.csv", as.is = TRUE)
str(samples15)

# Convert Date
samples15$Bulk.Date <- as.Date(samples15$Bulk.Date, format = "%Y-%m-%d")
samples15$Bulk.End.Date <- as.Date(samples15$Bulk.End.Date, format = "%Y-%m-%d")
range(samples15$Bulk.Date)

# Create Temporal Strata
samples15$Strata <- ifelse(samples15$Bulk.Date >= as.Date("2015-06-01") & samples15$Bulk.Date <= as.Date("2015-07-03"),
                           "Early",
                           ifelse(samples15$Bulk.Date >= as.Date("2015-07-04") & samples15$Bulk.Date <= as.Date("2015-08-01"),
                                  "Middle",
                                  ifelse(samples15$Bulk.Date >= as.Date("2015-08-02") & samples15$Bulk.Date <= as.Date("2015-08-29"),
                                         "Late",
                                         NA)))

samples15$Strata <- factor(samples15$Strata, levels = c("Early", "Middle", "Late"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samples16 <- read.csv(file = "Harvest/Bulk Tissue Inventory 2016.csv", as.is = TRUE)
str(samples16)

# Convert Date
samples16$Bulk.Date <- as.Date(samples16$Bulk.Date, format = "%Y-%m-%d")
samples16$Bulk.End.Date <- as.Date(samples16$Bulk.End.Date, format = "%Y-%m-%d")
range(samples16$Bulk.Date)

# Create Temporal Strata
samples16$Strata <- ifelse(samples16$Bulk.Date >= as.Date("2016-06-01") & samples16$Bulk.Date <= as.Date("2016-06-27"),
                           "Early",
                           ifelse(samples16$Bulk.Date >= as.Date("2016-06-28") & samples16$Bulk.Date <= as.Date("2016-07-25"),
                                  "Middle",
                                  ifelse(samples16$Bulk.Date >= as.Date("2016-07-26") & samples16$Bulk.Date <= as.Date("2016-08-29"),
                                         "Late",
                                         NA)))

samples16$Strata <- factor(samples16$Strata, levels = c("Early", "Middle", "Late"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

s2014 <- t(data.frame(cast(aggregate(Field.Count ~ Silly.Code + Strata, data = samples14, sum), Silly.Code ~ Strata, value = "Field.Count"))[c(4,5,3,2,1), ])
s2014[is.na(s2014)] <- 0
write.xlsx(x = s2014, file = "V:/Documents/4_Westward/Sockeye/KMA 2014-2016/KMA Mixtures FMS/Table 1.xlsx",
           col.names = FALSE, sheetName = "2014 Samples", append = TRUE)

s2015 <- t(data.frame(cast(aggregate(Field.Count ~ Silly.Code + Strata, data = samples15, sum), Silly.Code ~ Strata, value = "Field.Count"))[c(5,6,4,2,1,3), ])
s2015[is.na(s2015)] <- 0
write.xlsx(x = s2015, file = "V:/Documents/4_Westward/Sockeye/KMA 2014-2016/KMA Mixtures FMS/Table 1.xlsx",
           col.names = FALSE, sheetName = "2015 Samples", append = TRUE)

s2016 <- t(data.frame(cast(aggregate(Field.Count ~ Silly.Code + Strata, data = samples16, sum), Silly.Code ~ Strata, value = "Field.Count"))[c(5,6,4,2,1,3), ])
s2016[is.na(s2016)] <- 0
write.xlsx(x = s2016, file = "V:/Documents/4_Westward/Sockeye/KMA 2014-2016/KMA Mixtures FMS/Table 1.xlsx",
           col.names = FALSE, sheetName = "2016 Samples", append = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Likelihood profiles Ayakulik Early/Late ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

### 103 Loci
## Locus Control
LocusControl <- dget(file = "Objects/LocusControl103.txt")
loci103 <- dget(file = "Objects/loci103.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get Populations
KMA473Pops <- dget(file = "Objects/KMA473Pops.txt")

require(beepr)
invisible(sapply(KMA473Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPopsloci103/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
length(objects(pattern = "\\.gcl"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get objects
KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("OLD", "Likelihood Profiles", "OriginalLocusControl.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 13 Kodiak RGs
dput(x = Kodiak57PopsGroupVec13, file = "Objects/Kodiak57PopsGroupVec13.txt")
Kodiak57PopsGroupVec13 <- dget(file = "Objects/Kodiak57PopsGroupVec13.txt")

dput(x = Groups13.Kodiak, file = "Objects/Groups13.Kodiak.txt")
Groups13.Kodiak <- dget(file = "Objects/Groups13.Kodiak.txt")

Colors13 <- c("darkred", "red", "brown", "cyan", "darkblue", "blue", "purple", "grey", "darkorange", "magenta", "yellow3", "green", "darkgreen")
dput(x = Colors13, file = "Objects/Colors13.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 46 loci

# Kodiak57Pops_13groups_46loci_Likelihood_Profile <- LeaveOneOutDist.GCL(sillyvec = Kodiak57Pops, loci = loci46, groupvec = Kodiak57PopsGroupVec13); beep(5)
# dput(x = Kodiak57Pops_13groups_46loci_Likelihood_Profile, file = "Likelihood Profiles/Kodiak57Pops_13groups_46loci_Likelihood_Profile.txt")
Kodiak57Pops_13groups_46loci_Likelihood_Profile <- dget(file = "Likelihood Profiles/Kodiak57Pops_13groups_46loci_Likelihood_Profile.txt")

# Kodiak57Pops_13groups_46loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = Kodiak57Pops_13groups_46loci_Likelihood_Profile, groupnames = Groups13.Kodiak, groupvec = Kodiak57PopsGroupVec13, sillyvec = Kodiak57Pops); beep(2)
# dput(x = Kodiak57Pops_13groups_46loci_Confusion, file = "Objects/Kodiak57Pops_13groups_46loci_Confusion.txt")
Kodiak57Pops_13groups_46loci_Confusion <- dget(file = "Objects/Kodiak57Pops_13groups_46loci_Confusion.txt")

# Group to Group
require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(Kodiak57Pops_13groups_46loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", 
          ylab = "Mean Genotype Likelihood", main = "loci46", at = seq(0, 1, length.out = 100), 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })
str(Kodiak57Pops_13groups_46loci_Confusion)
abline(v = 1, col = "red")
apply(Kodiak57Pops_13groups_46loci_Confusion[[1]], 2, sum)



# Pop to Pop
low <- c(sapply(seq(Groups13.Kodiak), function(grp) {which.max(sort(Kodiak57PopsGroupVec13) == grp)} ), length(Kodiak57PopsGroupVec13) + 1) - 0.5

levelplot(Kodiak57Pops_13groups_46loci_Confusion[[3]], col.regions = new.colors, xlab = "Known Origin", 
          ylab = "Mean Genotype Likelihood", main = "loci46", at = seq(0, 1, length.out = 100), 
          scales = list(x = list(rot = 0), labels = 1:57),
          panel = function(...) {
            panel.levelplot(...)
            panel.segments(x0 = low[1:13], x1 = low[1:13], y0 = low[1:13], y1 = low[2:14], col = Colors13, lwd = 3)
            panel.segments(x0 = low[2:14], x1 = low[2:14], y0 = low[1:13], y1 = low[2:14], col = Colors13, lwd = 3)
            panel.segments(x0 = low[1:13], x1 = low[2:14], y0 = low[1:13], y1 = low[1:13], col = Colors13, lwd = 3)
            panel.segments(x0 = low[1:13], x1 = low[2:14], y0 = low[2:14], y1 = low[2:14], col = Colors13, lwd = 3)
          })



rownames(Kodiak57Pops_13groups_46loci_Confusion[[3]])

# Assignment into Boxplots
suppressWarnings(invisible(lapply(Kodiak57Pops_13groups_46loci_Likelihood_Profile[[1]], function(lst) {boxplot(lst, pch=16, notch=TRUE, col=Colors13[sort(Kodiak57PopsGroupVec13)], ylab="Probability")} )))  # Regional Flat Prior

# Jim's New
# Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW <- LeaveOneOutLikeProfile.GCL(popvec = Kodiak57Pops, loci = loci46, groupvec = Kodiak57PopsGroupVec13, groupnames = Groups13.Kodiak, ncores = 4)
# dput(x = Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW, file = "Likelihood Profiles/Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW.txt")
Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW <- dget(file = "Likelihood Profiles/Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW.txt")
str(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW)
head(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$IndividualByGroup)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Heatmap of individuals and which RG they go to
which(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromGroup == "Frazer")

RGpops <- Kodiak57Pops[which(Kodiak57PopsGroupVec13 == which(Groups13.Kodiak == "Frazer"))]

low <- c(sapply(RGpops, function(pop) {which.max(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromPop == pop)} ) - which.max(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromPop == RGpops[1]),
         length(which(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromGroup == "Frazer")) - 1) + 0.5

levelplot(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$IndividualByGroup[which(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromGroup == "Frazer"), ],
          col.regions = new.colors, xlab = "Individuals", ylab = "Genotype Likelihood", main = "loci46", aspect = "fill",
          panel = function(...) {
            panel.levelplot(...)
            panel.segments(x0 = low[1:(length(low)-1)], x1 = low[1:(length(low)-1)], y0 = low[1:(length(low)-1)], y1 = low[2:length(low)], col = "red", lwd = 3)
            panel.segments(x0 = low[2:length(low)], x1 = low[2:length(low)], y0 = low[1:(length(low)-1)], y1 = low[2:length(low)], col = "red", lwd = 3)
            panel.segments(x0 = low[1:(length(low)-1)], x1 = low[2:length(low)], y0 = low[1:(length(low)-1)], y1 = low[1:(length(low)-1)], col = "red", lwd = 3)
            panel.segments(x0 = low[1:(length(low)-1)], x1 = low[2:length(low)], y0 = low[2:length(low)], y1 = low[2:length(low)], col = "red", lwd = 3)
          })

# Heatmap of individuals and which Pop they go to
levelplot(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$IndividualByPop[which(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromGroup == "Frazer"), ],
          col.regions = new.colors, xlab = "Individuals", ylab = "Genotype Likelihood", main = "loci46", aspect = "fill")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






# Heatmap of individuals and which Pop they go to
sapply(Groups13.Kodiak, function(Group) {
  myplot <- levelplot(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$IndividualByPop[which(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromGroup == Group), ],
                      col.regions = new.colors, xlab = "Individuals", ylab = "Genotype Likelihood", main = paste("loci46: Individuals from ", Group, sep = ''), aspect = "fill")
  print(myplot)
})







sapply(Groups13.Kodiak, function (Group) {boxplot(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$IndividualByGroup[which(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromGroup == Group), ], pch = 16, notch = TRUE, col = Colors13, main = Group)} )
sapply(Groups13.Kodiak, function (Group) {boxplot(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$IndividualByPop[which(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromGroup == Group), ], pch = 16, notch = TRUE, col = Colors13[Kodiak57PopsGroupVec13], main = Group)} )



PlotLikeProfile.GCL(likeprof = Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW, popvec = Kodiak57Pops, loci = loci46, groupvec = Kodiak57PopsGroupVec13, groupnames = Groups13.Kodiak, dir = "Likelihood Profiles", filename = "Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW", col = Colors13)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 89 loci
# Kodiak57Pops_13groups_89loci_Likelihood_Profile <- LeaveOneOutDist.GCL(sillyvec = Kodiak57Pops, loci = loci89, groupvec = Kodiak57PopsGroupVec13); beep(5)
# dput(x = Kodiak57Pops_13groups_89loci_Likelihood_Profile, file = "Likelihood Profiles/Kodiak57Pops_13groups_89loci_Likelihood_Profile.txt")
Kodiak57Pops_13groups_89loci_Likelihood_Profile <- dget(file = "Likelihood Profiles/Kodiak57Pops_13groups_89loci_Likelihood_Profile.txt")

# Kodiak57Pops_13groups_89loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = Kodiak57Pops_13groups_89loci_Likelihood_Profile, groupnames = Groups13.Kodiak, groupvec = Kodiak57PopsGroupVec13, sillyvec = Kodiak57Pops); beep(2)
# dput(x = Kodiak57Pops_13groups_89loci_Confusion, file = "Objects/Kodiak57Pops_13groups_89loci_Confusion.txt")
Kodiak57Pops_13groups_89loci_Confusion <- dget(file = "Objects/Kodiak57Pops_13groups_89loci_Confusion.txt")

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(Kodiak57Pops_13groups_89loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci89", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))
str(Kodiak57Pops_13groups_89loci_Confusion)
apply(Kodiak57Pops_13groups_89loci_Confusion[[1]], 2, sum)

low <- c(sapply(seq(Groups13.Kodiak), function(grp) {which.max(sort(Kodiak57PopsGroupVec13) == grp)} ), length(Kodiak57PopsGroupVec13) + 1) - 0.5

levelplot(Kodiak57Pops_13groups_89loci_Confusion[[3]], col.regions = new.colors, xlab = "Known Origin", 
          ylab = "Mean Genotype Likelihood", main = "loci89", at = seq(0, 1, length.out = 100), 
          scales = list(x = list(rot = 0), labels = 1:57),
          panel = function(...) {
            panel.levelplot(...)
            panel.segments(x0 = low[1:13], x1 = low[1:13], y0 = low[1:13], y1 = low[2:14], col = Colors13, lwd = 3)
            panel.segments(x0 = low[2:14], x1 = low[2:14], y0 = low[1:13], y1 = low[2:14], col = Colors13, lwd = 3)
            panel.segments(x0 = low[1:13], x1 = low[2:14], y0 = low[1:13], y1 = low[1:13], col = Colors13, lwd = 3)
            panel.segments(x0 = low[1:13], x1 = low[2:14], y0 = low[2:14], y1 = low[2:14], col = Colors13, lwd = 3)
          })

rownames(Kodiak57Pops_13groups_89loci_Confusion[[3]])

# Assignment into Boxplots
suppressWarnings(invisible(lapply(Kodiak57Pops_13groups_89loci_Likelihood_Profile[[1]], function(lst) {boxplot(lst, pch=16, notch=TRUE, col=Colors13[sort(Kodiak57PopsGroupVec13)], ylab="Probability")} )))  # Regional Flat Prior


# Kodiak57Pops_13groups_89loci_Likelihood_ProfileNEW <- LeaveOneOutLikeProfile.GCL(popvec = Kodiak57Pops, loci = loci89, groupvec = Kodiak57PopsGroupVec13, groupnames = Groups13.Kodiak, ncores = 4)
# dput(x = Kodiak57Pops_13groups_89loci_Likelihood_ProfileNEW, file = "Likelihood Profiles/Kodiak57Pops_13groups_89loci_Likelihood_ProfileNEW.txt")
Kodiak57Pops_13groups_89loci_Likelihood_ProfileNEW <- dget(file = "Likelihood Profiles/Kodiak57Pops_13groups_89loci_Likelihood_ProfileNEW.txt")

sapply(Groups13.Kodiak, function (Group) {boxplot(Kodiak57Pops_13groups_89loci_Likelihood_ProfileNEW$IndividualByGroup[which(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromGroup == Group), ], pch = 16, notch = TRUE, col = Colors13, main = Group)} )
sapply(Groups13.Kodiak, function (Group) {boxplot(Kodiak57Pops_13groups_89loci_Likelihood_ProfileNEW$IndividualByPop[which(Kodiak57Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromGroup == Group), ], pch = 16, notch = TRUE, col = Colors13[Kodiak57PopsGroupVec13], main = Group)} )



PlotLikeProfileKS.GCL(likeprof = Kodiak57Pops_13groups_89loci_Likelihood_ProfileNEW, popvec = Kodiak57Pops, loci = loci89, groupvec = Kodiak57PopsGroupVec13, groupnames = Groups13.Kodiak, dir = "Likelihood Profiles", filename = "Kodiak57Pops_13groups_89loci_Likelihood_ProfileNEW", col = Colors13)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 46 loci No Ayakulik Weir Samples!!!!
weirpops <- c(grep(pattern = "SAYAK00.SAYAK08L.SAYAK11L.SAYAK12L", x = Kodiak57Pops), grep(pattern = "SAYAK12E", x = Kodiak57Pops))

Kodiak55Pops <- Kodiak57Pops[-weirpops]
Kodiak55PopsGroupVec13 <- Kodiak57PopsGroupVec13[-weirpops]


# Kodiak55Pops_13groups_46loci_Likelihood_Profile <- LeaveOneOutDist.GCL(sillyvec = Kodiak55Pops, loci = loci46, groupvec = Kodiak55PopsGroupVec13); beep(5)
# dput(x = Kodiak55Pops_13groups_46loci_Likelihood_Profile, file = "Likelihood Profiles/Kodiak55Pops_13groups_46loci_Likelihood_Profile.txt")
Kodiak55Pops_13groups_46loci_Likelihood_Profile <- dget(file = "Likelihood Profiles/Kodiak55Pops_13groups_46loci_Likelihood_Profile.txt")

# Kodiak55Pops_13groups_46loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = Kodiak55Pops_13groups_46loci_Likelihood_Profile, groupnames = Groups13.Kodiak, groupvec = Kodiak55PopsGroupVec13, sillyvec = Kodiak55Pops); beep(2)
# dput(x = Kodiak55Pops_13groups_46loci_Confusion, file = "Objects/Kodiak55Pops_13groups_46loci_Confusion.txt")
Kodiak55Pops_13groups_46loci_Confusion <- dget(file = "Objects/Kodiak55Pops_13groups_46loci_Confusion.txt")

# Pop to Pop
low <- c(sapply(seq(Groups13.Kodiak), function(grp) {which.max(sort(Kodiak55PopsGroupVec13) == grp)} ), length(Kodiak55PopsGroupVec13) + 1) - 0.5

levelplot(Kodiak55Pops_13groups_46loci_Confusion[[3]], col.regions = new.colors, xlab = "Known Origin", 
          ylab = "Mean Genotype Likelihood", main = "loci46", at = seq(0, 1, length.out = 100), 
          scales = list(x = list(rot = 0), labels = 1:55),
          panel = function(...) {
            panel.levelplot(...)
            panel.segments(x0 = low[1:13], x1 = low[1:13], y0 = low[1:13], y1 = low[2:14], col = Colors13, lwd = 3)
            panel.segments(x0 = low[2:14], x1 = low[2:14], y0 = low[1:13], y1 = low[2:14], col = Colors13, lwd = 3)
            panel.segments(x0 = low[1:13], x1 = low[2:14], y0 = low[1:13], y1 = low[1:13], col = Colors13, lwd = 3)
            panel.segments(x0 = low[1:13], x1 = low[2:14], y0 = low[2:14], y1 = low[2:14], col = Colors13, lwd = 3)
          })


# Jim's
# Kodiak55Pops_13groups_46loci_Likelihood_ProfileNEW <- LeaveOneOutLikeProfile.GCL(popvec = Kodiak55Pops, loci = loci46, groupvec = Kodiak55PopsGroupVec13, groupnames = Groups13.Kodiak, ncores = 4)
# dput(x = Kodiak55Pops_13groups_46loci_Likelihood_ProfileNEW, file = "Likelihood Profiles/Kodiak55Pops_13groups_46loci_Likelihood_ProfileNEW.txt")
Kodiak55Pops_13groups_46loci_Likelihood_ProfileNEW <- dget(file = "Likelihood Profiles/Kodiak55Pops_13groups_46loci_Likelihood_ProfileNEW.txt")
PlotLikeProfile.GCL(likeprof = Kodiak55Pops_13groups_46loci_Likelihood_ProfileNEW, popvec = Kodiak55Pops, loci = loci46, groupvec = Kodiak55PopsGroupVec13, groupnames = Groups13.Kodiak, dir = "Likelihood Profiles", filename = "Kodiak55Pops_13groups_46loci_Likelihood_ProfileNEW", col = Colors13)

# Heatmap of individuals and which Pop they go to
sapply(Groups13.Kodiak, function(Group) {
  myplot <- levelplot(Kodiak55Pops_13groups_46loci_Likelihood_ProfileNEW$IndividualByPop[which(Kodiak55Pops_13groups_46loci_Likelihood_ProfileNEW$Attributes$FromGroup == Group), ],
                      col.regions = new.colors, xlab = "Individuals", ylab = "Genotype Likelihood", main = paste("loci46: Individuals from ", Group, sep = ''), aspect = "fill")
  print(myplot)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Looking into HWE and Fis for Ayakulik Weir samples
KMA473Pops[186:187]


HWE <- ReadGenepopHWE.GCL(file = "Genepop/KMA473Pops_93nuclearloci_MCMC.txt.P")

# Histogram of Fis for SAYAK12E
hist(HWE$DataByPop[as.numeric(HWE$DataByPop$Pop) == 187, "WC Fis"], col=8, xlim = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.01))

# How does distribution of Fis compare to other Frazer/Ayakulik pops
invisible(sapply(which(KMA473PopsGroupVec17 %in% c(5,6,7)), function(pop) {
  hist(HWE$DataByPop[as.numeric(HWE$DataByPop$Pop) == pop, "WC Fis"], 
       col=8, xlim = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.01), main = KMA473Pops[pop])
}) )


# How do HWE p-values compare to other Frazer/Ayakulik pops
invisible(sapply(which(KMA473PopsGroupVec17 %in% c(5,6,7)), function(pop) {
  hist(HWE$DataByPop[as.numeric(HWE$DataByPop$Pop) == pop, "PValue"], 
       col=8, xlim = c(0, 1), breaks = seq(0, 1, 0.05), main = KMA473Pops[pop])
}) )

# Do not appear to see any p-value anomolies



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Quick look at Karluk / Frazer / Ayakulik PostQC Collections MDS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 96 Loci
## Locus Control
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
loci96 <- dget(file = "Objects/loci96.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### SSUMMM09
LOKI2R.GCL(sillyvec = "SSUMMM09", username = "krshedd", password = password)
rm(password)

str(SSUMMM09.gcl)

dput(x = SSUMMM09.gcl, file = "Raw genotypes/OriginalCollections/SSUMMM09.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC SSUMMM09
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Check loci
## Get sample size by locus
Original_SSUMMM09_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = "SSUMMM09", loci = loci96)
min(Original_SSUMMM09_SampleSizebyLocus)  ## 85/95
apply(Original_SSUMMM09_SampleSizebyLocus, 1, min) / apply(Original_SSUMMM09_SampleSizebyLocus, 1, max)  ## Good, 0.90

Original_SSUMMM09_PercentbyLocus <- apply(Original_SSUMMM09_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_SSUMMM09_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_SSUMMM09_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


#### Check individuals
## View Histogram of Failure Rate by Strata
invisible(sapply("SSUMMM09", function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  failure <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(ind == "0") / length(ind)} )
  hist(x = failure, main = mix, xlab = "Failure Rate", col = 8, xlim = c(0, 1), ylim = c(0, 20), breaks = seq(from = 0, to = 1, by = 0.02))
  abline(v = 0.2, lwd = 3)
}))



SSUMMM09_SampleSizes <- matrix(data = NA, nrow = length("SSUMMM09"), ncol = 5, 
                                         dimnames = list("SSUMMM09", c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_SSUMMM09_ColSize <- sapply(paste("SSUMMM09", ".gcl", sep = ''), function(x) get(x)$n)
SSUMMM09_SampleSizes[, "Genotyped"] <- Original_SSUMMM09_ColSize


### Alternate
## Indentify alternate species individuals
SSUMMM09_Alternate <- FindAlternateSpecies.GCL(sillyvec = "SSUMMM09", species = "sockeye")

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport = SSUMMM09_Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)

## Get number of individuals per silly after removing alternate species individuals
ColSize_SSUMMM09_PostAlternate <- sapply(paste("SSUMMM09", ".gcl", sep = ''), function(x) get(x)$n)
SSUMMM09_SampleSizes[, "Alternate"] <- Original_SSUMMM09_ColSize-ColSize_SSUMMM09_PostAlternate


### Missing
## Remove individuals with >20% missing data
SSUMMM09_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = "SSUMMM09", proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SSUMMM09_PostMissLoci <- sapply(paste("SSUMMM09", ".gcl", sep = ''), function(x) get(x)$n)
SSUMMM09_SampleSizes[, "Missing"] <- ColSize_SSUMMM09_PostAlternate-ColSize_SSUMMM09_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
SSUMMM09_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = "SSUMMM09", loci = loci96, quantile = NULL, minproportion = 0.95)
SSUMMM09_DuplicateCheckReportSummary <- sapply("SSUMMM09", function(x) SSUMMM09_DuplicateCheck95MinProportion[[x]]$report)
SSUMMM09_DuplicateCheckReportSummary

## Remove duplicate individuals
SSUMMM09_RemovedDups <- RemoveDups.GCL(SSUMMM09_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SSUMMM09_PostDuplicate <- sapply(paste("SSUMMM09", ".gcl", sep = ''), function(x) get(x)$n)
SSUMMM09_SampleSizes[, "Duplicate"] <- ColSize_SSUMMM09_PostMissLoci-ColSize_SSUMMM09_PostDuplicate


### Final
SSUMMM09_SampleSizes[, "Final"] <- ColSize_SSUMMM09_PostDuplicate
SSUMMM09_SampleSizes

# Lost two fish for missing.

dput(x = SSUMMM09.gcl, file = "Raw genotypes/PostQCCollections/SSUMMM09.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get PostQC Collections
# FrazerKarlukAyakulik41Collections <- readClipboard()
# dput(x = FrazerKarlukAyakulik41Collections, file = "Objects/FrazerKarlukAyakulik41Collections.txt")
FrazerKarlukAyakulik41Collections <- dget(file = "Objects/FrazerKarlukAyakulik41Collections.txt")

FrazerKarlukAyakulik42Collections <- c("SSUMMM09", FrazerKarlukAyakulik41Collections)
# dput(x = FrazerKarlukAyakulik42Collections, file = "Objects/FrazerKarlukAyakulik42Collections.txt")


require(beepr)
invisible(sapply(FrazerKarlukAyakulik42Collections, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCCollections/", silly, ".txt", sep = "")), pos = 1)})); beep(2)
length(objects(pattern = "\\.gcl"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get objects
KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("OLD", "Likelihood Profiles", "OriginalLocusControl.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FrazerKarlukAyakulik41CollectionsGroupVec3 <- as.numeric(readClipboard())
# dput(x = FrazerKarlukAyakulik41CollectionsGroupVec3, file = "Objects/FrazerKarlukAyakulik41CollectionsGroupVec3.txt")
FrazerKarlukAyakulik41CollectionsGroupVec3 <- dget(file = "Objects/FrazerKarlukAyakulik41CollectionsGroupVec3.txt")

FrazerKarlukAyakulik42CollectionsGroupVec3 <- c(5, FrazerKarlukAyakulik41CollectionsGroupVec3)
# dput(x = FrazerKarlukAyakulik41CollectionsGroupVec3, file = "Objects/FrazerKarlukAyakulik41CollectionsGroupVec3.txt")


# FrazerKarlukAyakulik41CollectionsCommonNames <- readClipboard()
# dput(x = FrazerKarlukAyakulik41CollectionsCommonNames, file = "Objects/FrazerKarlukAyakulik41CollectionsCommonNames.txt")
FrazerKarlukAyakulik41CollectionsCommonNames <- dget(file = "Objects/FrazerKarlukAyakulik41CollectionsCommonNames.txt")

FrazerKarlukAyakulik42CollectionsCommonNames <- c("Frazer Lake - Summit Creek", FrazerKarlukAyakulik41CollectionsCommonNames)
# dput(x = FrazerKarlukAyakulik42CollectionsCommonNames, file = "Objects/FrazerKarlukAyakulik42CollectionsCommonNames.txt")


Groups3 <- c("Frazer", "Ayakulik", "Karluk")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Allele Frequencies
# Counts
FrazerKarlukAyakulik42CollectionsPostRemovalsAlleleCounts <- FreqPop.GCL(sillyvec = FrazerKarlukAyakulik42Collections, loci = loci96)
str(FrazerKarlukAyakulik42CollectionsPostRemovalsAlleleCounts)

# Frequencies
FrazerKarlukAyakulik42CollectionsPostRemovalsFreqs <- FrazerKarlukAyakulik42CollectionsPostRemovalsAlleleCounts[,,"Allele 1"] / (FrazerKarlukAyakulik42CollectionsPostRemovalsAlleleCounts[,,"Allele 2"] + FrazerKarlukAyakulik42CollectionsPostRemovalsAlleleCounts[,,"Allele 1"])
str(FrazerKarlukAyakulik42CollectionsPostRemovalsFreqs)

## Plot points and smoothing line
getwd()
pdf(file = "FreqPlots/FrazerKarlukAyakulik42Collections_96SNPs_GroupVec3_FreqPlots.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
par(mar = c(2.1, 4.1, 4.1, 4.6))
for(locus in loci96){
  plot(FrazerKarlukAyakulik42CollectionsPostRemovalsFreqs[, locus], main = locus, col = Colors15[FrazerKarlukAyakulik42CollectionsGroupVec3], pch = 19, ylim = c(0, 1), ylab = "Frequency", xlab = "Collections", xaxt = 'n', cex = 0.7)
  lines(supsmu(seq(length(FrazerKarlukAyakulik42Collections)), FrazerKarlukAyakulik42CollectionsPostRemovalsFreqs[, locus]), lwd = 2)
  mtext(text = "Collections", side = 1, line = 0.7)
  text(x = seq(FrazerKarlukAyakulik42Collections), y = FrazerKarlukAyakulik42CollectionsPostRemovalsFreqs[, locus], labels = seq(FrazerKarlukAyakulik42Collections), cex = 0.3, col = c("black", "white", "white")[FrazerKarlukAyakulik42CollectionsGroupVec3 - 4])
  legend(x = 44, y = 0.9, legend = Groups3, fill = Colors15[5:7], xpd = TRUE, bty = "n", cex = 0.7)
}; rm(locus)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Combine Loci
CombineLoci.GCL(sillyvec = FrazerKarlukAyakulik42Collections, markerset = c("One_CO1", "One_Cytb_17", "One_Cytb_26"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = FrazerKarlukAyakulik42Collections, markerset = c("One_Tf_ex10-750", "One_Tf_ex3-182"), update = TRUE); beep(2)

loci89

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MDS 89 loci
# gcl2Genepop.GCL(sillyvec = FrazerKarlukAyakulik42Collections, loci = loci89, path = "Genepop/FrazerKarlukAyakulik42Collections_89loci.gen", VialNums = TRUE)

detach("package:adegenet", unload = TRUE)
require(package = adegenet, lib.loc = "C:/Users/krshedd/Documents/R/win-library/3.1")

genind <- read.genepop(file = "Genepop/FrazerKarlukAyakulik42Collections_89loci.gen")

genpop <- genind2genpop(genind)

# AdegenetNei42Col89loci <- dist.genpop(genpop, method = 1, diag = TRUE, upper = TRUE)
# dput(x = AdegenetNei42Col89loci, file = "Trees/AdegenetNei42Col89loci.txt")
AdegenetNei42Col89loci <- dget(file = "Trees/AdegenetNei42Col89loci.txt")
str(AdegenetNei42Col89loci)

require(ape)
Nei42NJtree <- nj(AdegenetNei42Col89loci)
str(Nei42NJtree)
Nei42NJtree$tip.label <- readClipboard()
plot.phylo(x = Nei42NJtree, cex = 0.5, no.margin = TRUE)

library('rgl')

MDS <- cmdscale(as.matrix(AdegenetNei42Col89loci), k = 3)
# MDS <- cmdscale(as.matrix(AdegenetNei42Col89loci), k = 40, eig = TRUE)  # Do with all possible dimensions
# dput(x = MDS, file = "Objects/MDSAdegenetNei42Col89loci.txt")
# dput(x = MDS, file = "Objects/MDSAdegenetNei42Col89loci_alldim.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei42Col89loci.txt")

x <- as.vector(MDS[, 1])   
y <- as.vector(MDS[, 2])
z <- as.vector(MDS[, 3])

FrazerKarlukAyakulik42CollectionsGroupVec3
Colors15 <- dget(file = "Objects/Colors15.txt")


open3d()
par(family = "serif")
plot3d(x, y, z + abs(range(z)[1]), xlab = '', ylab = '', zlab = '', aspect = FALSE, col = Colors15[FrazerKarlukAyakulik42CollectionsGroupVec3], size = 0.5, type = 's', axes = TRUE, box = TRUE, top = TRUE, cex = 1)
box3d()
# plot3d(x, y, z + abs(range(z)[1]), aspect = FALSE, col = "black", size = 0.5, type = 'h', box = TRUE, axes = FALSE, top = FALSE, add = TRUE)  # adds pins to spheres
# texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.8, 0.8), text = seq(FrazerKarlukAyakulik42Collections), font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)
texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.2, 0.2), text = FrazerKarlukAyakulik42Collections, font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)

rgl.snapshot("MDS/MDSAdegenetNei42ColFrazerKarlukAyakulik89loci2.png", fmt="png", top=TRUE )
# SSUMMM09 looks very similar SMIDWM08


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create gsi_sim input file
dir.create("gsi_sim")
LocusControl <- dget(file = "Objects/LocusControl103.txt")
write_snps_gsi_sim.GCL(collections = FrazerKarlukAyakulik42Collections, loci = c(loci89[-which(LocusControl$ploidy[loci89] == 1)], "One_Tf_ex3-182"), path.name = "gsi_sim/FrazerKarlukAyakulik42Collections88diploidloci.txt")



setwd("gsi_sim")  # Where you want your output to be saved

setwd("C:/Users/krshedd/Documents/R/gsi_sim")  # Where you want your output to be saved
dir.create("AyakulikFrazerKarluk")
file.copy(from = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/gsi_sim/FrazerKarlukAyakulik42Collections88diploidloci.txt",
          to = "AyakulikFrazerKarluk")  # copy baseline from V: drive
setwd("AyakulikFrazerKarluk")

system2(command = "C:/Users/krshedd/Documents/R/gsi_sim/gsi_sim-MINGW64_NT-6.1",  # Where the executable is
        args = c("-b FrazerKarlukAyakulik42Collections88diploidloci.txt",  # Where the baseline file is
                 "--self-assign"), 
        stdout = "dumpfile")  # What to name the big output file

# Splitting output with Eric's code, may not work, depends on baseline format and naming conventions!
system2(command = shQuote("awk -F";" 'BEGIN {print "ID TopPop Score"} /SELF_ASSIGN_A_LA_GC_CSV:/ {print $1, $2, $3}' FrazerKarlukAyakulik42Collections88diploidloci_Output.txt | sed 's/SELF_ASSIGN_A_LA_GC_CSV:\///g;'"),
        stdout = "self-ass-results.txt")


x <- paste0("awk -F\";\" \'BEGIN {print \"ID TopPop Score\"} /SELF_ASSIGN_A_LA_GC_CSV:/ {print $1, $2, $3}\' dumpfile | sed \'s/SELF_ASSIGN_A_LA_GC_CSV:\\///g;\'")
writeLines(x)
system2(command = writeLines(x),
        stdout = "self-ass-results.txt")


# Typed into command line bash and that worked...

library(dplyr)
library(stringr)

AssignmentTally <- read.table("self-ass-results.txt", header = TRUE) %>%
  tbl_df %>%
  mutate(FromPop = paste(str_split_fixed(ID, "_", 3)[,1], str_split_fixed(ID, "_", 3)[,2], sep ="_")) %>%
  group_by(FromPop, TopPop) %>%
  tally %>%
  arrange(FromPop, desc(n))
str(AssignmentTally)
as.data.frame(AssignmentTally)[1:25,]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dog Salmon Weir Collections 2011-2012 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 96 Loci
## Locus Control
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
loci96 <- dget(file = "Objects/loci96.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Dog Salmon Weir Collections
DogSalmonCollections5 <- c("SDOGSC11E", "SDOGSC11EJ", "SDOGSC11W", "SDOGSC12E", "SDOGSC12W")

LOKI2R.GCL(sillyvec = DogSalmonCollections5, username = "krshedd", password = password)
rm(password)


## Save unaltered .gcl's as back-up:
invisible(sapply(DogSalmonCollections5, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/" , silly, ".txt", sep = ''))} )); beep(5)

## Original sample sizes by SILLY
collection.size.original <- sapply(DogSalmonCollections5, function(silly) get(paste(silly, ".gcl", sep = ""))$n)
collection.size.original

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC DogSalmonCollections5
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Check loci
## Get sample size by locus
Original_DogSalmonCollections5_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = DogSalmonCollections5, loci = loci96)
min(Original_DogSalmonCollections5_SampleSizebyLocus)  ## 68
apply(Original_DogSalmonCollections5_SampleSizebyLocus, 1, min) / apply(Original_DogSalmonCollections5_SampleSizebyLocus, 1, max)  ## Good, 0.91

Original_DogSalmonCollections5_PercentbyLocus <- apply(Original_DogSalmonCollections5_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_DogSalmonCollections5_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_DogSalmonCollections5_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


#### Check individuals
## View Histogram of Failure Rate by Strata
invisible(sapply(DogSalmonCollections5, function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  failure <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(ind == "0") / length(ind)} )
  hist(x = failure, main = mix, xlab = "Failure Rate", col = 8, xlim = c(0, 1), ylim = c(0, 20), breaks = seq(from = 0, to = 1, by = 0.02))
  abline(v = 0.2, lwd = 3)
}))



DogSalmonCollections5_SampleSizes <- matrix(data = NA, nrow = length(DogSalmonCollections5), ncol = 5, 
                               dimnames = list(DogSalmonCollections5, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_DogSalmonCollections5_ColSize <- sapply(paste(DogSalmonCollections5, ".gcl", sep = ''), function(x) get(x)$n)
DogSalmonCollections5_SampleSizes[, "Genotyped"] <- Original_DogSalmonCollections5_ColSize


### Alternate
## Indentify alternate species individuals
DogSalmonCollections5_Alternate <- FindAlternateSpecies.GCL(sillyvec = DogSalmonCollections5, species = "sockeye")

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport = DogSalmonCollections5_Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)

## Get number of individuals per silly after removing alternate species individuals
ColSize_DogSalmonCollections5_PostAlternate <- sapply(paste(DogSalmonCollections5, ".gcl", sep = ''), function(x) get(x)$n)
DogSalmonCollections5_SampleSizes[, "Alternate"] <- Original_DogSalmonCollections5_ColSize-ColSize_DogSalmonCollections5_PostAlternate


### Missing
## Remove individuals with >20% missing data
DogSalmonCollections5_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = DogSalmonCollections5, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_DogSalmonCollections5_PostMissLoci <- sapply(paste(DogSalmonCollections5, ".gcl", sep = ''), function(x) get(x)$n)
DogSalmonCollections5_SampleSizes[, "Missing"] <- ColSize_DogSalmonCollections5_PostAlternate-ColSize_DogSalmonCollections5_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
DogSalmonCollections5_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = DogSalmonCollections5, loci = loci96, quantile = NULL, minproportion = 0.95)
DogSalmonCollections5_DuplicateCheckReportSummary <- sapply(DogSalmonCollections5, function(x) DogSalmonCollections5_DuplicateCheck95MinProportion[[x]]$report)
DogSalmonCollections5_DuplicateCheckReportSummary

## Remove duplicate individuals
DogSalmonCollections5_RemovedDups <- RemoveDups.GCL(DogSalmonCollections5_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_DogSalmonCollections5_PostDuplicate <- sapply(paste(DogSalmonCollections5, ".gcl", sep = ''), function(x) get(x)$n)
DogSalmonCollections5_SampleSizes[, "Duplicate"] <- ColSize_DogSalmonCollections5_PostMissLoci-ColSize_DogSalmonCollections5_PostDuplicate


### Final
DogSalmonCollections5_SampleSizes[, "Final"] <- ColSize_DogSalmonCollections5_PostDuplicate
DogSalmonCollections5_SampleSizes

# Lost two fish for missing.

invisible(sapply(DogSalmonCollections5, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostQCCollections/" , silly, ".txt", sep = ''))} )); beep(5)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check HWE in collections (only use diploid loci)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mito.loci <- dget(file = "Objects/mito.loci.txt")

gcl2Genepop.GCL(sillyvec = DogSalmonCollections5, loci = loci96[-mito.loci], path = "Genepop/DogSalmonCollections5_93nuclearloci.txt", VialNums = TRUE); beep(2)

## Check HWE in Genepop
# Originally done with Exact test and subsequently re-done with MCMC
HWE <- ReadGenepopHWE.GCL(file = "Genepop/DogSalmonCollections5_93nuclearloci.txt.P")
str(HWE)


sapply(as.character(unique(HWE$DataByPop$Pop)), function(pop) {
  x <- subset(x = HWE$DataByPop, subset = Pop == pop)
  plot(sort(x[, "WC Fis"]), type = "h", lwd = 5, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40", main = pop)
  abline(h = 0, lwd = 5)
  x[x$PValue[!is.na(x$PValue)] < 0.05, ]
}, USE.NAMES = TRUE, simplify = FALSE
)

HWE$SummaryPValues


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add these weir collections to MDS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FrazerKarlukAyakulik41Collections <- dget(file = "Objects/FrazerKarlukAyakulik41Collections.txt")
FrazerKarlukAyakulik42Collections <- c("SSUMMM09", FrazerKarlukAyakulik41Collections)

invisible(sapply(FrazerKarlukAyakulik42Collections, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCCollections/", silly, ".txt", sep = "")), pos = 1)})); beep(2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get objects
KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("OLD", "Likelihood Profiles", "OriginalLocusControl.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Combine Loci
CombineLoci.GCL(sillyvec = c(DogSalmonCollections5, FrazerKarlukAyakulik42Collections), markerset = c("One_CO1", "One_Cytb_17", "One_Cytb_26"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = c(DogSalmonCollections5, FrazerKarlukAyakulik42Collections), markerset = c("One_Tf_ex10-750", "One_Tf_ex3-182"), update = TRUE); beep(2)

loci89

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MDS 89 loci
 gcl2Genepop.GCL(sillyvec = c(DogSalmonCollections5, FrazerKarlukAyakulik42Collections), loci = loci89, path = "Genepop/FrazerKarlukAyakulik47Collections_89loci.gen", VialNums = TRUE)

detach("package:adegenet", unload = TRUE)
require(package = adegenet, lib.loc = "C:/Users/krshedd/Documents/R/win-library/3.1")

genind <- read.genepop(file = "Genepop/FrazerKarlukAyakulik47Collections_89loci.gen")

genpop <- genind2genpop(genind)

# AdegenetNei47Col89loci <- dist.genpop(genpop, method = 1, diag = TRUE, upper = TRUE)
# dput(x = AdegenetNei47Col89loci, file = "Trees/AdegenetNei47Col89loci.txt")
AdegenetNei47Col89loci <- dget(file = "Trees/AdegenetNei47Col89loci.txt")
str(AdegenetNei47Col89loci)

require(ape)
Nei47NJtree <- nj(AdegenetNei47Col89loci)
str(Nei47NJtree)
Nei47NJtree$tip.label <- readClipboard()
plot.phylo(x = Nei47NJtree, cex = 0.5, no.margin = TRUE)

library('rgl')

MDS <- cmdscale(as.matrix(AdegenetNei47Col89loci), k = 3)
# MDS <- cmdscale(as.matrix(AdegenetNei47Col89loci), k = 40, eig = TRUE)  # Do with all possible dimensions
# dput(x = MDS, file = "Objects/MDSAdegenetNei47Col89loci.txt")
# dput(x = MDS, file = "Objects/MDSAdegenetNei47Col89loci_alldim.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei47Col89loci.txt")

x <- as.vector(MDS[, 1])   
y <- as.vector(MDS[, 2])
z <- as.vector(MDS[, 3])


FrazerKarlukAyakulik47Collections <- c(DogSalmonCollections5, FrazerKarlukAyakulik42Collections)
FrazerKarlukAyakulik42CollectionsGroupVec3 <- c(5, FrazerKarlukAyakulik41CollectionsGroupVec3)
FrazerKarlukAyakulik47CollectionsGroupVec3 <- c(rep(5, 5), FrazerKarlukAyakulik42CollectionsGroupVec3)
Colors15 <- dget(file = "Objects/Colors15.txt")


open3d()
par(family = "serif")
plot3d(x, y, z + abs(range(z)[1]), xlab = '', ylab = '', zlab = '', aspect = FALSE, col = Colors15[FrazerKarlukAyakulik47CollectionsGroupVec3], size = 0.5, type = 's', axes = TRUE, box = TRUE, top = TRUE, cex = 1)
box3d()
# plot3d(x, y, z + abs(range(z)[1]), aspect = FALSE, col = "black", size = 0.5, type = 'h', box = TRUE, axes = FALSE, top = FALSE, add = TRUE)  # adds pins to spheres
# texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.8, 0.8), text = seq(FrazerKarlukAyakulik47Collections), font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)
texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.2, 0.2), text = FrazerKarlukAyakulik47Collections, font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)
rgl.snapshot("MDS/MDSAdegenetNei47ColFrazerKarlukAyakulik89loci.png", fmt="png", top=TRUE )






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Frazer/Ayakulik Structure ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 96 Loci
## Locus Control
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
loci96 <- dget(file = "Objects/loci96.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Get genotypes

SWKodiakCollections46 <- c("SPINNM08", "SSTUM08", "SCOUR08", "SMIDWS08", "SHOLFS08", "SMIDWM08", "SSUMMM09", "SLINDM08", "SVALA08", "SOUTS08", "SDOGSC08", "SCAIDA14",
                           "SDOGSC11E", "SDOGSC11EJ", "SDOGSC11W", "SDOGSC12E", "SDOGSC12W",
                           "SREDSS11", "SREDSWS11", "SREDWS12", "SREDNWS11", "SREDNES11", "SREDCRY11", "SREDCON11", "SAYAK00", "SAYAK08L",
                           "SAYAK11", "SAYAK12",
                           "SFAL99E", "SCAN99E", "SCAS99E", "SUTHU99E", "SUTHU00E", "SSAL99E", "SHAL01E", "SGRA99E", "SCOT99E", "SMOR99E",
                           "SOMALL99", "SKARLSE11", "SKARLSE99L", "SLTHUM99", "STHUS99L", "SKARLW99L", "SKARLE99L",
                           "SKARL01L")

LOKI2R.GCL(sillyvec = sillyvec, username = "krshedd", password = password)
sapply(SWKodiakCollections46, function(silly) {LOKI2R.GCL(sillyvec = silly, username = "krshedd", password = password)} )  # looping through due to heap space error when getting all sillys at once

rm(password)

str(SSUMMM09.gcl)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Check loci
## Get sample size by locus
Original_SWKodiakCollections46_SampleSizebyLocus <- SampSizeByLocus.GCL(sillyvec = SWKodiakCollections46, loci = loci96)
min(Original_SWKodiakCollections46_SampleSizebyLocus)  ## 52
apply(Original_SWKodiakCollections46_SampleSizebyLocus, 1, min) / apply(Original_SWKodiakCollections46_SampleSizebyLocus, 1, max)  ## Good, 0.91

Original_SWKodiakCollections46_PercentbyLocus <- apply(Original_SWKodiakCollections46_SampleSizebyLocus, 1, function(row) {row / max(row)} )
which(apply(Original_SWKodiakCollections46_PercentbyLocus, 2, min) < 0.8)  # no re-runs!

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_SWKodiakCollections46_PercentbyLocus), 
          col.regions = new.colors, 
          at = seq(from = 0, to = 1, length.out = 100), 
          main = "% Genotyped", xlab = "SILLY", ylab = "Locus", 
          scales = list(x = list(rot = 90)), 
          aspect = "fill")  # aspect = "iso" will make squares


#### Check individuals
## View Histogram of Failure Rate by Strata
invisible(sapply(SWKodiakCollections46, function(mix) {
  my.gcl <- get(paste(mix, ".gcl", sep = ''))
  failure <- apply(my.gcl$scores[, , 1], 1, function(ind) {sum(ind == "0") / length(ind)} )
  hist(x = failure, main = mix, xlab = "Failure Rate", col = 8, xlim = c(0, 1), ylim = c(0, 20), breaks = seq(from = 0, to = 1, by = 0.02))
  abline(v = 0.2, lwd = 3)
}))



SWKodiakCollections46_SampleSizes <- matrix(data = NA, nrow = length(SWKodiakCollections46), ncol = 5, 
                                            dimnames = list(SWKodiakCollections46, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_SWKodiakCollections46_ColSize <- sapply(paste(SWKodiakCollections46, ".gcl", sep = ''), function(x) get(x)$n)
SWKodiakCollections46_SampleSizes[, "Genotyped"] <- Original_SWKodiakCollections46_ColSize


### Alternate
## Indentify alternate species individuals
SWKodiakCollections46_Alternate <- FindAlternateSpecies.GCL(sillyvec = SWKodiakCollections46, species = "sockeye")

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport = SWKodiakCollections46_Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5)

## Get number of individuals per silly after removing alternate species individuals
ColSize_SWKodiakCollections46_PostAlternate <- sapply(paste(SWKodiakCollections46, ".gcl", sep = ''), function(x) get(x)$n)
SWKodiakCollections46_SampleSizes[, "Alternate"] <- Original_SWKodiakCollections46_ColSize-ColSize_SWKodiakCollections46_PostAlternate


### Missing
## Remove individuals with >20% missing data
SWKodiakCollections46_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = SWKodiakCollections46, proportion = 0.8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_SWKodiakCollections46_PostMissLoci <- sapply(paste(SWKodiakCollections46, ".gcl", sep = ''), function(x) get(x)$n)
SWKodiakCollections46_SampleSizes[, "Missing"] <- ColSize_SWKodiakCollections46_PostAlternate-ColSize_SWKodiakCollections46_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
SWKodiakCollections46_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = SWKodiakCollections46, loci = loci96, quantile = NULL, minproportion = 0.95)
SWKodiakCollections46_DuplicateCheckReportSummary <- sapply(SWKodiakCollections46, function(x) SWKodiakCollections46_DuplicateCheck95MinProportion[[x]]$report)
SWKodiakCollections46_DuplicateCheckReportSummary

## Remove duplicate individuals
SWKodiakCollections46_RemovedDups <- RemoveDups.GCL(SWKodiakCollections46_DuplicateCheck95MinProportion)

## Get number of individuals per silly after removing duplicate individuals
ColSize_SWKodiakCollections46_PostDuplicate <- sapply(paste(SWKodiakCollections46, ".gcl", sep = ''), function(x) get(x)$n)
SWKodiakCollections46_SampleSizes[, "Duplicate"] <- ColSize_SWKodiakCollections46_PostMissLoci-ColSize_SWKodiakCollections46_PostDuplicate


### Final
SWKodiakCollections46_SampleSizes[, "Final"] <- ColSize_SWKodiakCollections46_PostDuplicate
SWKodiakCollections46_SampleSizes


### Dput post-QC fish
dir.create(path = "Raw genotypes/PostQCCollectionsSWKodiak46")
invisible(sapply(SWKodiakCollections46, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostQCCollectionsSWKodiak46/" , silly, ".txt", sep = ''))} )); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check HWE in collections (only use diploid loci)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mito.loci <- dget(file = "Objects/mito.loci.txt")

gcl2Genepop.GCL(sillyvec = SWKodiakCollections46, loci = loci96[-mito.loci], path = "Genepop/SWKodiakCollections46_93nuclearloci.txt", VialNums = TRUE); beep(2)

## Check HWE in Genepop
# Originally done with Exact test and subsequently re-done with MCMC
HWE <- ReadGenepopHWE.GCL(file = "Genepop/SWKodiakCollections46_93nuclearloci.txt.P")
str(HWE)


sapply(as.character(unique(HWE$DataByPop$Pop)), function(pop) {
  x <- subset(x = HWE$DataByPop, subset = Pop == pop)
  plot(sort(x[, "WC Fis"]), type = "h", lwd = 5, ylab = "WC Fis", xlab = "Loci (sorted)", col = "grey40", main = pop)
  abline(h = 0, lwd = 5)
  x[x$PValue[!is.na(x$PValue)] < 0.05, ]
}, USE.NAMES = TRUE, simplify = FALSE
)

str(HWE$SummaryPValues)
HWE$SummaryPValues["Overall Loci", ]

sapply(SWKodiakCollections46, function(silly) {
  hist(HWE$SummaryPValues[1:93, grep(pattern = silly, colnames(HWE$SummaryPValues))], breaks = seq(from = 0, to = 1, by = 0.05), col = 8, main = silly, ylim = c(0, 50), xlab = "HWE P-values")
  abline(h = 5, lwd = 3, col = "red")
  })

apply(HWE$SummaryPValues[1:93, seq(SWKodiakCollections46)], 2, function(silly) {sum(silly < 0.05, na.rm = TRUE)})
hist(apply(HWE$SummaryPValues[1:93, seq(SWKodiakCollections46)], 2, function(silly) {sum(silly < 0.05, na.rm = TRUE)}), main = "# loci out of HWE per collection", xlab = "# loci with HWE P-value < 0.05", breaks = 0:14, col = 8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create Genepop for Convert for Structure
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

loci89 <- dget(file = "Objects/loci89.txt")

loci96[!loci96 %in% loci89]
loci2remove <- c("One_ACBP-79", "One_c3-98", "One_GPDH2", "One_MHC2_251", "One_Tf_ex10-750")
loci91 <- loci96[!loci96 %in% loci2remove]
gcl2Genepop.GCL(sillyvec = SWKodiakCollections46, loci = loci91, path = "Genepop/SWKodiakCollections46_loci91.txt", VialNums = TRUE); beep(2)


loci2remove <- c("One_ACBP-79", "One_c3-98", "One_CO1", "One_Cytb_17", "One_Cytb_26", "One_GPDH2", "One_MHC2_251", "One_Tf_ex10-750")
loci88 <- loci96[!loci96 %in% loci2remove]
gcl2Genepop.GCL(sillyvec = SWKodiakCollections46, loci = loci88, path = "Genepop/SWKodiakCollections46_loci88.txt", VialNums = TRUE); beep(2)



# Collection Size
sum(sapply(SWKodiakCollections46, function(silly) {get(paste(silly, ".gcl", sep = ''))$n}))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clean workspace, get post-QC collections
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 96 Loci
## Locus Control
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
loci96 <- dget(file = "Objects/loci96.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Get genotypes

SWKodiakCollections46 <- c("SPINNM08", "SSTUM08", "SCOUR08", "SMIDWS08", "SHOLFS08", "SMIDWM08", "SSUMMM09", "SLINDM08", "SVALA08", "SOUTS08", "SDOGSC08", "SCAIDA14",
                           "SDOGSC11E", "SDOGSC11EJ", "SDOGSC11W", "SDOGSC12E", "SDOGSC12W",
                           "SREDSS11", "SREDSWS11", "SREDWS12", "SREDNWS11", "SREDNES11", "SREDCRY11", "SREDCON11", "SAYAK00", "SAYAK08L",
                           "SAYAK11", "SAYAK12",
                           "SFAL99E", "SCAN99E", "SCAS99E", "SUTHU99E", "SUTHU00E", "SSAL99E", "SHAL01E", "SGRA99E", "SCOT99E", "SMOR99E",
                           "SOMALL99", "SKARLSE11", "SKARLSE99L", "SLTHUM99", "STHUS99L", "SKARLW99L", "SKARLE99L",
                           "SKARL01L")

invisible(sapply(SWKodiakCollections46, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCCollectionsSWKodiak46/", silly, ".txt", sep = "")), pos = 1)} )); beep(4)

# Collection Size
sapply(SWKodiakCollections46, function(silly) {get(paste(silly, ".gcl", sep = ''))$n})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stratify SAYAK12, Ayakulik Weir 2012
table(SAYAK12.gcl$attributes$CAPTURE_DATE)

sort(unique(SAYAK12.gcl$attributes$CAPTURE_DATE))


# June 5-11
AyakulikWeir_June5_11_IDs <- AttributesToIDs.GCL(silly = "SAYAK12", attribute = "CAPTURE_DATE", matching = sort(unique(SAYAK12.gcl$attributes$CAPTURE_DATE))[1:2])
AyakulikWeir_June5_11_IDs <- list(na.omit(AyakulikWeir_June5_11_IDs))
names(AyakulikWeir_June5_11_IDs) <- "SAYAK12"

PoolCollections.GCL("SAYAK12", loci = loci96, IDs = AyakulikWeir_June5_11_IDs, newname = "AyakulikWeir_June5_11")
AyakulikWeir_June5_11.gcl$n


# June 17-25
AyakulikWeir_June17_25_IDs <- AttributesToIDs.GCL(silly = "SAYAK12", attribute = "CAPTURE_DATE", matching = sort(unique(SAYAK12.gcl$attributes$CAPTURE_DATE))[3:4])
AyakulikWeir_June17_25_IDs <- list(na.omit(AyakulikWeir_June17_25_IDs))
names(AyakulikWeir_June17_25_IDs) <- "SAYAK12"

PoolCollections.GCL("SAYAK12", loci = loci96, IDs = AyakulikWeir_June17_25_IDs, newname = "AyakulikWeir_June17_25")
AyakulikWeir_June17_25.gcl$n


# July 2-8
AyakulikWeir_July2_8_IDs <- AttributesToIDs.GCL(silly = "SAYAK12", attribute = "CAPTURE_DATE", matching = sort(unique(SAYAK12.gcl$attributes$CAPTURE_DATE))[5:6])
AyakulikWeir_July2_8_IDs <- list(na.omit(AyakulikWeir_July2_8_IDs))
names(AyakulikWeir_July2_8_IDs) <- "SAYAK12"

PoolCollections.GCL("SAYAK12", loci = loci96, IDs = AyakulikWeir_July2_8_IDs, newname = "AyakulikWeir_July2_8")
AyakulikWeir_July2_8.gcl$n


# July 16-19
AyakulikWeir_July16_19_IDs <- AttributesToIDs.GCL(silly = "SAYAK12", attribute = "CAPTURE_DATE", matching = sort(unique(SAYAK12.gcl$attributes$CAPTURE_DATE))[7:8])
AyakulikWeir_July16_19_IDs <- list(na.omit(AyakulikWeir_July16_19_IDs))
names(AyakulikWeir_July16_19_IDs) <- "SAYAK12"

PoolCollections.GCL("SAYAK12", loci = loci96, IDs = AyakulikWeir_July16_19_IDs, newname = "AyakulikWeir_July16_19")
AyakulikWeir_July16_19.gcl$n


# July 26-August 2
AyakulikWeir_July26_August2_IDs <- AttributesToIDs.GCL(silly = "SAYAK12", attribute = "CAPTURE_DATE", matching = sort(unique(SAYAK12.gcl$attributes$CAPTURE_DATE))[9:10])
AyakulikWeir_July26_August2_IDs <- list(na.omit(AyakulikWeir_July26_August2_IDs))
names(AyakulikWeir_July26_August2_IDs) <- "SAYAK12"

PoolCollections.GCL("SAYAK12", loci = loci96, IDs = AyakulikWeir_July26_August2_IDs, newname = "AyakulikWeir_July26_August2")
AyakulikWeir_July26_August2.gcl$n


# August 9-19
AyakulikWeir_August9_19_IDs <- AttributesToIDs.GCL(silly = "SAYAK12", attribute = "CAPTURE_DATE", matching = sort(unique(SAYAK12.gcl$attributes$CAPTURE_DATE))[11:12])
AyakulikWeir_August9_19_IDs <- list(na.omit(AyakulikWeir_August9_19_IDs))
names(AyakulikWeir_August9_19_IDs) <- "SAYAK12"

PoolCollections.GCL("SAYAK12", loci = loci96, IDs = AyakulikWeir_August9_19_IDs, newname = "AyakulikWeir_August9_19")
AyakulikWeir_August9_19.gcl$n


# August 23
AyakulikWeir_August23_IDs <- AttributesToIDs.GCL(silly = "SAYAK12", attribute = "CAPTURE_DATE", matching = sort(unique(SAYAK12.gcl$attributes$CAPTURE_DATE))[13])
AyakulikWeir_August23_IDs <- list(na.omit(AyakulikWeir_August23_IDs))
names(AyakulikWeir_August23_IDs) <- "SAYAK12"

PoolCollections.GCL("SAYAK12", loci = loci96, IDs = AyakulikWeir_August23_IDs, newname = "AyakulikWeir_August23")
AyakulikWeir_August23.gcl$n


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FreqPlot


#### All Ayakulik

Ayakulik14Collections <- c(SWKodiakCollections46[23:24], 
                         "AyakulikWeir_June5_11",
                         "AyakulikWeir_June17_25",
                         "AyakulikWeir_July2_8",
                         "AyakulikWeir_July16_19",
                         "AyakulikWeir_July26_August2",
                         "AyakulikWeir_August9_19",
                         "AyakulikWeir_August23",
                         SWKodiakCollections46[18:22])


sampsizes <- sapply(Ayakulik14Collections, function(silly) {get(paste(silly, ".gcl", sep = ''))$n})


# Counts
Ayakulik14CollectionsPostRemovalsAlleleCounts <- FreqPop.GCL(sillyvec = Ayakulik14Collections, loci = loci96)
str(Ayakulik14CollectionsPostRemovalsAlleleCounts)

# Frequencies
Ayakulik14CollectionsPostRemovalsFreqs <- Ayakulik14CollectionsPostRemovalsAlleleCounts[,,"Allele 1"] / (Ayakulik14CollectionsPostRemovalsAlleleCounts[,,"Allele 2"] + Ayakulik14CollectionsPostRemovalsAlleleCounts[,,"Allele 1"])
str(Ayakulik14CollectionsPostRemovalsFreqs)

# 95% CI
Ayakulik14CollectionsPostRemovalsFreqs95CI
sampsizebylocus <- apply(Ayakulik14CollectionsPostRemovalsAlleleCounts, c(1, 2), sum)

Ayakulik14CollectionsPostRemovalsFreqs5CI <- qbinom(p = 0.05, size = sampsizebylocus, prob = Ayakulik14CollectionsPostRemovalsFreqs) / sampsizebylocus
Ayakulik14CollectionsPostRemovalsFreqs95CI <- qbinom(p = 0.95, size = sampsizebylocus, prob = Ayakulik14CollectionsPostRemovalsFreqs) / sampsizebylocus



#### Create plots
## Create groupvec
Ayakulik14CollectionsGroupVec3 <- c(rep(1, 2), rep(2, 7), rep(3, 5))


## Plot points and smoothing line
getwd()
pdf(file = "FreqPlots/Ayakulik14Collections_96SNPs_GroupVec2_FreqPlots.pdf",width=11, height=8.5, family="Times",pointsize=20)   
par(mar = c(2.1, 4.1, 4.1, 4.6))
for(locus in loci96){
  plot(Ayakulik14CollectionsPostRemovalsFreqs[, locus], main = locus, col = c("blue", "skyblue", "blue4")[Ayakulik14CollectionsGroupVec3], pch = 19, ylim = c(0, 1), ylab = "Frequency", xlab = "Collections", xaxt = 'n', cex = 1)
  arrows(x0 = seq(Ayakulik14Collections), x1 = seq(Ayakulik14Collections), y0 = Ayakulik14CollectionsPostRemovalsFreqs5CI[, locus], y1 = Ayakulik14CollectionsPostRemovalsFreqs95CI[, locus], angle = 90, code = 3, lwd = 1.5, length = 0.1)
  lines(supsmu(seq(length(Ayakulik14Collections)), Ayakulik14CollectionsPostRemovalsFreqs[, locus]), lwd = 2)
  points(Ayakulik14CollectionsPostRemovalsFreqs[, locus], col = c("blue", "skyblue", "blue4")[Ayakulik14CollectionsGroupVec3], pch = 19, cex = 1)
  mtext(text = "Collections", side = 1, line = 0.7)
  text(x = seq(Ayakulik14Collections), y = Ayakulik14CollectionsPostRemovalsFreqs[, locus], labels = seq(Ayakulik14Collections), cex = 0.4, col = "white")
  text(x = seq(Ayakulik14Collections), y = -0.07, labels = sampsizes, cex = 0.6, xpd = TRUE)
  legend(x = length(Ayakulik14CollectionsGroupVec3) + length(Ayakulik14CollectionsGroupVec3)/40, y = 0.9, legend = c("Early Baseline", "2012 Weir", "Late Baseline"), fill = c("blue", "skyblue", "blue4"), xpd = TRUE, bty = "n", cex = 0.7)
}; rm(locus)
dev.off()

t(t(Ayakulik14Collections))









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stratify SDOGSC12, Dog Salmon Weir 2012
table(SDOGSC12E.gcl$attributes$CAPTURE_DATE)
table(SDOGSC12W.gcl$attributes$CAPTURE_DATE)

PoolCollections.GCL(collections = c("SDOGSC12E", "SDOGSC12W"), loci = loci96, newname = "SDOGSC12")
str(SDOGSC12.gcl)
table(SDOGSC12.gcl$attributes$CAPTURE_DATE)


# June 25-26
DogSalmonWeir_June25_26_IDs <- AttributesToIDs.GCL(silly = "SDOGSC12", attribute = "CAPTURE_DATE", matching = sort(unique(SDOGSC12.gcl$attributes$CAPTURE_DATE))[1:2])
DogSalmonWeir_June25_26_IDs <- list(na.omit(DogSalmonWeir_June25_26_IDs))
names(DogSalmonWeir_June25_26_IDs) <- "SDOGSC12"

PoolCollections.GCL("SDOGSC12", loci = loci96, IDs = DogSalmonWeir_June25_26_IDs, newname = "DogSalmonWeir_June25_26")
DogSalmonWeir_June25_26.gcl$n


# July 5
DogSalmonWeir_July5_IDs <- AttributesToIDs.GCL(silly = "SDOGSC12", attribute = "CAPTURE_DATE", matching = sort(unique(SDOGSC12.gcl$attributes$CAPTURE_DATE))[3])
DogSalmonWeir_July5_IDs <- list(na.omit(DogSalmonWeir_July5_IDs))
names(DogSalmonWeir_July5_IDs) <- "SDOGSC12"

PoolCollections.GCL("SDOGSC12", loci = loci96, IDs = DogSalmonWeir_July5_IDs, newname = "DogSalmonWeir_July5")
DogSalmonWeir_July5.gcl$n


# July 10-11
DogSalmonWeir_July10_11_IDs <- AttributesToIDs.GCL(silly = "SDOGSC12", attribute = "CAPTURE_DATE", matching = sort(unique(SDOGSC12.gcl$attributes$CAPTURE_DATE))[4:5])
DogSalmonWeir_July10_11_IDs <- list(na.omit(DogSalmonWeir_July10_11_IDs))
names(DogSalmonWeir_July10_11_IDs) <- "SDOGSC12"

PoolCollections.GCL("SDOGSC12", loci = loci96, IDs = DogSalmonWeir_July10_11_IDs, newname = "DogSalmonWeir_July10_11")
DogSalmonWeir_July10_11.gcl$n


# July 16-17
DogSalmonWeir_July16_17_IDs <- AttributesToIDs.GCL(silly = "SDOGSC12", attribute = "CAPTURE_DATE", matching = sort(unique(SDOGSC12.gcl$attributes$CAPTURE_DATE))[6:7])
DogSalmonWeir_July16_17_IDs <- list(na.omit(DogSalmonWeir_July16_17_IDs))
names(DogSalmonWeir_July16_17_IDs) <- "SDOGSC12"

PoolCollections.GCL("SDOGSC12", loci = loci96, IDs = DogSalmonWeir_July16_17_IDs, newname = "DogSalmonWeir_July16_17")
DogSalmonWeir_July16_17.gcl$n


# July 23-24
DogSalmonWeir_July23_24_IDs <- AttributesToIDs.GCL(silly = "SDOGSC12", attribute = "CAPTURE_DATE", matching = sort(unique(SDOGSC12.gcl$attributes$CAPTURE_DATE))[8:9])
DogSalmonWeir_July23_24_IDs <- list(na.omit(DogSalmonWeir_July23_24_IDs))
names(DogSalmonWeir_July23_24_IDs) <- "SDOGSC12"

PoolCollections.GCL("SDOGSC12", loci = loci96, IDs = DogSalmonWeir_July23_24_IDs, newname = "DogSalmonWeir_July23_24")
DogSalmonWeir_July23_24.gcl$n


# July 30-August 7
DogSalmonWeir_July30_August7_IDs <- AttributesToIDs.GCL(silly = "SDOGSC12", attribute = "CAPTURE_DATE", matching = sort(unique(SDOGSC12.gcl$attributes$CAPTURE_DATE))[10:11])
DogSalmonWeir_July30_August7_IDs <- list(na.omit(DogSalmonWeir_July30_August7_IDs))
names(DogSalmonWeir_July30_August7_IDs) <- "SDOGSC12"

PoolCollections.GCL("SDOGSC12", loci = loci96, IDs = DogSalmonWeir_July30_August7_IDs, newname = "DogSalmonWeir_July30_August7")
DogSalmonWeir_July30_August7.gcl$n


# August 13-17
DogSalmonWeir_August13_17_IDs <- AttributesToIDs.GCL(silly = "SDOGSC12", attribute = "CAPTURE_DATE", matching = sort(unique(SDOGSC12.gcl$attributes$CAPTURE_DATE))[12])
DogSalmonWeir_August13_17_IDs <- list(na.omit(DogSalmonWeir_August13_17_IDs))
names(DogSalmonWeir_August13_17_IDs) <- "SDOGSC12"

PoolCollections.GCL("SDOGSC12", loci = loci96, IDs = DogSalmonWeir_August13_17_IDs, newname = "DogSalmonWeir_August13_17")
DogSalmonWeir_August13_17.gcl$n


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FreqPlot


#### All Frazer

Frazer19Collections <- c(SWKodiakCollections46[1:12], 
                         "DogSalmonWeir_June25_26",
                         "DogSalmonWeir_July5",
                         "DogSalmonWeir_July10_11",
                         "DogSalmonWeir_July16_17",
                         "DogSalmonWeir_July23_24",
                         "DogSalmonWeir_July30_August7",
                         "DogSalmonWeir_August13_17")
                         

sampsizes <- sapply(Frazer19Collections, function(silly) {get(paste(silly, ".gcl", sep = ''))$n})


# Counts
Frazer19CollectionsPostRemovalsAlleleCounts <- FreqPop.GCL(sillyvec = Frazer19Collections, loci = loci96)
str(Frazer19CollectionsPostRemovalsAlleleCounts)

# Frequencies
Frazer19CollectionsPostRemovalsFreqs <- Frazer19CollectionsPostRemovalsAlleleCounts[,,"Allele 1"] / (Frazer19CollectionsPostRemovalsAlleleCounts[,,"Allele 2"] + Frazer19CollectionsPostRemovalsAlleleCounts[,,"Allele 1"])
str(Frazer19CollectionsPostRemovalsFreqs)

# 95% CI
Frazer19CollectionsPostRemovalsFreqs95CI
sampsizebylocus <- apply(Frazer19CollectionsPostRemovalsAlleleCounts, c(1, 2), sum)

Frazer19CollectionsPostRemovalsFreqs5CI <- qbinom(p = 0.05, size = sampsizebylocus, prob = Frazer19CollectionsPostRemovalsFreqs) / sampsizebylocus
Frazer19CollectionsPostRemovalsFreqs95CI <- qbinom(p = 0.95, size = sampsizebylocus, prob = Frazer19CollectionsPostRemovalsFreqs) / sampsizebylocus



#### Create plots
## Create groupvec
Frazer19CollectionsGroupVec2 <- c(rep(1, 12), rep(2, 7))


## Plot points and smoothing line
getwd()
pdf(file = "FreqPlots/Frazer19Collections_96SNPs_GroupVec2_FreqPlots.pdf",width=11, height=8.5, family="Times",pointsize=20)   
par(mar = c(2.1, 4.1, 4.1, 4.6))
for(locus in loci96){
  plot(Frazer19CollectionsPostRemovalsFreqs[, locus], main = locus, col = c("cyan", "skyblue")[Frazer19CollectionsGroupVec2], pch = 19, ylim = c(0, 1), ylab = "Frequency", xlab = "Collections", xaxt = 'n', cex = 1)
  arrows(x0 = seq(Frazer19Collections), x1 = seq(Frazer19Collections), y0 = Frazer19CollectionsPostRemovalsFreqs5CI[, locus], y1 = Frazer19CollectionsPostRemovalsFreqs95CI[, locus], angle = 90, code = 3, lwd = 1.5, length = 0.1)
  lines(supsmu(seq(length(Frazer19Collections)), Frazer19CollectionsPostRemovalsFreqs[, locus]), lwd = 2)
  points(Frazer19CollectionsPostRemovalsFreqs[, locus], col = c("cyan", "skyblue")[Frazer19CollectionsGroupVec2], pch = 19, cex = 1)
  mtext(text = "Collections", side = 1, line = 0.7)
  text(x = seq(Frazer19Collections), y = Frazer19CollectionsPostRemovalsFreqs[, locus], labels = seq(Frazer19Collections), cex = 0.4)
  text(x = seq(Frazer19Collections), y = -0.07, labels = sampsizes, cex = 0.6, xpd = TRUE)
  legend(x = length(Frazer19CollectionsGroupVec2) + length(Frazer19CollectionsGroupVec2)/20, y = 0.9, legend = c("Baseline", "2012 Weir"), fill = c("cyan", "skyblue"), xpd = TRUE, bty = "n", cex = 0.7)
}; rm(locus)
dev.off()

t(t(Frazer19Collections))









#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FreqPlot


#### All Frazer and Ayakulik

FrazerAyakulik33Collections <- c(Frazer19Collections, Ayakulik14Collections)


sampsizes <- sapply(FrazerAyakulik33Collections, function(silly) {get(paste(silly, ".gcl", sep = ''))$n})


# Counts
FrazerAyakulik33CollectionsPostRemovalsAlleleCounts <- FreqPop.GCL(sillyvec = FrazerAyakulik33Collections, loci = loci96)
str(FrazerAyakulik33CollectionsPostRemovalsAlleleCounts)

# Frequencies
FrazerAyakulik33CollectionsPostRemovalsFreqs <- FrazerAyakulik33CollectionsPostRemovalsAlleleCounts[,,"Allele 1"] / (FrazerAyakulik33CollectionsPostRemovalsAlleleCounts[,,"Allele 2"] + FrazerAyakulik33CollectionsPostRemovalsAlleleCounts[,,"Allele 1"])
str(FrazerAyakulik33CollectionsPostRemovalsFreqs)

# 95% CI
FrazerAyakulik33CollectionsPostRemovalsFreqs95CI
sampsizebylocus <- apply(FrazerAyakulik33CollectionsPostRemovalsAlleleCounts, c(1, 2), sum)

FrazerAyakulik33CollectionsPostRemovalsFreqs5CI <- qbinom(p = 0.05, size = sampsizebylocus, prob = FrazerAyakulik33CollectionsPostRemovalsFreqs) / sampsizebylocus
FrazerAyakulik33CollectionsPostRemovalsFreqs95CI <- qbinom(p = 0.95, size = sampsizebylocus, prob = FrazerAyakulik33CollectionsPostRemovalsFreqs) / sampsizebylocus



#### Create plots
## Create groupvec
FrazerAyakulik33CollectionsGroupVec3 <- c(Frazer19CollectionsGroupVec2, Ayakulik14CollectionsGroupVec3 + 2)
colors5 <- c("cyan", "skyblue3", "blue", "skyblue", "blue4")

## Plot points and smoothing line
getwd()
pdf(file = "FreqPlots/FrazerAyakulik33Collections_96SNPs_GroupVec2_FreqPlots.pdf",width=11, height=8.5, family="Times",pointsize=20)   
par(mar = c(2.1, 4.1, 4.1, 4.6))
for(locus in loci96){
  plot(FrazerAyakulik33CollectionsPostRemovalsFreqs[, locus], main = locus, col = colors5[FrazerAyakulik33CollectionsGroupVec3], pch = 19, ylim = c(0, 1), ylab = "Frequency", xlab = "Collections", xaxt = 'n', cex = 1)
  arrows(x0 = seq(FrazerAyakulik33Collections), x1 = seq(FrazerAyakulik33Collections), y0 = FrazerAyakulik33CollectionsPostRemovalsFreqs5CI[, locus], y1 = FrazerAyakulik33CollectionsPostRemovalsFreqs95CI[, locus], angle = 90, code = 3, lwd = 1.5, length = 0.1)
  lines(supsmu(seq(length(FrazerAyakulik33Collections)), FrazerAyakulik33CollectionsPostRemovalsFreqs[, locus]), lwd = 2)
  points(FrazerAyakulik33CollectionsPostRemovalsFreqs[, locus], col = colors5[FrazerAyakulik33CollectionsGroupVec3], pch = 19, cex = 1)
  mtext(text = "Collections", side = 1, line = 0.7)
  text(x = seq(FrazerAyakulik33Collections), y = FrazerAyakulik33CollectionsPostRemovalsFreqs[, locus], labels = seq(FrazerAyakulik33Collections), cex = 0.4, col = "white")
  text(x = seq(FrazerAyakulik33Collections), y = -0.07, labels = sampsizes, cex = 0.5, xpd = TRUE)
  legend(x = length(FrazerAyakulik33CollectionsGroupVec3) + length(FrazerAyakulik33CollectionsGroupVec3)/40, y = 0.9, legend = c("Frazer Baseline", "Frazer Weir", "Ayakulik Early", "Ayakulik Weir", "Ayakulik Late"), fill = colors5, xpd = TRUE, bty = "n", cex = 0.6)
}; rm(locus)
dev.off()

t(t(FrazerAyakulik33Collections))








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Combine Loci
CombineLoci.GCL(sillyvec = FrazerAyakulik33Collections, markerset = c("One_CO1", "One_Cytb_17", "One_Cytb_26"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = FrazerAyakulik33Collections, markerset = c("One_Tf_ex10-750", "One_Tf_ex3-182"), update = TRUE); beep(2)

loci89 <- dget(file = "Objects/loci89.txt")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MDS

gcl2Genepop.GCL(sillyvec = FrazerAyakulik33Collections, loci = loci89, path = "Genepop/FrazerAyakulik33Collections_89loci.gen", VialNums = TRUE)

detach("package:adegenet", unload = TRUE)
require(package = adegenet, lib.loc = "C:/Users/krshedd/Documents/R/win-library/3.1")

genind <- read.genepop(file = "Genepop/FrazerAyakulik33Collections_89loci.gen")

genpop <- genind2genpop(genind)

# AdegenetNei33Col89loci <- dist.genpop(genpop, method = 1, diag = TRUE, upper = TRUE)
# dput(x = AdegenetNei33Col89loci, file = "Trees/AdegenetNei33Col89loci.txt")
AdegenetNei33Col89loci <- dget(file = "Trees/AdegenetNei33Col89loci.txt")
str(AdegenetNei33Col89loci)

require(ape)
Nei33NJtree <- nj(AdegenetNei33Col89loci)
str(Nei33NJtree)
Nei33NJtree$tip.label <- readClipboard()
plot.phylo(x = Nei33NJtree, cex = 0.5, no.margin = TRUE)

library('rgl')

MDS <- cmdscale(as.matrix(AdegenetNei33Col89loci), k = 3)
# MDS <- cmdscale(as.matrix(AdegenetNei33Col89loci), k = 40, eig = TRUE)  # Do with all possible dimensions
# dput(x = MDS, file = "Objects/MDSAdegenetNei33Col89loci.txt")
# dput(x = MDS, file = "Objects/MDSAdegenetNei33Col89loci_alldim.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei33Col89loci.txt")

x <- as.vector(MDS[, 1])   
y <- as.vector(MDS[, 2])
z <- as.vector(MDS[, 3])

FrazerAyakulik33CollectionsGroupVec5 <- c(Frazer19CollectionsGroupVec2, Ayakulik14CollectionsGroupVec3 + 2)
colors5 <- c("cyan", "red", "blue", "green", "blue4")


open3d()
par(family = "serif")
plot3d(x, y, z + abs(range(z)[1]), xlab = '', ylab = '', zlab = '', aspect = FALSE, col = colors5[FrazerAyakulik33CollectionsGroupVec5], size = 0.5, type = 's', axes = TRUE, box = TRUE, top = TRUE, cex = 1)
box3d()
# plot3d(x, y, z + abs(range(z)[1]), aspect = FALSE, col = "black", size = 0.5, type = 'h', box = TRUE, axes = FALSE, top = FALSE, add = TRUE)  # adds pins to spheres
# texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.8, 0.8), text = seq(FrazerAyakulik33Collections), font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)
texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.2, 0.2), text = FrazerAyakulik33Collections, font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)

rgl.snapshot("MDS/MDSAdegenetNei33ColFrazerKarlukAyakulik89loci2.png", fmt="png", top=TRUE )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Likelihood profile

FrazerAyakulik33CollectionsGroupVec5 <- c(Frazer19CollectionsGroupVec2, Ayakulik14CollectionsGroupVec3 + 2)
FrazerAyakulik33CollectionsGroupVec2 <- c(rep(1, length(Frazer19CollectionsGroupVec2)), rep(2, length(Ayakulik14CollectionsGroupVec3)))

FrazerAyakulik33CollectionsLeaveOneOutLikeProfile <- LeaveOneOutLikeProfile.GCL(popvec = FrazerAyakulik33Collections, loci = loci89, groupvec = FrazerAyakulik33CollectionsGroupVec2, groupnames = c("Frazer", "Ayakulik"), groupcomps = NULL, ncores = 5)
str(FrazerAyakulik33CollectionsLeaveOneOutLikeProfile)

PlotLikeProfile.GCL(likeprof = FrazerAyakulik33CollectionsLeaveOneOutLikeProfile, popvec = FrazerAyakulik33Collections, loci = loci89, groupvec = FrazerAyakulik33CollectionsGroupVec2, groupnames = c("Frazer", "Ayakulik"), dir = "Likelihood Profiles", filename = "FrazerAyakulik33CollectionsLeaveOneOutLikeProfile.pdf", col = c("cyan", "blue"))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Escapement test for Frazer ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 96 Loci
## Locus Control
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
loci96 <- dget(file = "Objects/loci96.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Get genotypes

SWKodiakCollections46 <- c("SPINNM08", "SSTUM08", "SCOUR08", "SMIDWS08", "SHOLFS08", "SMIDWM08", "SSUMMM09", "SLINDM08", "SVALA08", "SOUTS08", "SDOGSC08", "SCAIDA14",
                           "SDOGSC11E", "SDOGSC11EJ", "SDOGSC11W", "SDOGSC12E", "SDOGSC12W",
                           "SREDSS11", "SREDSWS11", "SREDWS12", "SREDNWS11", "SREDNES11", "SREDCRY11", "SREDCON11", "SAYAK00", "SAYAK08L",
                           "SAYAK11", "SAYAK12",
                           "SFAL99E", "SCAN99E", "SCAS99E", "SUTHU99E", "SUTHU00E", "SSAL99E", "SHAL01E", "SGRA99E", "SCOT99E", "SMOR99E",
                           "SOMALL99", "SKARLSE11", "SKARLSE99L", "SLTHUM99", "STHUS99L", "SKARLW99L", "SKARLE99L",
                           "SKARL01L")

invisible(sapply(SWKodiakCollections46, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCCollectionsSWKodiak46/", silly, ".txt", sep = "")), pos = 1)} )); beep(4)


## FreqFis Plots
SWKodiakCollections46GroupVec8 <- c(rep(1, 12), rep(2, 5), rep(3, 5), rep(4, 2), rep(5, 4), rep(6, 10), rep(7, 7), 8)
FreqFisPlot_SWKodiakCollections46GroupVec8 <- FreqFisPlot4SNPs.GCL(sillyvec = SWKodiakCollections46, loci = loci96, groupvec = SWKodiakCollections46GroupVec8, alpha = 0.05, groupcol = c("cyan", "red", "blue", "blue4", "green", "purple", "grey", "orange"), file = "HWE/FreqFisPlot_SWKodiakCollections46GroupVec8.pdf")
# sillyvec = SWKodiakCollections46; loci = loci96; groupvec = SWKodiakCollections46GroupVec8; alpha = 0.05; groupcol = c("cyan", "red", "blue", "blue4", "green", "purple", "grey", "orange"); file = "HWE/FreqFisPlot_SWKodiakCollections46GroupVec8.pdf"
str(FreqFisPlot_SWKodiakCollections46GroupVec8)

## Pool
PoolCollections.GCL(collections = c("SDOGSC12E", "SDOGSC12W"), loci = loci96, newname = "SDOGSC12")
str(SDOGSC12.gcl)


## Combine Loci
CombineLoci.GCL(sillyvec = "SDOGSC12", markerset = c("One_CO1", "One_Cytb_17", "One_Cytb_26"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = "SDOGSC12", markerset = c("One_CO1", "One_Cytb_26"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = "SDOGSC12", markerset = c("One_Tf_ex10-750", "One_Tf_ex3-182"), update = TRUE); beep(2)

loci89 <- dget(file = "Objects/loci89.txt")
loci46 <- dget(file = "Objects/loci46.txt")




# Starting with a Regionally flat prior, then rolling prior afterwards
KMA473Pops15FlatPrior <- dget(file = "Objects/KMA473Pops15FlatPrior.txt")
KMA473Pops <- dget(file = "Objects/KMA473Pops.txt")
KMA473PopsGroupVec15 <- dget(file = "Objects/KMA473PopsGroupVec15.txt")
KMA473PopsInits <- dget(file = "Objects/KMA473PopsInits.txt")
KMA47346Baseline <- dget(file = "Objects/KMA47346Baseline.txt")
WASSIPSockeyeSeeds <- dget(file = "v:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects/WASSIPSockeyeSeeds.txt")


## Dumping Mixture files
dir.create(path = "BAYES/Escapement Tests")
dir.create(path = "BAYES/Escapement Tests/Control")
dir.create(path = "BAYES/Escapement Tests/Mixture")
dir.create(path = "BAYES/Escapement Tests/Output")


KMA47346MixtureFormat <- CreateMixture.GCL(sillys = "SDOGSC12", loci = loci46, IDs = NULL, mixname = "SDOGSC12",
                                           dir = "BAYES/Escapement Tests/Mixture", type = "BAYES", PT = FALSE)

## Dumping Control files
CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = "SDOGSC12", basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                      groupvec = KMA473PopsGroupVec15, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits, dir = "BAYES/Escapement Tests/Control",
                      seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")

## Create output directory
dir.create(path = "BAYES/Escapement Tests/Output/SDOGSC12")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Read in output
KMA15GroupsPC <- dget(file = "Objects/PCGroups15.txt")
KMA15GroupsPC2Rows <- dget(file = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects/KMA15GroupsPC2Rows.txt")


SDOGSC12_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC ,
  maindir = "BAYES/Escapement Tests/Output/", 
  mixvec = "SDOGSC12", prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(SDOGSC12_Estimates)

QuickBarplot <- dget(file = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects/QuickBarplot.txt")
QuickBarplot(mixvec = "SDOGSC12", estimatesstats = SDOGSC12_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = setNames(object = "Dog Salmon Weir 2012", nm = "SDOGSC12"))
# Escapement test works well, but this is likely due to the "bayesian" pull effect...
