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
KMAobjects <- KMAobjects[!KMAobjects %in% c("EASSIP-LateAugust2014", "OriginalLocusControl96.txt", "OriginalLocusControl48.txt")]
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
KMAobjects <- KMAobjects[!KMAobjects %in% c("EASSIP-LateAugust2014", "OriginalLocusControl96.txt", "OriginalLocusControl48.txt", "LocusControl98.txt", "LocusControl99.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)


## Get post-QC, stratified, combined loci, mixtures
invisible(sapply(KMA2014_2015Strata, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollections_Strata_PostQC_CombinedLoci/", silly, ".txt", sep = "")), pos = 1)} )); beep(4)
objects(pattern = "\\.gcl")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get MSA Objects ####
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
# WASSIPSockeyeSeeds <- dget("V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt")
# KMA473CommonNames <- dget(file = "CommonNames473.txt")
# KMA15Colors <- dget(file = "Colors15.txt")
# setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# 
# ## Copy these baseline objects and put them in the Mixtures/Objects directory
# 
# setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/Objects")
# file.copy(from = c("KMA473Pops15FlatPrior.txt", "KMA473PopsInits.txt", "KMA473PopsGroupVec15.txt", "KMA473Pops.txt", "PCGroups15.txt", "KMA47346Baseline.txt", "KMA47389Baseline.txt", "CommonNames473.txt", "Colors15.txt",
#                    "V:/Analysis/5_Coastwide/Sockeye/WASSIP/Mixture/Objects/WASSIPSockeyeSeeds.txt"), 
#           to = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects")
# setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# 
# 
# KMA15GroupsPC2Rows <- PCGroups15Rows2 <- c("West of\nChignik", "Black\nLake", "Chignik\nLake", "U. Station\nAkalura", "Frazer\n", "Ayakulik\n", "Karluk\n", "Uganik\n", "Northwest\nKodiak", "Afognak\n", "Eastside\nKodiak", "Saltery\n", "Cook\nInlet", "PWS\n", "South of\nCape Suckling")
# dput(x = KMA15GroupsPC2Rows, file = "Objects/KMA15GroupsPC2Rows.txt")
# # Note, I changed the names for CommonNames473.txt -> KMA473CommonNames.txt & PCGroups15 -> KMA15GroupsPC.txt

dir.create(path = "BAYES/2014-2015 Mixtures 46loci")
sapply(c("Control", "Mixture", "Output"), function(folder) {dir.create(path = paste(getwd(), "BAYES/2014-2015 Mixtures 46loci", folder, sep = "/"))} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 MSA files for BAYES 2014 Early Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Starting with a Regionally flat prior, then rolling prior afterwards
KMA473Pops15FlatPrior


KMA2014Strata_1_Early <- grep(pattern = "1_Early", x = KMA2014Strata, value = TRUE)
Round1Mixtures_2014 <- c(KMA2014Strata_1_Early, "SALITC14_2_Middle")
dput(x = Round1Mixtures_2014, file = "Objects/Round1Mixtures_2014.txt")

## Dumping Mixture files
KMA47346MixtureFormat <- CreateMixture.GCL(sillys = KMA2014Strata_1_Early[1], loci = loci46, IDs = NULL, mixname = KMA2014Strata_1_Early[1],
                                           dir = "BAYES/2014-2015 Mixtures 46loci/Mixture", type = "BAYES", PT = FALSE)
dput(KMA47346MixtureFormat, file = "Objects/KMA47346MixtureFormat.txt")

sapply(Round1Mixtures_2014, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round1Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round1Mixtures_2014, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                              maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
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
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
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


PlotPosterior <- function(mixvec = NULL, output, header, groups, colors = NULL, set.mfrow, thin = 10){
  if(is.null(colors)) {colors <- rep("black", length(groups))}
  if(is.null(mixvec)) {mixvec <- names(output)}
  
  par(mfrow = set.mfrow, mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 1.1))
  
  invisible(sapply(mixvec, function(Mix) {
    invisible(sapply(seq(groups), function(i) {
      RG <- output[[Mix]][, i]
      RG <- RG[seq(1, length(RG), thin)]
      plot(RG, type = "l", ylim = c(0,1), xlab = "", ylab = "")
      abline(v = seq(0, length(RG), length(RG)/5), xpd = FALSE)
      text(x = length(RG)/2, y = 0.96, labels = groups[i], col = colors[i], cex = 1.2, font = 2)} ))
    mtext(text = header[Mix], side = 3, outer = TRUE, cex = 1.5)
    mtext(text = "Iteration (5 chains each)", side = 1, outer = TRUE, cex = 1.5, line = 1.5)
    mtext(text = "Posterior", side = 2, outer = TRUE, cex = 1.5, line = 1.5)
  }))
  
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))
}
dput(x = PlotPosterior, file = "Objects/PlotPosterior.txt")

PlotPosterior(mixvec = Round1Mixtures_2014[c(3, 4, 2, 1, 5)], output = Round1Mixtures_2014_Estimates$Output, 
              groups = KMA15GroupsPC, colors = KMA15Colors, 
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

QuickPlot(mixvec = Round1Mixtures_2014, estimatesstats = Round1Mixtures_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round1Mixtures_2014_Header)


## Barplots
QuickBarplot <- function(mixvec, estimatesstats, groups, groups2rows = NULL, header) {
  while(!require(gplots, quietly = TRUE)){install.packages("vioplot")}
  
  # if(is.null(groups2rows)) {groups2rows <- groups}
  if(length(mixvec) != length(header)) {stop("Header is not the same length as mixvec!")}
  for(i in seq(mixvec)) {
    header[i] <- paste(header[i], "\nn = ", get(paste(mixvec[i], ".gcl", sep = ''))$n, sep = '')
  }
  if("Stats" %in% names(estimatesstats)) {estimatesstats <- estimatesstats$Stats}
  
  par(mfrow = c(1, 1), mar = c(3.1, 5.1, 4.1, 2.1), oma = rep(0, 4))
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
      text(x = Barplot[, 1], y = -1, labels = groups, srt = 90, adj =  1, xpd = TRUE, cex = 0.5)
    } else {
      mtext(text = groups2rows, side = 1, line = 1, at = Barplot[, 1], adj = 0.5, cex = 0.6)
    }
  })
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}
dput(x = QuickBarplot, file = "Objects/QuickBarplot.txt")
QuickBarplot(mixvec = Round1Mixtures_2014, estimatesstats = Round1Mixtures_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round1Mixtures_2014_Header)


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


ViolinPlot(estimates = Round1Mixtures_2014_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round1Mixtures_2014_Header)
rm(Round1Mixtures_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Looking into Alitak Middle Within Chain Convergence Issues ####
# Clearly, Frazer and Ayakulik reporting groups fail to converge on a stable mode
# The distribution is bimodal, as the posterior jumps around
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Pulling entire posterior with no burn-in
Round1Alitak_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                              maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
                                                              mixvec = "SALITC14_2_Middle", prior = "",  
                                                              ext = "RGN", nchains = 5, burn = 0, alpha = 0.1, PosteriorOutput = TRUE)

# Plotting entire posterior with no burn-in
nchains <- 5
RG <- "Frazer"
# RG <- "Ayakulik"
par(mfrow = c(5, 1), mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 0))
sapply(seq(nchains), function(chain) {
  posterior.length <- dim(Round1Alitak_2014_Estimates$Output$SALITC14_2_Middle)[1] / nchains
  plot(Round1Alitak_2014_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA15GroupsPC == RG)], 
       type = "l", ylim = c(0, 1), xlab = '', ylab = '')
  abline(v = posterior.length / 2, lwd = 2)
  text(x = posterior.length / 4, y = 1, labels = "Burn-in", pos = 1)
  text(x = posterior.length / 4 * 3, y = 1, labels = "Posterior", pos = 1)} )
mtext(text = paste("Alitak 2014 Middle Strata: ", RG, sep = ''), side = 3, outer = TRUE, cex = 1.5)
mtext(text = "Iteration (5 chains each)", side = 1, outer = TRUE, cex = 1.5, line = 1.5)
mtext(text = "Posterior", side = 2, outer = TRUE, cex = 1.5, line = 1.5)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Looking at different ways to measure MCMC convergence, both within and between chains
# Following methods from <http://www.people.fas.harvard.edu/~plam/teaching/methods/convergence/convergence_print.pdf>
# All done with the "coda" package; code is taken from CustomCombine.
groupvec = seq(KMA15GroupsPC)
groupnames = KMA15GroupsPC
maindir = "BAYES/2014-2015 Mixtures 46loci/Output"
mixvec = "SALITC14_2_Middle"
mix = "SALITC14_2_Middle"
prior = ""
ext = "RGN"
nchains = 5
burn = 0.5
alpha = 0.1
PosteriorOutput = TRUE



require(coda)

G <- max(groupvec)

C <- length(groupvec)

nummix <- length(mixvec)

results <- setNames(vector("list",nummix),mixvec)

Output <- setNames(vector("list",nummix),mixvec)

filenames <- paste(maindir,"\\",mix,"\\",mix,prior,"Chain",1:nchains,ext,".",ext,sep="")

files <- lapply(filenames,function(filename){mcmc(as.matrix(read.table(filename)[,-1]))})

end <- sapply(files,nrow)

if(length(unique(end))>1){stop("Chains must be the same length!!!")}

end <- end[1]

begin <- floor(burn*end)+1

files4GR <- vector("list",nchains)

for(chain in seq(nchains)){
  
  files4GR[[chain]] <- as.mcmc(t(rowsum(t(files[[chain]][begin:end,]),group=groupvec)))
  
}#chain

files4GR <- as.mcmc.list(files4GR)    
str(files4GR)

# Gelman-Ruben (between chain)
GR <- gelman.diag(files4GR,multivariate=FALSE,transform=TRUE)
GR
gelman.plot(files4GR, ylim = c(1, 2))

# Raftery-Lewis (within chain)
RL <- raftery.diag(files4GR, q = 0.025, r = 0.005, s = 0.95, converge.eps = 0.001) 
RL
RL.BAYES <- raftery.diag(files4GR, q = 0.975, r = 0.02, s = 0.95, converge.eps = 0.001) 
RL.BAYES

RL.BAYES.Median <- raftery.diag(files4GR, q = 0.5, r = 0.02, s = 0.95, converge.eps = 0.001) 
RL.BAYES.Median

summary(files4GR)

# Plot trace and density
plot(files4GR)

# Plot running mean
par(mfrow = c(5, 1), mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 0))
lapply(files4GR, function(chain) {
  plot(sapply(seq(dim(chain)[1]), function(iter) {mean(chain[seq(from = 1, to = iter), 5])}), type = "l", ylim = c(0.45, 0.6))
} )
mtext(text = paste("Alitak 2014 Middle Strata: Frazer", sep = ''), side = 3, outer = TRUE, cex = 1.5)
mtext(text = "Iteration", side = 1, outer = TRUE, cex = 1.5, line = 1.5)
mtext(text = "Running Mean Posterior", side = 2, outer = TRUE, cex = 1.5, line = 1.5)

# Ultimately, we decided to run this out much farther and see if it stabilizes on one of the modes (200K vs. 40K)
# Also going to re-do the proof test such that we grab only early-Karluk/early-Ayakulik vs. Frazer and then again with late-Karluk/late-Frazer


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 1 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA473PopsGroupVec31 <- as.numeric(readClipboard())
dput(x = KMA473PopsGroupVec31, file = "Objects/KMA473PopsGroupVec31.txt")
KMA31GroupsPC <- readClipboard()
dput(x = KMA31GroupsPC, file = "Objects/KMA31GroupsPC.txt")

Round1Mixtures_2014_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci/Output/40K Iterations", 
                                                                   mixvec = Round1Mixtures_2014, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round1Mixtures_2014_31RG_Estimates)
QuickBarplot(mixvec = Round1Mixtures_2014[c(3, 4, 2, 1, 5)], estimatesstats = Round1Mixtures_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round1Mixtures_2014_Header[c(3, 4, 2, 1, 5)])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 Repeated MSA files for BAYES Alitak 2014 Middle Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping out a new control file to re-analyze "SALITC14_2_Middle" with a 200K posterior as opposed to our standard 40K
## Perhaps it will converge on one of the modes?

# Just use original mixture, no need to create a new one!
# CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = "SALITC14_2_Middle", dir = "BAYES/2014-2015 Mixtures 46loci/Mixture", type = "BAYES", PT = FALSE)

## Dumping a New Control file that goes out further!

CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = "SALITC14_2_Middle", basename = "KMA473Pops46Markers", suffix = "", nreps = 200000, nchains = 5,
                      groupvec = KMA473PopsGroupVec15, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                      seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")


## Create output directory
dir.create("BAYES/2014-2015 Mixtures 46loci/Output/SALITC14_2_Middle")
# Will rename later, but leave as is for now for BayesCopyPaste.GCL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pulling entire posterior with no burn-in
Round1Alitak_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                            maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
                                                            mixvec = "SALITC14_2_Middle", prior = "",  
                                                            ext = "RGN", nchains = 5, burn = 0, alpha = 0.1, PosteriorOutput = TRUE)

# Plotting entire posterior with no burn-in
nchains <- 5
RG <- "Frazer"
# RG <- "Ayakulik"
par(mfrow = c(5, 1), mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 0))
sapply(seq(nchains), function(chain) {
  posterior.length <- dim(Round1Alitak_2014_Estimates$Output$SALITC14_2_Middle)[1] / nchains
  plot(Round1Alitak_2014_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA15GroupsPC == RG)], 
       type = "l", ylim = c(0, 1), xlab = '', ylab = '')
  abline(v = posterior.length / 2, lwd = 2)
  text(x = posterior.length / 4, y = 1, labels = "Burn-in", pos = 1)
  text(x = posterior.length / 4 * 3, y = 1, labels = "Posterior", pos = 1)} )
mtext(text = paste("Alitak 2014 Middle Strata: ", RG, sep = ''), side = 3, outer = TRUE, cex = 1.5)
mtext(text = "Iteration (5 chains each)", side = 1, outer = TRUE, cex = 1.5, line = 1.5)
mtext(text = "Posterior", side = 2, outer = TRUE, cex = 1.5, line = 1.5)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Looking at different ways to measure MCMC convergence, both within and between chains
# Following methods from <http://www.people.fas.harvard.edu/~plam/teaching/methods/convergence/convergence_print.pdf>
# All done with the "coda" package; code is taken from CustomCombine.
groupvec = seq(KMA15GroupsPC)
groupnames = KMA15GroupsPC
maindir = "BAYES/2014-2015 Mixtures 46loci/Output"
mixvec = "SALITC14_2_Middle"
mix = "SALITC14_2_Middle"
prior = ""
ext = "RGN"
nchains = 5
burn = 0.5
alpha = 0.1
PosteriorOutput = TRUE



require(coda)

G <- max(groupvec)

C <- length(groupvec)

nummix <- length(mixvec)

results <- setNames(vector("list",nummix),mixvec)

Output <- setNames(vector("list",nummix),mixvec)

filenames <- paste(maindir,"\\",mix,"\\",mix,prior,"Chain",1:nchains,ext,".",ext,sep="")

files <- lapply(filenames,function(filename){mcmc(as.matrix(read.table(filename)[,-1]))})

end <- sapply(files,nrow)

if(length(unique(end))>1){stop("Chains must be the same length!!!")}

end <- end[1]

begin <- floor(burn*end)+1

files4GR <- vector("list",nchains)

for(chain in seq(nchains)){
  
  files4GR[[chain]] <- as.mcmc(t(rowsum(t(files[[chain]][begin:end,]),group=groupvec)))
  
}#chain

files4GR <- as.mcmc.list(files4GR)    
str(files4GR)

# Gelman-Ruben (between chain)
GR <- gelman.diag(files4GR,multivariate=FALSE,transform=TRUE)
GR
gelman.plot(files4GR, ylim = c(1, 2))

# Raftery-Lewis (within chain)
RL <- raftery.diag(files4GR, q = 0.025, r = 0.005, s = 0.95, converge.eps = 0.001) 
RL
RL.BAYES <- raftery.diag(files4GR, q = 0.975, r = 0.02, s = 0.95, converge.eps = 0.001) 
RL.BAYES

RL.BAYES.Median <- raftery.diag(files4GR, q = 0.5, r = 0.02, s = 0.95, converge.eps = 0.001) 
RL.BAYES.Median

summary(files4GR)

# Plot trace and density
plot(files4GR)

# Plot running mean
par(mfrow = c(5, 1), mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 0))
lapply(files4GR, function(chain) {
  plot(sapply(seq(dim(chain)[1]), function(iter) {mean(chain[seq(from = 1, to = iter), 5])}), type = "l", ylim = c(0.45, 0.6))
} )
mtext(text = paste("Alitak 2014 Middle Strata: Frazer", sep = ''), side = 3, outer = TRUE, cex = 1.5)
mtext(text = "Iteration", side = 1, outer = TRUE, cex = 1.5, line = 1.5)
mtext(text = "Running Mean Posterior", side = 2, outer = TRUE, cex = 1.5, line = 1.5)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Re-doing some baselin proof tests with only early-Ayakulik and only early-Karluk ####
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
# New GroupVec with Ayakulik and Karluk split into Early/Late
KMA473PopsGroupVec17 <- as.numeric(readClipboard())
dput(x = KMA473PopsGroupVec17, file = "Objects/KMA473PopsGroupVec17.txt")
# New Groups
KMA17GroupsPC <- c(PCGroups15[1:5], "Ayakulik Early", "Ayakulik Late", "Karluk Early", "Karluk Late", PCGroups15[8:15])
dput(x = KMA17GroupsPC, file = "Objects/KMA17GroupsPC.txt")
# New Regionally Flat Prior
KMA473Pops17FlatPrior <- Prior.GCL(groupvec = KMA473PopsGroupVec17, groupweights = rep(1 / length(KMA17GroupsPC), length(KMA17GroupsPC)), minval = 0.01)
dput(x = KMA473Pops17FlatPrior, file = "Objects/KMA473Pops17FlatPrior.txt")

KMA473PopsInits
CommonNames473
KMA473Pops

# Population sample sizes
KMA473Pops.SampleSize <- sapply(paste(KMA473Pops.named, ".gcl", sep = ''), function(x) get(x)$n)

# RG samples sizes
KMA473Pops.17RG.SampleSize <- setNames(object = sapply(seq(KMA17GroupsPC), function(RG) {sum(KMA473Pops.SampleSize[which(KMA473PopsGroupVec17 == RG)])} ), nm = KMA17GroupsPC)

## Mixture proof test
# Get fishery scenarios from managers
MixtureProofTestProportions17RG <- read.table(file = "FisheryProofTestScenariosRG17.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = KMA17GroupsPC)
dput(x = t(data.matrix(MixtureProofTestProportions17RG[, -1])), file = "Objects/MixtureProofTestProportions17RG.txt")
MixtureProofTestProportions17RG <- dget(file = "Objects/MixtureProofTestProportions17RG.txt")


FisheryProofTestScenarioNames17RG <- rownames(MixtureProofTestProportions17RG)
dput(x = FisheryProofTestScenarioNames17RG, file = "Objects/FisheryProofTestScenarioNames17RG.txt")
FisheryProofTestScenarioNames17RG <- dget(file = "Objects/FisheryProofTestScenarioNames17RG.txt")


# Sample sizes are adjusted to get whole numbers
MixtureProofTest17RG.SampleSize <- t(apply(MixtureProofTestProportions17RG, 1, function(scenario) {floor(scenario * 400)} ))
apply(MixtureProofTest17RG.SampleSize, 1, sum)
dput(x = MixtureProofTest17RG.SampleSize, file = "Objects/MixtureProofTest17RG.SampleSize.txt")


# Create directories
dir.create("BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit")
invisible(sapply(c("baseline", "control", "mixture", "output"), function(fldr) {dir.create(path = paste("BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.", fldr, sep = ""))} ))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# # ReProofTest????
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ### Fishery Scenario Tests
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KMA473PopsGroups15RepeatedMixProofTests
# 
# #~~~~~~~~~~~~~~~~~~
# ## Get Proof Objects
# for(proof in KMA473PopsGroups15RepeatedMixProofTests) {
#   assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("MixtureProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
# }; rm(proof)
# 
# #~~~~~~~~~~~~~~~~~~
# ## Create BAYES files
# for (proof in KMA473PopsGroups15RepeatedMixProofTests){
#   ReProofTestKS.GCL(sillyvec = KMA473Pops, loci = loci46, groupnames = Groups15Short, groupvec = KMA473PopsGroupVec15, 
#                     ProofTestIDs = get(paste(proof, "Proof", sep = "")), prefix = proof, dir = "BAYES/Mixture Proof Tests/loci46", 
#                     suffix = "", nreps = 40000, nchains = 1, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits[, 1, drop = FALSE], 
#                     type = "BAYES", thin = c(1, 1, 100), switches = "F T F T T T F")
# }
# 
# # Create BAYES Output folders
# dir.create("BAYES/Mixture Proof Tests/loci46/BAYES.output")
# invisible(sapply(KMA473PopsGroups15RepeatedMixProofTests, function(proof) {dir <- paste(getwd(), "/BAYES/Mixture Proof Tests/loci46/BAYES.output/", proof, sep = "")
# dir.create(dir) } ))
# # Decided agains this approach due to complexity and the fact that we ran multiple repeats, but only 1 chain.
# # We really want to comfirm that we get both within and between chains.

# Create mixture proof tests
sapply(FisheryProofTestScenarioNames17RG, function(Scenario) {sapply(1:5, function(Rpt) {
  assign(x = paste(Scenario, Rpt, "Proof", sep = ''), 
         value = ProofTest.GCL(sillyvec = KMA473Pops, loci = loci46, groupnames = KMA17GroupsPC, groupvec = KMA473PopsGroupVec17,
                               samplesize = MixtureProofTest17RG.SampleSize[Scenario, ], prefix = paste(Scenario, Rpt, sep = ''), 
                               dir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit", prprtnl = TRUE, type = "BAYES",
                               suffix = '', nreps = 200000, nchains = 5, priorvec = KMA473Pops17FlatPrior, initmat = KMA473PopsInits, 
                               thin = c(1, 1, 100), switches = "F T F T T T F"), pos = 1)
} )} ); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dput proof test objects
objects(pattern = "Proof$")
KMA473PopsGroups17RepeatedMixProofTests <- paste(rep(FisheryProofTestScenarioNames17RG, each = 5), 1:5, sep = '')
dput(x = KMA473PopsGroups17RepeatedMixProofTests, file = "Objects/KMA473PopsGroups17RepeatedMixProofTests.txt")

dir.create("MixtureProofTests objects/RG17")
invisible(sapply(KMA473PopsGroups17RepeatedMixProofTests, function(proof) {dput(x = get(paste(proof, "Proof", sep = "")), file = paste("MixtureProofTests objects/RG17/", proof, "Proof.txt", sep = ""))} ))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create BAYES.output folders
invisible(sapply(KMA473PopsGroups17RepeatedMixProofTests, function(proof) {dir <- paste(getwd(), "/BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.output/", proof, sep = "")
dir.create(dir) }))

rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
# This sources all of the new GCL functions to this workspace
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("H:/R Source Scripts/Functions.GCL_KS.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 MSA files for BAYES 2014 Middle Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Indetify Strata to Run
KMA2014Strata_2_Middle <- grep(pattern = "2_Middle", x = KMA2014Strata, value = TRUE)
Round2Mixtures_2014 <- KMA2014Strata_2_Middle[-which(KMA2014Strata_2_Middle == "SALITC14_2_Middle")]  # Remove SALITC14_2_Middle as I want to run that out further
dput(x = Round2Mixtures_2014, file = "Objects/Round2Mixtures_2014.txt")

# Create rolling prior based on Round 1 estimates
Round2Mixtures_2014_Prior <- sapply(Round1Mixtures_2014_EstimatesStats[1:4], function(Mix) {Prior.GCL(groupvec = KMA473PopsGroupVec15, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round2Mixtures_2014_Prior) <- gsub(pattern = "1_Early", replacement = "2_Middle", x = names(Round2Mixtures_2014_Prior))  # This changes the names
dput(x = Round2Mixtures_2014_Prior, file = "Objects/Round2Mixtures_2014_Prior.txt")
str(Round2Mixtures_2014_Prior)


## Dumping Mixture files
sapply(Round2Mixtures_2014, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round2Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round2Mixtures_2014_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2014, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                              maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
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
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
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
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round2Mixtures_2014_Header, set.mfrow = c(5, 3), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round2Mixtures_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round2Mixtures_2014_Header)

## Quick barplot
QuickBarplot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round2Mixtures_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round2Mixtures_2014_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimates = Round2Mixtures_2014_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round2Mixtures_2014_Header)
rm(Round2Mixtures_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 Repeated MSA files for BAYES 80K Iterations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping out a new control file to re-analyze all Round 2 mixtures with a 80K posterior as opposed to our standard 40K
## Perhaps it will converge on one of the modes?

# Just use original mixture, no need to create a new one!

## Dumping a New Control file that goes out further!

sapply(Round2Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round2Mixtures_2014_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2014, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})
# Will rename later, but leave as is for now for BayesCopyPaste.GCL


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures80K_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                                 maindir = "BAYES/2014-2015 Mixtures 46loci/Output/80K Iterations", 
                                                                 mixvec = Round2Mixtures_2014, prior = "",  
                                                                 ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures80K_2014_Estimates, file = "Estimates objects/Round2Mixtures80K_2014_Estimates.txt")
dput(Round2Mixtures80K_2014_Estimates$Stats, file = "Estimates objects/Round2Mixtures80K_2014_EstimatesStats.txt")

Round2Mixtures80K_2014_Estimates <- dget(file = "Estimates objects/Round2Mixtures80K_2014_Estimates.txt")
Round2Mixtures80K_2014_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures80K_2014_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures80K_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures80K_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SAYAKC14_2_Middle and SUYAKC14_2_
sapply(Round2Mixtures_2014[c(3, 4, 2, 1)], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round2Mixtures80K_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures80K_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures80K_2014_Estimates$Output)
Round2Mixtures_2014_Header <- setNames(object = c("Ayakulik Section Middle June 28-July 25, 2014",
                                                  "Karluk Section Middle June 28-July 25, 2014",
                                                  "Uganik Section Middle June 28-July 25, 2014",
                                                  "Uyak Section Middle June 28-July 25, 2014"), 
                                       nm = Round2Mixtures_2014)
dput(x = Round2Mixtures_2014_Header, file = "Objects/Round2Mixtures_2014_Header.txt")

PlotPosterior(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], output = Round2Mixtures80K_2014_Estimates$Output, 
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round2Mixtures_2014_Header, set.mfrow = c(5, 3), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round2Mixtures80K_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round2Mixtures_2014_Header)

## Quick barplot
QuickBarplot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round2Mixtures80K_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round2Mixtures_2014_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimates = Round2Mixtures80K_2014_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round2Mixtures_2014_Header)
rm(Round2Mixtures80K_2014_Estimates)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
sapply(Round2Mixtures_2014[c(3, 4, 2, 1)], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round2Mixtures_2014_EstimatesStats[[Mix]][, "median"],
                                                 Round2Mixtures80K_2014_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round2Mixtures_2014_EstimatesStats[[Mix]][, "5%"],
                                               Round2Mixtures80K_2014_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round2Mixtures_2014_EstimatesStats[[Mix]][, "95%"],
                                               Round2Mixtures80K_2014_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round2Mixtures_2014_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA15GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from June 40K to July 80K
Mixtures_2014_Header <- setNames(object = c("Ayakulik Section 2014",
                                            "Karluk Section 2014",
                                            "Uganik Section 2014",
                                            "Uyak Section 2014"), 
                                 nm = KMA2014[2:5])

sapply(KMA2014[c(4, 5, 3, 2)], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round1Mixtures_2014_EstimatesStats[[paste(Mix, "_1_Early", sep = '')]][, "median"],
                                                 Round2Mixtures80K_2014_EstimatesStats[[paste(Mix, "_2_Middle", sep = '')]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA15GroupsPC, c("Early", "Middle")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round1Mixtures_2014_EstimatesStats[[paste(Mix, "_1_Early", sep = '')]][, "5%"],
                                               Round2Mixtures80K_2014_EstimatesStats[[paste(Mix, "_2_Middle", sep = '')]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("Early", "Middle")))),
                      ci.u = t(matrix(data = c(Round1Mixtures_2014_EstimatesStats[[paste(Mix, "_1_Early", sep = '')]][, "95%"],
                                               Round2Mixtures80K_2014_EstimatesStats[[paste(Mix, "_2_Middle", sep = '')]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("Early", "Middle")))),
                      ylim = c(0, 100), col = colorpanel(n = 2, low = "blue", high = "white"), yaxt = 'n', xaxt = 'n', 
                      main = Mixtures_2014_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA15GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("Early - June 1-27", "Middle - June 28 - July 25"), fill = colorpanel(n = 2, low = "blue", high = "white"), bty = 'n', cex = 1.5)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 2 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA473PopsGroupVec31 <- as.numeric(readClipboard())
KMA31GroupsPC <- readClipboard()

Round2Mixtures_2014_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci/Output/80K Iterations", 
                                                                   mixvec = Round2Mixtures_2014, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round2Mixtures_2014_31RG_Estimates)
QuickBarplot(mixvec = Round2Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round2Mixtures_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round2Mixtures_2014_Header[c(3, 4, 2, 1)])
QuickBarplot(mixvec = Round1Mixtures_2014[5], estimatesstats = Round1Mixtures_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round1Mixtures_2014_Header[5])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 2 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SKARLC14_2_Middle and SUGANC14_2_Middle as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from SAYAKC14_2_Middle and SUYAKC14_2_Middle

sapply(Round2Mixtures_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SAYAKC14_2_Middle and SUYAKC14_2_
sapply(Round2Mixtures80K_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SAYAKC14_2_Middle and SUYAKC14_2_

str(Round2Mixtures_2014_Estimates)

Round2Mixtures_2014_Estimates_Final <- list(Stats = c(Round2Mixtures80K_2014_Estimates$Stats[1],
                                                         Round2Mixtures_2014_Estimates$Stats[2:3],
                                                         Round2Mixtures80K_2014_Estimates$Stats[4]),
                                            Output = c(Round2Mixtures80K_2014_Estimates$Output[1],
                                                       Round2Mixtures_2014_Estimates$Output[2:3],
                                                       Round2Mixtures80K_2014_Estimates$Output[4]))
str(Round2Mixtures_2014_Estimates_Final)
dput(x = Round2Mixtures_2014_Estimates_Final, file = "Estimates objects/Round2Mixtures_2014_Estimates_Final.txt")
dput(x = Round2Mixtures_2014_Estimates_Final$Stats, file = "Estimates objects/Round2Mixtures_2014_EstimatesStats_Final.txt")

Round2Mixtures_2014_Estimates_Final <- dget(file = "Estimates objects/Round2Mixtures_2014_Estimates_Final.txt")
Round2Mixtures_2014_EstimatesStats_Final <- dget(file = "Estimates objects/Round2Mixtures_2014_EstimatesStats_Final.txt")

sapply(Round2Mixtures_2014_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SAYAKC14_2_Middle and SUYAKC14_2_


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 MSA files for BAYES 2014 Middle Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Indetify Strata to Run
KMA2014Strata_3_Late <- grep(pattern = "3_Late", x = KMA2014Strata, value = TRUE)
Round3Mixtures_2014 <- KMA2014Strata_3_Late[-which(KMA2014Strata_3_Late == "SALITC14_3_Late")]  # Remove SALITC14_3_Late as I want to run that out further
dput(x = Round3Mixtures_2014, file = "Objects/Round3Mixtures_2014.txt")

# Create rolling prior based on Round 2 estimates
Round3Mixtures_2014_Prior <- sapply(Round2Mixtures_2014_EstimatesStats_Final[1:4], function(Mix) {Prior.GCL(groupvec = KMA473PopsGroupVec15, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round3Mixtures_2014_Prior) <- gsub(pattern = "2_Middle", replacement = "3_Late", x = names(Round3Mixtures_2014_Prior))  # This changes the names
dput(x = Round3Mixtures_2014_Prior, file = "Objects/Round3Mixtures_2014_Prior.txt")
str(Round3Mixtures_2014_Prior)


## Dumping Mixture files
sapply(Round3Mixtures_2014, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round3Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round3Mixtures_2014_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round3Mixtures_2014, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 3 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2014_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
  mixvec = Round3Mixtures_2014, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round3Mixtures_2014_Estimates, file = "Estimates objects/Round3Mixtures_2014_Estimates.txt")
dput(Round3Mixtures_2014_Estimates$Stats, file = "Estimates objects/Round3Mixtures_2014_EstimatesStats.txt")

Round3Mixtures_2014_Estimates <- dget(file = "Estimates objects/Round3Mixtures_2014_Estimates.txt")
Round3Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures_2014_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round3Mixtures_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round3Mixtures_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC14_3_Late
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round3Mixtures_2014[c(3, 4, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round3Mixtures_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round3Mixtures_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round3Mixtures_2014_Estimates$Output)
Round3Mixtures_2014_Header <- setNames(object = c("Ayakulik Section Middle July 26-August 29, 2014",
                                                  "Karluk Section Middle July 26-August 29, 2014",
                                                  "Uganik Section Middle July 26-August 29, 2014",
                                                  "Uyak Section Middle July 26-August 29, 2014"), 
                                       nm = Round3Mixtures_2014)
dput(x = Round3Mixtures_2014_Header, file = "Objects/Round3Mixtures_2014_Header.txt")

PlotPosterior(mixvec = Round3Mixtures_2014[c(3, 4, 2, 1)], output = Round3Mixtures_2014_Estimates$Output, 
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round3Mixtures_2014_Header, set.mfrow = c(5, 3), thin = 10)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 3 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round3Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round3Mixtures_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round3Mixtures_2014_Header)

## Quick barplot
QuickBarplot(mixvec = Round3Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round3Mixtures_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round3Mixtures_2014_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round3Mixtures_2014[c(3, 4, 2, 1)], estimates = Round3Mixtures_2014_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round3Mixtures_2014_Header)
rm(Round3Mixtures_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 Repeated MSA files for BAYES 80K Iterations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping out a new control file to re-analyze SUGANC14_3_Late Round 3 mixtures with a 80K posterior as opposed to our standard 40K
## Perhaps it will converge on one of the modes for South of Cape Suckling?

# Just use original mixture, no need to create a new one!

## Dumping a New Control file that goes out further!

sapply(Round3Mixtures_2014[3], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round3Mixtures_2014_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round3Mixtures_2014[3], function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})
# Will rename later, but leave as is for now for BayesCopyPaste.GCL


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot results KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")
Round2Mixtures_2014_EstimatesStats_Final <- dget("Estimates objects/Round2Mixtures_2014_EstimatesStats_Final.txt")
Round3Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round3Mixtures_2014_EstimatesStats.txt")

KMA2014Strata_EstimatesStats <- c(Round1Mixtures_2014_EstimatesStats, Round2Mixtures_2014_EstimatesStats_Final, Round3Mixtures_2014_EstimatesStats)
str(KMA2014Strata_EstimatesStats)


# plot.ci.extract <- function(x, mixture, group, stat.col) {
#   list.index <- grep(pattern = mixture, x = names(x))  
#   sapply(list.index, function(i) {x[[i]][group, stat.col]} )
# }


# emf(file = "V:/Presentations/Regional/4_Westward/Sockeye/COMFISH/KarlukScenarioFigure.emf", width = 7, height = 7, family = "sans", bg = "white")

# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             1, 3,
                             1, 4,
                             5, 6), nrow = 4, ncol = 2, byrow = TRUE)

TempMix <- c("_1_Early", "_2_Middle", "_3_Late")
TempLegend14 <- c("June 1-27", "June 28-July 25", "July 26-August 29")
TempLegend15 <- c("June 1-July 3", "July 4-August 1", "August 2-August 29")
TempLegend16 <- c("June 1-27", "June 28-July 25", "July 26-August 29")
GeoMix <- c("SALITC", "SAYAKC", "SKARLC", "SUGANC", "SUYAKC")
GeoHeader <- setNames(object = c(paste("Cape Alitak/Humpy Deadman 257-10,20,50,60,70", sep = ''),
                                 paste("Ayakulik/Halibut Bay 256-10", "\u2013", "256-30", sep = ''),
                                 paste("Karluk/Sturgeon 255-10", "\u2013", "255-20; 256-40", sep = ''),
                                 paste("Uganik/Kupreanof 253", sep = ''),
                                 paste("Uyak Bay 254", sep = '')),
                      nm = GeoMix)
darkcolor <- "blue"
Estimates <- KMA2014Strata_EstimatesStats
Groups <- KMA15GroupsPC
Groups2Rows <- KMA15GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.5
cex.xaxis <- 0.7
cex.main <- 2
cex.leg <- 1.5
ci.lwd <- 2.5


# layout(mat = layoutmat, widths = c(0.1, 1, 1), heights = c(0.9, 0.9, 1, 0.1))
# par(mar = rep(0, 4))

# Y-axis label
# plot.new()
# text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2014
# par(mar = c(1, 1, 1, 1))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Temp - Comment this section out for 3 panel plot!!!
par(mar = c(3.6, 5.1, 3.1, 1.1))
sapply(GeoMix[2:5], function(geomix) {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Barplot1 <- barplot2(height = t(sapply(TempMix, function(tempmix) {Estimates[[paste(geomix, "14", tempmix, sep = '')]][, "median"]})) * 100, 
                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                       ci.l = t(sapply(TempMix, function(tempmix) {Estimates[[paste(geomix, "14", tempmix, sep = '')]][, "5%"]})) * 100, 
                       ci.u = t(sapply(TempMix, function(tempmix) {Estimates[[paste(geomix, "14", tempmix, sep = '')]][, "95%"]})) * 100, 
                       ylim = c(0, 100), col = colorpanel(length(TempMix), low = darkcolor, high = "white"), cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14, x = "topleft", fill = colorpanel(length(TempMix), low = darkcolor, high = "white"), border = "black", bty = "n", cex = cex.leg, title="2014")
  abline(h = 0)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Temp - Comment this section out for 3 panel plot!!! 
  mtext(text = "Percentage of Catch", side = 2, cex = cex.yaxis, line = 3)
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot1, 2, mean), adj = 0.5, cex = cex.xaxis)
  mtext(text = GeoHeader[geomix], side = 3, cex = cex.main)
})
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015
# par(mar = c(1, 1, 1, 1))
# Barplot2 <- barplot2(height = t(sapply(TempMix, function(tempmix) {Estimates[[paste(GeoMix, "15", tempmix, sep = '')]][, "median"]})) * 100, 
#                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
#                      ci.l = t(sapply(TempMix, function(tempmix) {Estimates[[paste(GeoMix, "15", tempmix, sep = '')]][, "5%"]})) * 100, 
#                      ci.u = t(sapply(TempMix, function(tempmix) {Estimates[[paste(GeoMix, "15", tempmix, sep = '')]][, "95%"]})) * 100, 
#                      ylim = c(0, 100), col = colorpanel(length(TempMix), low = "blue", high = "white"), cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
# axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
# legend(legend = TempLegend15, x = "topleft", fill = colorpanel(length(TempMix), low = "blue", high = "white"), border = "black", bty = "n", cex = 1, title="2015")
# abline(h = 0)
plot.new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2016
# par(mar = c(2, 1, 1, 1))
# Barplot3 <- barplot2(height = t(sapply(TempMix, function(tempmix) {Estimates[[paste(GeoMix, "16", tempmix, sep = '')]][, "median"]})) * 100, 
#                      beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
#                      ci.l = t(sapply(TempMix, function(tempmix) {Estimates[[paste(GeoMix, "16", tempmix, sep = '')]][, "5%"]})) * 100, 
#                      ci.u = t(sapply(TempMix, function(tempmix) {Estimates[[paste(GeoMix, "16", tempmix, sep = '')]][, "95%"]})) * 100, 
#                      ylim = c(0, 100), col = colorpanel(length(TempMix), low = "blue", high = "white"), cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
# axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
# legend(legend = TempLegend16, x = "topleft", fill = colorpanel(length(TempMix), low = "blue", high = "white"), border = "black", bty = "n", cex = 1, title="2016")
# abline(h = 0)
mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot3, 2, mean), adj = 0.5, cex = cex.xaxis)
plot.new()

## Blank Corner
par(mar = rep(0, 4))
plot.new()

## x-axis label
par(mar = rep(0, 4))
plot.new()
text(x = 0.5, y = 0.5, labels = "Reporting Group", cex = cex.lab)


dev.off()











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