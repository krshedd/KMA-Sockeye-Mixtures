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
  while(!require(gplots, quietly = TRUE)){install.packages("gplots")}
  
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
Round1Mixtures200K_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                                  maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
                                                                  mixvec = "SALITC14_2_Middle", prior = "",  
                                                                  ext = "RGN", nchains = 5, burn = 0, alpha = 0.1, PosteriorOutput = TRUE)


# Verify that Gelman-Rubin < 1.2
sapply(Round1Mixtures200K_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round1Mixtures200K_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SAYAKC14_2_Middle and SUYAKC14_2_
sapply(Round1Mixtures_2014[5], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round1Mixtures200K_2014_Estimates$Stats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round1Mixtures200K_2014_Estimates$Stats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round1Mixtures200K_2014_Estimates$Output)

PlotPosterior(mixvec = Round1Mixtures_2014[5], output = Round1Mixtures200K_2014_Estimates$Output, 
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round1Mixtures_2014_Header[5], set.mfrow = c(5, 3), thin = 100)

# Plotting entire posterior with no burn-in
nchains <- 5
# RG <- "Frazer"
RG <- "Ayakulik"
par(mfrow = c(5, 1), mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 0))
sapply(seq(nchains), function(chain) {
  posterior.length <- dim(Round1Mixtures200K_2014_Estimates$Output$SALITC14_2_Middle)[1] / nchains
  plot(Round1Mixtures200K_2014_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA15GroupsPC == RG)], 
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
#### Sneak Peak Round 1 Middle Alitak Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures200K_2014_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                       maindir = "BAYES/2014-2015 Mixtures 46loci/Output/200K Iterations", 
                                                                       mixvec = Round1Mixtures_2014[5], prior = "",  
                                                                       ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round1Mixtures200K_2014_31RG_Estimates)
QuickBarplot(mixvec = Round1Mixtures_2014[5], estimatesstats = Round1Mixtures200K_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round1Mixtures_2014_Header[5])

ViolinPlot(mixvec = Round1Mixtures_2014[5], estimates = Round1Mixtures200K_2014_31RG_Estimates, groups = KMA31GroupsPC, colors = KMA15Colors[c(1,1,1,1,2,3,4,4,4,5,6,6,7,7,8,9,10,11,12,13,13,13,13,13,13,13,13,13,14,14,15)], header = Round1Mixtures_2014_Header[5])

# Still won't converge with 200K iterations/chain
# Appears to be bouncing between Frazer and Ayakulik Late

Round1Mixtures200KNoBurn_2014_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                             maindir = "BAYES/2014-2015 Mixtures 46loci/Output/200K Iterations", 
                                                                             mixvec = Round1Mixtures_2014[5], prior = "",  
                                                                             ext = "BOT", nchains = 5, burn = 0, alpha = 0.1, PosteriorOutput = TRUE)

# Plotting trace entire posterior with no burn-in
nchains <- 5
par(mfrow = c(5, 1), mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 0))
sapply(seq(nchains), function(chain) {
  posterior.length <- dim(Round1Mixtures200KNoBurn_2014_31RG_Estimates$Output$SALITC14_2_Middle)[1] / nchains
  plot(Round1Mixtures200KNoBurn_2014_31RG_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA31GroupsPC == "Frazer")], 
       type = "l", ylim = c(0, 1), xlab = '', ylab = '', col = "cyan")
  points(Round1Mixtures200KNoBurn_2014_31RG_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA31GroupsPC == "Ayakulik Late")],
         type = "l", add = TRUE, col = "blue")
  points(Round1Mixtures200KNoBurn_2014_31RG_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA31GroupsPC == "Ayakulik Early")],
         type = "l", add = TRUE, col = "black")
  abline(v = posterior.length / 2, lwd = 2)
  text(x = posterior.length / 4, y = 1, labels = "Burn-in", pos = 1)
  text(x = posterior.length / 4 * 3, y = 1, labels = "Posterior", pos = 1)} )
mtext(text = "Alitak 2014 Middle Strata", side = 3, outer = TRUE, cex = 1.5)
mtext(text = "Iteration (5 chains each)", side = 1, outer = TRUE, cex = 1.5, line = 1.5)
mtext(text = "Posterior", side = 2, outer = TRUE, cex = 1.5, line = 1.5)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))


# Plotting histogram entire posterior with no burn-in
nchains <- 5
par(mfrow = c(5, 1), mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 0))
sapply(seq(nchains), function(chain) {
  posterior.length <- dim(Round1Mixtures200KNoBurn_2014_31RG_Estimates$Output$SALITC14_2_Middle)[1] / nchains
  hist(Round1Mixtures200KNoBurn_2014_31RG_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA31GroupsPC == "Frazer")], 
       xlim = c(0, 1), xlab = '', ylab = '', col = "cyan", breaks = seq(0 , 1, 0.01), main = '')
  hist(Round1Mixtures200KNoBurn_2014_31RG_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA31GroupsPC == "Ayakulik Late")],
         add = TRUE, col = "blue", breaks = seq(0 , 1, 0.01), xlab = '', main = '')
  # hist(Round1Mixtures200KNoBurn_2014_31RG_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA31GroupsPC == "Ayakulik Early")],
  #        add = TRUE, col = "black", breaks = seq(0 , 1, 0.01), xlab = '')
  if(chain == 1) {legend("topright", legend = c("Frazer", "Ayakulik Late"), fill = c("cyan", "blue"), bty = 'n', cex = 2)}
})
mtext(text = "Alitak 2014 Middle Strata", side = 3, outer = TRUE, cex = 1.5)
mtext(text = "Posterior", side = 1, outer = TRUE, cex = 1.5, line = 1.5)
mtext(text = "Iteration (5 chains each)", side = 2, outer = TRUE, cex = 1.5, line = 1.5)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))


# Plotting histogram original 40K, burn first half
nchains <- 5
par(mfrow = c(5, 1), mar = c(2.1, 2.1, 1.1, 1.1), oma = c(4.1, 4.1, 3.1, 0))
sapply(seq(nchains), function(chain) {
  posterior.length <- dim(Round1Mixtures_2014_Estimates$Output$SALITC14_2_Middle)[1] / nchains
  hist(Round1Mixtures_2014_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA15GroupsPC == "Frazer")], 
       xlim = c(0, 1), xlab = '', ylab = '', col = "cyan", breaks = seq(0 , 1, 0.01), main = '')
  hist(Round1Mixtures_2014_Estimates$Output$SALITC14_2_Middle[seq(from = (chain - 1) * posterior.length + 1, length.out = posterior.length), which(KMA15GroupsPC == "Ayakulik")],
       add = TRUE, col = "blue", breaks = seq(0 , 1, 0.01), xlab = '', main = '')
  if(chain == 1) {legend("topright", legend = c("Frazer", "Ayakulik"), fill = c("cyan", "blue"), bty = 'n', cex = 2)}
})
mtext(text = "Alitak 2014 Middle Strata", side = 3, outer = TRUE, cex = 1.5)
mtext(text = "Posterior", side = 1, outer = TRUE, cex = 1.5, line = 1.5)
mtext(text = "Iteration (5 chains each)", side = 2, outer = TRUE, cex = 1.5, line = 1.5)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))

# These fish are clearly bouncing back and forth between Frazer and Ayakulik Late


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 Middle Alitak Final Decisions ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round1Mixtures200K_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                                  maindir = "BAYES/2014-2015 Mixtures 46loci/Output/200K Iterations/", 
                                                                  mixvec = "SALITC14_2_Middle", prior = "",  
                                                                  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

str(Round1Mixtures200K_2014_Estimates)

Round1Mixtures200K_2014_EstimatesStats <- Round1Mixtures200K_2014_Estimates$Stats
Round1Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")

sapply(Round1Mixtures_2014[5], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round1Mixtures_2014_EstimatesStats[[Mix]][, "median"],
                                                 Round1Mixtures200K_2014_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round1Mixtures_2014_EstimatesStats[[Mix]][, "5%"],
                                               Round1Mixtures200K_2014_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round1Mixtures_2014_EstimatesStats[[Mix]][, "95%"],
                                               Round1Mixtures200K_2014_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round1Mixtures_2014_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA15GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "200K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})

# No difference in summary statistics (median and 90% CI) between 40K and 200K iterations
# Keep 40K output as results


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Re-doing some baseline proof tests with only early-Ayakulik and only early-Karluk ####
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

# Confirm that groupvec is correct
sapply(KMA17GroupsPC[5:9], function(group) {KMA473Pops[which(KMA473PopsGroupVec17 == which(KMA17GroupsPC == group))]} )

# New Regionally Flat Prior
KMA473Pops17FlatPrior <- Prior.GCL(groupvec = KMA473PopsGroupVec17, groupweights = rep(1 / length(KMA17GroupsPC), length(KMA17GroupsPC)), minval = 0.01)
dput(x = KMA473Pops17FlatPrior, file = "Objects/KMA473Pops17FlatPrior.txt")

KMA473PopsInits
CommonNames473
KMA473Pops

# Population sample sizes
KMA473Pops.SampleSize <- sapply(paste(KMA473Pops.named, ".gcl", sep = ''), function(x) get(x)$n)

# RG samples sizes
KMA473Pops.15RG.SampleSize <- setNames(object = sapply(seq(KMA15GroupsPC), function(RG) {sum(KMA473Pops.SampleSize[which(KMA473PopsGroupVec15 == RG)])} ), nm = KMA15GroupsPC)

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

# Doesn't work, but should...
invisible(sapply(KMA473PopsGroups17RepeatedMixProofTests, function(proof) {
  assign(x = dget(file = paste("MixtureProofTests objects/RG17/", proof, "Proof.txt", sep = "")), 
         value = paste(proof, "Proof", sep = ""), pos = 1)
} ))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create BAYES.output folders
invisible(sapply(KMA473PopsGroups17RepeatedMixProofTests, function(proof) {dir <- paste(getwd(), "/BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.output/", proof, sep = "")
dir.create(dir) }))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Summarize Output


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# early.Ayakulik1

# BOT file
KMA473PopsGroups17RepeatedMixProofTests_early.Ayakulik1_Estimates_BOT_Stats <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_early.Ayakulik1_Estimates_BOT_Stats.txt")


# Summarize pop means
KMA473PopsGroups17RepeatedMixProofTests_early.Ayakulik1_Estimates_BOT_Stats$early.Ayakulik1[KMA473PopsGroupVec17 %in% c(5,6,7), ]

# Determine true pop means
early.Ayakulik1Proof <- dget(file = "MixtureProofTests objects/RG17/early.Ayakulik1Proof.txt")
sapply(c("Frazer", "Ayakulik Early", "Ayakulik Late"), function(RG) {sapply(early.Ayakulik1Proof[[RG]], function(pop) {length(pop) / length(unlist(early.Ayakulik1Proof))} )} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# middle.Alitak.early1

# BOT file
KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.early1_Estimates_BOT_Stats <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA473Pops), groupnames = KMA473Pops, 
  maindir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.output", 
  mixvec = paste(FisheryProofTestScenarioNames17RG[2], 1, sep = ''), prior = "",  
  ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE); beep(5)

dput(x = KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.early1_Estimates_BOT_Stats,
     file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.early1_Estimates_BOT_Stats.txt")

KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.early1_Estimates_BOT_Stats <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.early1_Estimates_BOT_Stats.txt")

str(KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.early1_Estimates_BOT_Stats)


# Summarize pop means
KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.early1_Estimates_BOT_Stats$middle.Alitak.early1[KMA473PopsGroupVec17 %in% c(5,6,7), ]

# Determine true pop means
middle.Alitak.early1Proof <- dget(file = "MixtureProofTests objects/RG17/middle.Alitak.early1Proof.txt")
sapply(c("Frazer", "Ayakulik Early", "Ayakulik Late"), function(RG) {sapply(middle.Alitak.early1Proof[[RG]], function(pop) {length(pop) / length(unlist(middle.Alitak.early1Proof))} )} )

# Both "SAYAK00.SAYAK08L.SAYAK11L.SAYAK12L" "SAYAK12E" are sucking fish!!!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# middle.Alitak.501

# BOT file
KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.501_Estimates_BOT_Stats <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.501_Estimates_BOT_Stats.txt")


# Summarize pop means
KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.501_Estimates_BOT_Stats$middle.Alitak.501[KMA473PopsGroupVec17 %in% c(5,6,7), ]

# Determine true pop means
middle.Alitak.501Proof <- dget(file = "MixtureProofTests objects/RG17/middle.Alitak.501Proof.txt")
sapply(c("Frazer", "Ayakulik Early", "Ayakulik Late"), function(RG) {sapply(middle.Alitak.501Proof[[RG]], function(pop) {length(pop) / length(unlist(middle.Alitak.501Proof))} )} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# middle.Alitak.late1

# BOT file
KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.late1_Estimates_BOT_Stats <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.late1_Estimates_BOT_Stats.txt")


# Summarize pop means
KMA473PopsGroups17RepeatedMixProofTests_middle.Alitak.late1_Estimates_BOT_Stats$middle.Alitak.late1[KMA473PopsGroupVec17 %in% c(5,6,7), ]

# Determine true pop means
middle.Alitak.late1Proof <- dget(file = "MixtureProofTests objects/RG17/middle.Alitak.late1Proof.txt")
sapply(c("Frazer", "Ayakulik Early", "Ayakulik Late"), function(RG) {sapply(middle.Alitak.late1Proof[[RG]], function(pop) {length(pop) / length(unlist(middle.Alitak.late1Proof))} )} )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# No Burnin
KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates_NoBurnin <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA17GroupsPC), groupnames = KMA17GroupsPC, 
  maindir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.output", 
  mixvec = paste(FisheryProofTestScenarioNames17RG, 1, sep = ''), prior = "",  
  ext = "RGN", nchains = 5, burn = 0, alpha = 0.1, PosteriorOutput = TRUE)
str(KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates_NoBurnin)

# Quick look at entire posterior with no burnin
KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Header <- setNames(
  object = c("Early Ayakulik Fishery Scenario", 
             "Middle Alitak Fishery Scenario: All Early Run",
             "Middle Alitak Fishery Scenario: 50/50 Early/Late Run",
             "Middle Alitak Fishery Scenario: All Late Run"), 
  nm = paste(FisheryProofTestScenarioNames17RG, 1, sep = ''))

Colors17 <- c("grey60", "black",  "brown",  "red",    "cyan",   "darkblue", "blue",   "purple", "purple4", 
              "darkorange2", "magenta", "yellow3", "green", "darkgreen", "pink2", "darkgoldenrod4", "turquoise4")
dput(x = Colors17, file = "Objects/Colors17.txt")

PlotPosterior(mixvec = paste(FisheryProofTestScenarioNames17RG, 1, sep = ''), 
              output = KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates_NoBurnin$Output, 
              groups = KMA17GroupsPC, 
              colors = Colors17,
              header = KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Header, set.mfrow = c(6, 3), thin = 100)



# Pull posterior with burnin
KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA17GroupsPC), groupnames = KMA17GroupsPC, 
  maindir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.output", 
  mixvec = paste(FisheryProofTestScenarioNames17RG, 1, sep = ''), prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates)


# Dput 1) estimates stats + posterior output & 2) estimates stats
dir.create("Estimates objects/loci46 KarlukAyakulikSplit")
dput(KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates, file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates.txt")
dput(KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates$Stats, file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeat1_EstimatesStats.txt")

KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates.txt")
KMA473PopsGroups17RepeatedMixProofTests_Repeat1_EstimatesStats <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeat1_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SAYAKC14_2_Middle and SUYAKC14_2_
sapply(paste(FisheryProofTestScenarioNames17RG, 1, sep = ''), function(Mix) {
  BarPlot <- barplot2(KMA473PopsGroups17RepeatedMixProofTests_Repeat1_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(KMA473PopsGroups17RepeatedMixProofTests_Repeat1_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA17GroupsPC, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates$Output)

PlotPosterior(mixvec = paste(FisheryProofTestScenarioNames17RG, 1, sep = ''), 
              output = KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates$Output, 
              groups = KMA17GroupsPC, 
              colors = Colors17, 
              header = KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Header, set.mfrow = c(6, 3), thin = 100)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot results
## Quick and dirty plot
QuickPlot(mixvec = paste(FisheryProofTestScenarioNames17RG, 1, sep = ''), estimatesstats = KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates, groups = KMA17GroupsPC, colors = Colors17, header = KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Header)

## Quick barplot
QuickBarplot(mixvec = paste(FisheryProofTestScenarioNames17RG, 1, sep = ''), estimatesstats = KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates, groups = KMA17GroupsPC, header = KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = paste(FisheryProofTestScenarioNames17RG, 1, sep = ''), estimates = KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates, groups = KMA17GroupsPC, colors = Colors17, header = KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Header)
rm(Round2Mixtures_2014_Estimates)

t(MixtureProofTestProportions17RG)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



require(devEMF)


# Barplots
estimatesstats <- KMA473PopsGroups17RepeatedMixProofTests_Repeat1_EstimatesStats
mixvec <- FisheryProofTestScenarioNames17RG
header <- KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Header
groups2rows <- NULL
groups <- KMA17GroupsPC
proportions <- MixtureProofTestProportions17RG

sapply(mixvec, function(mix) {
  
  emf(file = paste("BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/Figures/Barplot_", mix, ".emf", sep = ''), width = 6.5, height = 7.5, family = "Times")
  
  Mix <- paste(mix, 1, sep = '')
  
  par(mfrow = c(1, 1), mar = c(6.1, 5.1, 4.1, 2.1), oma = rep(0, 4))
  
  Barplot <- barplot2(height = estimatesstats[[Mix]][, "median"] * 100,
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 1,
                      ci.l = estimatesstats[[Mix]][, "5%"] * 100,
                      ci.u = estimatesstats[[Mix]][, "95%"] * 100,
                      ylim = c(0, 100), col = "blue", yaxt = 'n', xaxt = 'n',
                      main = header[Mix], ylab = "Percent of Mixture", cex.lab = 1.5, cex.main = 1.5)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  
  if(is.null(groups2rows)) {
    text(x = Barplot[, 1], y = -1, labels = groups, srt = 90, adj =  1, xpd = TRUE, cex = 0.6)
  } else {
    mtext(text = groups2rows, side = 1, line = 1, at = Barplot[, 1], adj = 0.5, cex = 0.7)
  }
  
  segments(x0 = Barplot - 0.5, x1 = Barplot + 0.5, y0 = proportions[mix, ] * 100, y1 = proportions[mix, ] * 100, lwd = 4, col = "red")
  
  dev.off()
})


# ViolinPlots
estimates <- KMA473PopsGroups17RepeatedMixProofTests_Repeat1_Estimates
colors <- Colors17
thin <- 100
wex <- 1


par(mfrow = c(1, 1), mar = c(6.1, 5.1, 4.1, 2.1), oma = rep(0, 4))
sapply(mixvec, function(mix) {
  
  emf(file = paste("BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/Figures/Violinplot_", mix, ".emf", sep = ''), width = 6.5, height = 7.5, family = "Times")
  
  Mix <- paste(mix, 1, sep = '')
  
  plot(estimates$Stats[[Mix]][, "median"], cex = 3, pch = 16, col = colors, ylab = "Proportion of Mixture", ylim = c(0, 1), xlab = "", axes = FALSE, main = header[[Mix]], cex.main = 1.5, cex.lab = 1.5)
  sapply(seq(groups), function(i) {vioplot2(estimates$Output[[Mix]][seq(from = 1, to = nrow(estimates$Output[[Mix]]), by = thin), i], at = i, horizontal = FALSE, col = colors[i], border = TRUE, drawRect = FALSE, rectCol = colors[i], add = TRUE, wex = wex, lwd = 2)})
  arrows(x0 = seq(groups), y0 = estimates$Stats[[Mix]][, "5%"], x1 = seq(groups), y1 = estimates$Stats[[Mix]][, "95%"], angle = 90, code = 3, length = 0.2, lwd = 2)
  points(estimates$Stats[[Mix]][, "median"], cex = 2, pch = 21, col = "white", bg = colors, lwd = 3)
  axis(side = 2, lwd = 3, cex.axis = 1.5)
  # text(x = (seq(groups)) - 0.35, y = 0, labels = groups, srt = 60, pos = 1, offset = 2.5, xpd = TRUE)
  axis(side = 1, labels = NA, at = seq(groups), pos = 0, lwd = 2, tick = FALSE)
  abline(h = 0, lwd = 3, xpd = FALSE)
  
  
  if(is.null(groups2rows)) {
    text(x = seq(groups), y = -.025, labels = groups, srt = 90, adj =  1, xpd = TRUE, cex = 0.7)
  } else {
    mtext(text = groups2rows, side = 1, line = 1, at = seq(groups), adj = 0.5, cex = 0.7)
  }
  
  
  segments(x0 = seq(groups) - 0.4, x1 = seq(groups) + 0.4, y0 = proportions[mix, ], y1 = proportions[mix, ], lwd = 5, col = "red")
  
  dev.off()
} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# More proof tests, 5 repeats, 1 chain each, 40K

# Create mixture proof tests
sapply(FisheryProofTestScenarioNames17RG, function(Scenario) {
  sapply(1:5, function(Rpt) {
    assign(x = paste(Scenario, Rpt, "40KProof", sep = ''), 
           value = ProofTest.GCL(sillyvec = KMA473Pops, loci = loci46, groupnames = KMA17GroupsPC, groupvec = KMA473PopsGroupVec17,
                                 samplesize = MixtureProofTest17RG.SampleSize[Scenario, ], prefix = paste(Scenario, Rpt, sep = ''), 
                                 dir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit", prprtnl = TRUE, type = "BAYES",
                                 suffix = '', nreps = 40000, nchains = 1, priorvec = KMA473Pops17FlatPrior, initmat = KMA473PopsInits, 
                                 thin = c(1, 1, 100), switches = "F T F T T T F"), pos = 1)
  } )
} ); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dput proof test objects
objects(pattern = "Proof$")
KMA473PopsGroups17RepeatedMixProofTests40K <- paste(rep(FisheryProofTestScenarioNames17RG, each = 5), 1:5, "40K", sep = '')
dput(x = KMA473PopsGroups17RepeatedMixProofTests40K, file = "Objects/KMA473PopsGroups17RepeatedMixProofTests40K.txt")

# dir.create("MixtureProofTests objects/RG17")
invisible(sapply(KMA473PopsGroups17RepeatedMixProofTests40K, function(proof) {dput(x = get(paste(proof, "Proof", sep = "")), file = paste("MixtureProofTests objects/RG17/", proof, "Proof.txt", sep = ""))} ))

# Grab proof test values
invisible(sapply(KMA473PopsGroups17RepeatedMixProofTests, function(proof) {
  assign(x = paste(proof, "40KProof", sep = ""), 
         value = dget(file = paste("MixtureProofTests objects/RG17/", proof, "40KProof.txt", sep = "")), pos = 1)
} ))


# Determine where fish are coming from for proof tests
str(early.Ayakulik540KProof, max.level = 1)

proof <- early.Ayakulik540KProof

sapply(proof[c("Frazer", "Ayakulik Early")], function(RG) {
  unlist(sapply(RG, function(pop) {
    if(!is.null(pop)) {length(pop)}
  } ))
} )



sapply(paste(KMA473PopsGroups17RepeatedMixProofTests, "40KProof", sep = ''), function(proof) {
  sapply(get(proof)[c("Frazer", "Ayakulik Early", "Ayakulik Late")], function(RG) {
    unlist(sapply(RG, function(pop) {
      if(!is.null(pop)) {length(pop)}
    } ))
  } )
}, simplify = FALSE )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create BAYES.output folders
invisible(sapply(KMA473PopsGroups17RepeatedMixProofTests, function(proof) {dir <- paste(getwd(), "/BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.output/", proof, sep = "")
dir.create(dir) }))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pull posterior with burnin
KMA473PopsGroups17RepeatedMixProofTests_Repeats5_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA17GroupsPC), groupnames = KMA17GroupsPC, 
  maindir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.output", 
  mixvec = KMA473PopsGroups17RepeatedMixProofTests, prior = "",  
  ext = "RGN", nchains = 1, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(KMA473PopsGroups17RepeatedMixProofTests_Repeats5_Estimates)


# Dput 1) estimates stats + posterior output & 2) estimates stats
# dir.create("Estimates objects/loci46 KarlukAyakulikSplit")
dput(KMA473PopsGroups17RepeatedMixProofTests_Repeats5_Estimates, file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5_Estimates.txt")
dput(KMA473PopsGroups17RepeatedMixProofTests_Repeats5_Estimates$Stats, file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5_EstimatesStats.txt")

KMA473PopsGroups17RepeatedMixProofTests_Repeats5_Estimates <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5_Estimates.txt")
KMA473PopsGroups17RepeatedMixProofTests_Repeats5_EstimatesStats <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2 ## No GR for 1 chain....

# Quick look at raw posterior output
str(KMA473PopsGroups17RepeatedMixProofTests_Repeats5_Estimates$Output)

PlotPosterior(mixvec = KMA473PopsGroups17RepeatedMixProofTests, 
              output = KMA473PopsGroups17RepeatedMixProofTests_Repeats5_Estimates$Output, 
              groups = KMA17GroupsPC, 
              colors = Colors17, 
              header = NULL, set.mfrow = c(6, 3), thin = 1, chains = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SEDM style proof plots



#~~~~~~~~~~~~~~~~~~
## Define layout structure
SixByThreeMatrix <- matrix(data=c(1,2,3,4,
                                  1,5,6,7,
                                  1,8,9,10,
                                  1,11,12,13,
                                  1,14,15,16,
                                  1,17,18,19,
                                  20,21,21,21),nrow=7,ncol=4,byrow=TRUE)
SixByThreeLayout <- layout(mat=SixByThreeMatrix,widths=c(0.2,rep(1,3)),heights=c(1,1,1,1,1,1,0.2))
layout.show(n=21)

#~~~~~~~~~~~~~~~~~~
# Create 5%, median, 95% extractor function
plot.ci.extract <- function(x, scenario, group, stat.col) {
  list.index <- grep(pattern = scenario, x = names(x))  
  sapply(list.index, function(i) {x[[i]][group, stat.col]} )
}

# Inputs for extractor function
x <- KMA473PopsGroups17RepeatedMixProofTests_Repeats5_EstimatesStats
scenario <- FisheryProofTestScenarioNames17RG[1]
groups <- KMA17GroupsPC
stat.col <- c("5%", "median", "95%")[2]

#~~~~~~~~~~~~~~~~~~
# Create a figure wrapup funtion
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")

library('gplots')
library('devEMF')

PlotProof.GCL <- function(x, scenario, groups, filedir, percent = TRUE) {
  
  MixtureProofTestProportions <- t(dget(file = "Objects/MixtureProofTestProportions17RG.txt"))
  
  if(percent){
    percent.multiplier <- 100
    x <- sapply(x, function(rpt) {rpt[, 1:5] * percent.multiplier}, simplify = FALSE)
    MixtureProofTestProportions <- MixtureProofTestProportions * percent.multiplier
    file = paste(filedir, "/", scenario, ".percent.emf", sep = "")
    y.label = "Percent"
  } else {
    percent.multiplier <- 1
    file = paste(filedir, "/", scenario, ".emf", sep = "")
    y.label = "Proportion"
  }
  
  emf(file = file, width = 6.5, height = 7.5, family = "Times")
  #png(file = paste(filedir, "/", scenario, ".png", sep = ""), width = 6.5, height = 7.5, units = "in", res = 1200, pointsize = 4, family = "Times")
  layout(mat = SixByThreeMatrix, widths = c(0.3, rep(1, 3)), heights = c(1, 1, 1, 1, 1, 1, 0.4))
  
  
  # Plot 1 - y-axis title
  par(mar = c(0, 0, 0, 0))
  par(srt = 90)
  plot.new()
  text(labels = y.label, x = 0.2, y = 0.5, adj = c(0.5, 0.5), font = 1, cex = 2)
  
  # Plot 2-18 - Groups
  sapply(seq(groups), function(i) {
    par(mar = c(0, 0, 0, 0))
    par(srt = 0)
    if(i %in% c(1, 4, 7, 10, 13)) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 5.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", cex = 1.5)
    }
    
    if(i == 16) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 5.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', cex = 1.5)
    }
    
    if(i == 17) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 5.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', yaxt = "n", cex = 1.5)
    }
    
    if(i %in% c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15)) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 5.5), ylim = c(0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", yaxt = "n", cex = 1.5)
    }
    
    abline(h = MixtureProofTestProportions[i, scenario], col = "red", lwd = 2)
    abline(h = 0, lwd = 1, lty = 2)
    text(KMA17GroupsPC[i], x = 5, y = 0.6 * percent.multiplier, adj = 1, cex = 1)
  } )
  
  # Plot 19 - Blank Corner
  plotCI(x = rep(50, 5),
         ui = rep(60, 5),
         li = rep(40, 5),
         xlim = c(0.5, 5.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', yaxt = "n", cex = 1.5, col = "white")  

  # Plot 20 - Blank Corner
  plot.new()
  
  ## Plot 21 - x-axis title
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(labels = "Replicate", x = 0.5, y = 0.25, adj = c(0.5, 0.5), font = 1, cex = 2)
  
  dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# loci46 KarlukAyakulikSplit
sapply(FisheryProofTestScenarioNames17RG, function(scenario){
  PlotProof.GCL(x = dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5_EstimatesStats.txt"), 
                scenario = scenario, groups = KMA17GroupsPC,
                filedir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/Figures",
                percent = TRUE)
})




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# More proof tests, 5 additional repeats for a total of 10, 1 chain each, 40K

# Create mixture proof tests
sapply(FisheryProofTestScenarioNames17RG, function(Scenario) {
  sapply(6:10, function(Rpt) {
    assign(x = paste(Scenario, Rpt, "40KProof", sep = ''), 
           value = ProofTest.GCL(sillyvec = KMA473Pops, loci = loci46, groupnames = KMA17GroupsPC, groupvec = KMA473PopsGroupVec17,
                                 samplesize = MixtureProofTest17RG.SampleSize[Scenario, ], prefix = paste(Scenario, Rpt, sep = ''), 
                                 dir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit", prprtnl = TRUE, type = "BAYES",
                                 suffix = '', nreps = 40000, nchains = 1, priorvec = KMA473Pops17FlatPrior, initmat = KMA473PopsInits, 
                                 thin = c(1, 1, 100), switches = "F T F T T T F"), pos = 1)
  } )
} ); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dput proof test objects
objects(pattern = "Proof$")
KMA473PopsGroups17RepeatedMixProofTests40KRound2 <- paste(rep(FisheryProofTestScenarioNames17RG, each = 5), 6:10, "40K", sep = '')
dput(x = KMA473PopsGroups17RepeatedMixProofTests40KRound2, file = "Objects/KMA473PopsGroups17RepeatedMixProofTests40KRound2.txt")

# dir.create("MixtureProofTests objects/RG17")
invisible(sapply(KMA473PopsGroups17RepeatedMixProofTests40KRound2, function(proof) {
  dput(x = get(paste(proof, "Proof", sep = "")), file = paste("MixtureProofTests objects/RG17/", proof, "Proof.txt", sep = ""))} ))

# Grab proof test values
invisible(sapply(KMA473PopsGroups17RepeatedMixProofTests, function(proof) {
  assign(x = paste(proof, "40KProof", sep = ""), 
         value = dget(file = paste("MixtureProofTests objects/RG17/", proof, "40KProof.txt", sep = "")), pos = 1)
} ))


# Determine where fish are coming from for proof tests
str(early.Ayakulik540KProof, max.level = 1)

proof <- early.Ayakulik540KProof

sapply(proof[c("Frazer", "Ayakulik Early")], function(RG) {
  unlist(sapply(RG, function(pop) {
    if(!is.null(pop)) {length(pop)}
  } ))
} )



sapply(paste(KMA473PopsGroups17RepeatedMixProofTests40K, "40KProof", sep = ''), function(proof) {
  sapply(get(proof)[c("Frazer", "Ayakulik Early", "Ayakulik Late")], function(RG) {
    unlist(sapply(RG, function(pop) {
      if(!is.null(pop)) {length(pop)}
    } ))
  } )
}, simplify = FALSE )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create BAYES.output folders
invisible(sapply(paste(rep(FisheryProofTestScenarioNames17RG, each = 5), 6:10, sep = ''), function(proof) {dir <- paste(getwd(), "/BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.output/", proof, sep = "")
dir.create(dir) }))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pull posterior with burnin
KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA17GroupsPC), groupnames = KMA17GroupsPC, 
  maindir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.output", 
  mixvec = paste(rep(FisheryProofTestScenarioNames17RG, each = 5), 6:10, sep = ''), prior = "",  
  ext = "RGN", nchains = 1, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_Estimates)


# Dput 1) estimates stats + posterior output & 2) estimates stats
# dir.create("Estimates objects/loci46 KarlukAyakulikSplit")
dput(KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_Estimates, file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_Estimates.txt")
dput(KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_Estimates$Stats, file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_EstimatesStats.txt")

KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_Estimates <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_Estimates.txt")
KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_EstimatesStats <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2 ## No GR for 1 chain....

# Quick look at raw posterior output
str(KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_Estimates$Output)

PlotPosterior(mixvec = paste(rep(FisheryProofTestScenarioNames17RG, each = 5), 6:10, sep = ''), 
              output = KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_Estimates$Output, 
              groups = KMA17GroupsPC, 
              colors = Colors17, 
              header = NULL, set.mfrow = c(6, 3), thin = 1, chains = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SEDM style proof plots



#~~~~~~~~~~~~~~~~~~
## Define layout structure
SixByThreeMatrix <- matrix(data=c(1,2,3,4,
                                  1,5,6,7,
                                  1,8,9,10,
                                  1,11,12,13,
                                  1,14,15,16,
                                  1,17,18,19,
                                  20,21,21,21),nrow=7,ncol=4,byrow=TRUE)
SixByThreeLayout <- layout(mat=SixByThreeMatrix,widths=c(0.2,rep(1,3)),heights=c(1,1,1,1,1,1,0.2))
layout.show(n=21)

#~~~~~~~~~~~~~~~~~~
# Create 5%, median, 95% extractor function
plot.ci.extract <- function(x, scenario, group, stat.col) {
  list.index <- grep(pattern = scenario, x = names(x))  
  sapply(list.index, function(i) {x[[i]][group, stat.col]} )
}

# Inputs for extractor function
x <- KMA473PopsGroups17RepeatedMixProofTests_Repeats5_EstimatesStats
scenario <- FisheryProofTestScenarioNames17RG[1]
groups <- KMA17GroupsPC
stat.col <- c("5%", "median", "95%")[2]

#~~~~~~~~~~~~~~~~~~
# Create a figure wrapup funtion
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")

library('gplots')
library('devEMF')

PlotProof.GCL <- function(x, scenario, groups, filedir, percent = TRUE) {
  
  MixtureProofTestProportions <- t(dget(file = "Objects/MixtureProofTestProportions17RG.txt"))
  
  if(percent){
    percent.multiplier <- 100
    x <- sapply(x, function(rpt) {rpt[, 1:5] * percent.multiplier}, simplify = FALSE)
    MixtureProofTestProportions <- MixtureProofTestProportions * percent.multiplier
    file = paste(filedir, "/", scenario, "10repeats.percent.emf", sep = "")
    y.label = "Percent"
  } else {
    percent.multiplier <- 1
    file = paste(filedir, "/", scenario, "10repeats.emf", sep = "")
    y.label = "Proportion"
  }
  
  emf(file = file, width = 6.5, height = 7.5, family = "Times")
  #png(file = paste(filedir, "/", scenario, ".png", sep = ""), width = 6.5, height = 7.5, units = "in", res = 1200, pointsize = 4, family = "Times")
  layout(mat = SixByThreeMatrix, widths = c(0.3, rep(1, 3)), heights = c(1, 1, 1, 1, 1, 1, 0.4))
  
  
  # Plot 1 - y-axis title
  par(mar = c(0, 0, 0, 0))
  par(srt = 90)
  plot.new()
  text(labels = y.label, x = 0.2, y = 0.5, adj = c(0.5, 0.5), font = 1, cex = 2)
  
  # Plot 2-18 - Groups
  sapply(seq(groups), function(i) {
    par(mar = c(0, 0, 0, 0))
    par(srt = 0)
    if(i %in% c(1, 4, 7, 10, 13)) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 10.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", cex = 1.5)
    }
    
    if(i == 16) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 10.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', cex = 1.5)
    }
    
    if(i == 17) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 10.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', yaxt = "n", cex = 1.5)
    }
    
    if(i %in% c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15)) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 10.5), ylim = c(0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", yaxt = "n", cex = 1.5)
    }
    
    abline(h = MixtureProofTestProportions[i, scenario], col = "red", lwd = 2)
    abline(h = 0, lwd = 1, lty = 2)
    text(KMA17GroupsPC[i], x = 1, y = 0.6 * percent.multiplier, adj = 0, cex = 1)
  } )
  
  # Plot 19 - Blank Corner
  plotCI(x = rep(50, 5),
         ui = rep(60, 5),
         li = rep(40, 5),
         xlim = c(0.5, 10.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', yaxt = "n", cex = 1.5, col = "white")  
  
  # Plot 20 - Blank Corner
  plot.new()
  
  ## Plot 21 - x-axis title
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(labels = "Replicate", x = 0.5, y = 0.25, adj = c(0.5, 0.5), font = 1, cex = 2)
  
  dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA473PopsGroups17RepeatedMixProofTests_Repeats5_EstimatesStats <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5_EstimatesStats.txt")
KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_EstimatesStats <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_EstimatesStats.txt")

KMA473PopsGroups17RepeatedMixProofTests_Repeats10_EstimatesStats <- c(KMA473PopsGroups17RepeatedMixProofTests_Repeats5_EstimatesStats, KMA473PopsGroups17RepeatedMixProofTests_Repeats5.2_EstimatesStats)
dput(x = KMA473PopsGroups17RepeatedMixProofTests_Repeats10_EstimatesStats, file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats10_EstimatesStats.txt")
KMA473PopsGroups17RepeatedMixProofTests_Repeats10_EstimatesStats <- dget(file = "Estimates objects/loci46 KarlukAyakulikSplit/KMA473PopsGroups17RepeatedMixProofTests_Repeats10_EstimatesStats.txt")

# loci46 KarlukAyakulikSplit
sapply(FisheryProofTestScenarioNames17RG, function(scenario){
  PlotProof.GCL(x = KMA473PopsGroups17RepeatedMixProofTests_Repeats10_EstimatesStats, 
                scenario = scenario, groups = KMA17GroupsPC,
                filedir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/Figures",
                percent = TRUE)
})



BiasRMSE.GCL <- function(stats, scenario, proportions, estimator = "median", percent = TRUE) {
  if(percent) {
    stats <- sapply(stats, function(rpt) {rpt[, 1:5] * 100}, simplify = FALSE)
    proportions <- proportions * 100
  }  
  stats.nums <- grep(pattern = scenario, x = names(stats))
  estimator.values <- sapply(stats[stats.nums], function(rpt) {rpt[, estimator]} )
  rmse <- sqrt(rowSums((estimator.values - proportions[, scenario])^2) / length(stats.nums))
  bias <- apply((estimator.values - proportions[, scenario]), 1, mean)
  average <- apply(estimator.values, 1, mean)
  ci.width <- apply(sapply(stats[stats.nums], function(rpt) {rpt[, "95%"]} ) - sapply(stats[stats.nums], function(rpt) {rpt[, "5%"]} ), 1, mean)
  
  return(cbind(average = average, bias = bias, rmse = rmse, ci.width = ci.width))
}




Proof.BiasRMSE.46FrazerAyakulikloci.list <- 
  sapply(FisheryProofTestScenarioNames17RG, function(scenario) {
    BiasRMSE.GCL(stats = KMA473PopsGroups17RepeatedMixProofTests_Repeats10_EstimatesStats, 
                 scenario = scenario, proportions = t(dget(file = "Objects/MixtureProofTestProportions17RG.txt")),
                 estimator = "median")}, simplify = FALSE )


# across all scenarios
x.array.46 <- array(unlist(Proof.BiasRMSE.46FrazerAyakulikloci.list), dim = c(17,4,4), dimnames = list(KMA17GroupsPC, c("average", "bias", "rmse", "ci.width"), FisheryProofTestScenarioNames17RG))
round(apply(X = x.array.46[, , ], MARGIN = 1:2, FUN = mean)[, 2:4], 1)
round(apply(X = x.array.46[, , ], MARGIN = 1:2, FUN = min)[, 2:4], 1)
round(apply(X = x.array.46[, , ], MARGIN = 1:2, FUN = max)[, 2:4], 1)

require(gplots)
lapply(FisheryProofTestScenarioNames17RG, function(scenario) {
  scenario.mat <- Proof.BiasRMSE.46FrazerAyakulikloci.list[[scenario]]
  barcolor <- ifelse(scenario.mat[, "bias"] > 0, "green", "red")
  locations <- barplot2(height = scenario.mat[, "bias"], col = barcolor, ylim = c(-5, 5), ylab = "Bias", xlab = "Reporting Group", xaxt = 'n', main = scenario, cex.main = 1.5)
  text(x = locations[, 1], y = -5.5, srt = 90, labels = KMA17GroupsPC, adj = 0)
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make dotplot of medians (green/red), with transparency with line showing mean of medians

str(KMA473PopsGroups17RepeatedMixProofTests_Repeats10_EstimatesStats)

Bias.GCL <- function(stats, scenario, proportions, estimator = "median", percent = TRUE) {
  if(percent) {
    stats <- sapply(stats, function(rpt) {rpt[, 1:5] * 100}, simplify = FALSE)
    proportions <- proportions * 100
  }  
  stats.nums <- grep(pattern = scenario, x = names(stats))
  estimator.values <- sapply(stats[stats.nums], function(rpt) {rpt[, estimator]} )
  bias <- estimator.values - proportions[, scenario]
  
  return(bias)
}

Proof.Bias.46FrazerAyakulikloci.list <- 
  sapply(FisheryProofTestScenarioNames17RG, function(scenario) {
    Bias.GCL(stats = KMA473PopsGroups17RepeatedMixProofTests_Repeats10_EstimatesStats, 
             scenario = scenario, proportions = t(dget(file = "Objects/MixtureProofTestProportions17RG.txt")),
             estimator = "median")}, simplify = FALSE )

str(Proof.Bias.46FrazerAyakulikloci.list)



estimates <- Proof.Bias.46FrazerAyakulikloci.list
scenario <- FisheryProofTestScenarioNames17RG[1]
proportions <- t(dget(file = "Objects/MixtureProofTestProportions17RG.txt"))
scenario.names.PC <- setNames(object = c("June Ayakulik (only early-run fish)",
                                         "July Alitak (only early-run fish)",
                                         "July Alitak (50/50 early-run/late-run fish)",
                                         "July Alitak (only late-run fish)"),
                              nm = FisheryProofTestScenarioNames17RG)

sapply(FisheryProofTestScenarioNames17RG, function(scenario) {
  png(file = paste("BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/Figures/BiasBoxplot_", scenario, ".png", sep = ''), width = 7, height = 7, units = "in", res = 600, family = "Times")
  # Make as vector
  bias.mat <- t(estimates[[scenario]])
  bias.vec <- c(bias.mat)
  range(bias.vec)
  
  # Assign colors
  bias.col <- rep("darkgreen", length(bias.vec))
  bias.col[bias.vec < 0] <- "darkred"
  bias.col <- adjustcolor(col = bias.col, alpha.f = 0.5)
  
  # Plot
  # par(mar = c(9.1, 4.1, 2.1, 1.1))
  # plot(x = rep(1:17, each = 10), y = bias.vec, pch = 16, col = bias.col, cex = 2.5, ylim = c(-10, 10), axes = FALSE, xlab = "", ylab = "")
  # axis(side = 1, at = 1:17, labels = NA)
  # axis(side = 2)
  # text(x = 1:17, y = rep(-11.5, 17), labels = KMA17GroupsPC, adj = 1, srt = 45, xpd = TRUE)
  # mtext(text = "Reporting Group", side = 1, line = 7, cex = 1.5)
  # mtext(text = "Bias", side = 2, line = 3, cex = 1.5)
  # mtext(text = scenario, side = 3, cex = 2)
  
  # Boxplot
  par(mar = c(9.1, 4.6, 2.1, 1.1))
  boxplot(bias.mat, axes = FALSE, ylim = c(-15, 10))
  abline(h = 0, lwd = 4, lty = 3)
  points(x = rep(1:17, each = 10), y = bias.vec, pch = 16, col = bias.col, cex = 2.5)
  axis(side = 1, at = 1:17, labels = NA)
  axis(side = 2)
  text(x = 1:17, y = rep(-17, 17), labels = KMA17GroupsPC, adj = 1, srt = 45, xpd = TRUE)
  text(x = 1:17, y = rep(-14.5, 17), labels = proportions[, scenario] * 100)
  text(x = 17/2, y = -13.5, labels = "True Percentage in Mixture Scenario")
  mtext(text = "Reporting Group", side = 1, line = 7, cex = 1.5)
  mtext(text = "Bias Over 10 Replicates (Percentage)", side = 2, line = 3, cex = 1.5)
  mtext(text = scenario.names.PC[scenario], side = 3, cex = 2)
  dev.off()
})



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make tables of Bias, RMSE, etc.

str(Proof.BiasRMSE.46FrazerAyakulikloci.list)

PCFisheryProofTestScenarioNames17RG <- c("June Ayakulik (early-run)", "July Alitak (early-run)", "July Alitak (50/50 early/late-run)", "July Alitak (late-run)")
dput(x = PCFisheryProofTestScenarioNames17RG, file = "Objects/PCFisheryProofTestScenarioNames17RG.txt")

MasterProofTableMaker <- function(loci = "loci89", BiasRMSElist, table.file, percent = TRUE) {
  
  require(xlsx)
  
  if(percent) {
    digits = 1
  } else {
    digits = 2
  }
  
  x <- NULL
  x.mat <- suppressWarnings(matrix(data = seq(FisheryProofTestScenarioNames17RG), ncol = 2, byrow = TRUE))
  scenario <- FisheryProofTestScenarioNames17RG
  proportions <- dget(file = "Objects/MixtureProofTestProportions17RG.txt")
  if(percent) {proportions <- proportions * 100}
  
  for(i in seq(nrow(x.mat))){
    x <- rbind(x, cbind('', paste("Hypothetical", PCFisheryProofTestScenarioNames17RG[x.mat[i, 1]], "Scenario"), '', '', '', '', '', paste("Hypothetical", PCFisheryProofTestScenarioNames17RG[x.mat[i, 2]], "Scenario"), '', '', '', ''),
               cbind("Reporting Group", "True", "Average", "Bias", "RMSE", "CI Width", "", "True", "Average", "Bias", "RMSE", "CI Width"),
               cbind(KMA17GroupsPC,
                     formatC(x = proportions[x.mat[i, 1], ], digits = digits, format = "f"),
                     formatC(x = BiasRMSElist[[scenario[x.mat[i, 1]]]], digits = digits, format = "f"),
                     rep('', length(KMA17GroupsPC)),
                     formatC(x = proportions[x.mat[i, 2], ], digits = digits, format = "f"),
                     formatC(x = BiasRMSElist[[scenario[x.mat[i, 2]]]], digits = digits, format = "f"))
    )
  }
  suppressWarnings(write.xlsx(x = x, file = table.file, sheetName = paste("Fishery Proof Summary", loci), col.names = FALSE, row.names = FALSE, append = TRUE))
  print(apply(MixtureProofTest.SampleSize, 1, sum))
  
}

MasterProofTableMaker(loci = "loci46 17RG", BiasRMSElist = Proof.BiasRMSE.46FrazerAyakulikloci.list, table.file = "Estimates tables/KMA473PopsGroups15RepeatedProofTestsTables.xlsx", percent = TRUE)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Composite reporting group

KMA14GroupsPC <- c(KMA17GroupsPC[1:4], "Frazer / Ayakulik", "Karluk", KMA17GroupsPC[10:17])
dput(x = KMA14GroupsPC, file = "Objects/KMA14GroupsPC.txt")

KMA473PopsGroupVec14 <- as.numeric(gsub(pattern = 6, replacement = 5, x = KMA473PopsGroupVec15))
KMA473PopsGroupVec14 <- c(KMA473PopsGroupVec14[KMA473PopsGroupVec14 < 6], KMA473PopsGroupVec14[KMA473PopsGroupVec14 > 6] - 1)
dput(x = KMA473PopsGroupVec14, file = "Objects/KMA473PopsGroupVec14.txt")


KMA473PopsGroups14RepeatedMixProofTests_Repeats10_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = KMA473PopsGroupVec14, groupnames = KMA14GroupsPC, 
  maindir = "BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/BAYES.output", 
  mixvec = paste(rep(FisheryProofTestScenarioNames17RG, each = 10), 1:10, sep = ''), prior = "",  
  ext = "BOT", nchains = 1, burn = 0.5, alpha = 0.1, PosteriorOutput = FALSE); beep(5)

str(KMA473PopsGroups14RepeatedMixProofTests_Repeats10_Estimates)

fishery.proportions <- t(dget(file = "Objects/MixtureProofTestProportions17RG.txt"))
fishery.proportions <- rbind(fishery.proportions[1:4, ],
                             "Frazer / Ayakulik" = colSums(fishery.proportions[5:7, ]),
                             "Karluk" = colSums(fishery.proportions[8:9, ]),
                             fishery.proportions[10:17, ])



Proof.BiasRMSE.46FrazerAyakulikloci.14RG.list <- 
  sapply(FisheryProofTestScenarioNames17RG, function(scenario) {
    BiasRMSE.GCL(stats = KMA473PopsGroups14RepeatedMixProofTests_Repeats10_Estimates, 
                 scenario = scenario, proportions = fishery.proportions,
                 estimator = "median")}, simplify = FALSE )

x.array.46 <- array(unlist(Proof.BiasRMSE.46FrazerAyakulikloci.14RG.list), dim = c(14,4,4), dimnames = list(KMA14GroupsPC, c("average", "bias", "rmse", "ci.width"), FisheryProofTestScenarioNames17RG))
round(apply(X = x.array.46[, , ], MARGIN = 1:2, FUN = mean)[, 2:4], 1)
round(apply(X = x.array.46[, , ], MARGIN = 1:2, FUN = min)[, 2:4], 1)
round(apply(X = x.array.46[, , ], MARGIN = 1:2, FUN = max)[, 2:4], 1)





Bias.GCL <- function(stats, scenario, proportions, estimator = "median", percent = TRUE) {
  if(percent) {
    stats <- sapply(stats, function(rpt) {rpt[, 1:5] * 100}, simplify = FALSE)
    proportions <- proportions * 100
  }  
  stats.nums <- grep(pattern = scenario, x = names(stats))
  estimator.values <- sapply(stats[stats.nums], function(rpt) {rpt[, estimator]} )
  bias <- estimator.values - proportions[, scenario]
  
  return(bias)
}

Proof.Bias.46FrazerAyakulikloci.14RG.list <- 
  sapply(FisheryProofTestScenarioNames17RG, function(scenario) {
    Bias.GCL(stats = KMA473PopsGroups14RepeatedMixProofTests_Repeats10_Estimates, 
             scenario = scenario, proportions = fishery.proportions,
             estimator = "median")}, simplify = FALSE )


estimates <- Proof.Bias.46FrazerAyakulikloci.14RG.list
scenario <- FisheryProofTestScenarioNames17RG[1]
proportions <- fishery.proportions
scenario.names.PC <- setNames(object = c("June Ayakulik (only early-run fish)",
                                         "July Alitak (only early-run fish)",
                                         "July Alitak (50/50 early-run/late-run fish)",
                                         "July Alitak (only late-run fish)"),
                              nm = FisheryProofTestScenarioNames17RG)

sapply(FisheryProofTestScenarioNames17RG, function(scenario) {
  png(file = paste("BAYES/Mixture Proof Tests/loci46 KarlukAyakulikSplit/Figures/BiasBoxplot_14RG_", scenario, ".png", sep = ''), width = 7, height = 7, units = "in", res = 600, family = "Times")
  # Make as vector
  bias.mat <- t(estimates[[scenario]])
  bias.vec <- c(bias.mat)
  range(bias.vec)
  
  # Assign colors
  bias.col <- rep("darkgreen", length(bias.vec))
  bias.col[bias.vec < 0] <- "darkred"
  bias.col <- adjustcolor(col = bias.col, alpha.f = 0.5)
  
  # Plot
  # par(mar = c(9.1, 4.1, 2.1, 1.1))
  # plot(x = rep(1:17, each = 10), y = bias.vec, pch = 16, col = bias.col, cex = 2.5, ylim = c(-10, 10), axes = FALSE, xlab = "", ylab = "")
  # axis(side = 1, at = 1:17, labels = NA)
  # axis(side = 2)
  # text(x = 1:17, y = rep(-11.5, 17), labels = KMA17GroupsPC, adj = 1, srt = 45, xpd = TRUE)
  # mtext(text = "Reporting Group", side = 1, line = 7, cex = 1.5)
  # mtext(text = "Bias", side = 2, line = 3, cex = 1.5)
  # mtext(text = scenario, side = 3, cex = 2)
  
  # Boxplot
  par(mar = c(9.1, 4.6, 2.1, 1.1))
  boxplot(bias.mat, axes = FALSE, ylim = c(-15, 10))
  abline(h = 0, lwd = 4, lty = 3)
  points(x = rep(1:14, each = 10), y = bias.vec, pch = 16, col = bias.col, cex = 2.5)
  axis(side = 1, at = 1:14, labels = NA)
  axis(side = 2)
  text(x = 1:14, y = rep(-17, 14), labels = KMA14GroupsPC, adj = 1, srt = 45, xpd = TRUE)
  text(x = 1:14, y = rep(-14.5, 14), labels = proportions[, scenario] * 100)
  text(x = 14/2, y = -13.5, labels = "True Percentage in Mixture Scenario")
  mtext(text = "Reporting Group", side = 1, line = 7, cex = 1.5)
  mtext(text = "Bias Over 10 Replicates (Percentage)", side = 2, line = 3, cex = 1.5)
  mtext(text = scenario.names.PC[scenario], side = 3, cex = 2)
  dev.off()
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
sapply(Round2Mixtures80K_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Resolved

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

sapply(Round2Mixtures_2014_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # No issues


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 MSA files for BAYES 2014 Late Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify Strata to Run
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
                                                  "Karluk Section Middle July 26-August 25, 2014",
                                                  "Uganik Section Middle July 26-August 25, 2014",
                                                  "Uyak Section Middle July 26-August 25, 2014"), 
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
#### Summarize Round 3 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures80K_2014_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                                 maindir = "BAYES/2014-2015 Mixtures 46loci/Output/80K Iterations", 
                                                                 mixvec = Round3Mixtures_2014[3], prior = "",  
                                                                 ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round3Mixtures80K_2014_Estimates, file = "Estimates objects/Round3Mixtures80K_2014_Estimates.txt")
dput(Round3Mixtures80K_2014_Estimates$Stats, file = "Estimates objects/Round3Mixtures80K_2014_EstimatesStats.txt")

Round3Mixtures80K_2014_Estimates <- dget(file = "Estimates objects/Round3Mixtures80K_2014_Estimates.txt")
Round3Mixtures80K_2014_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures80K_2014_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round3Mixtures80K_2014_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round3Mixtures80K_2014_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Resolved
require(gplots)
sapply(Round3Mixtures_2014[c(3, 4, 2, 1)], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round3Mixtures80K_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round3Mixtures80K_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round3Mixtures80K_2014_Estimates$Output)

PlotPosterior(mixvec = Round3Mixtures_2014[3], output = Round3Mixtures80K_2014_Estimates$Output, 
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round3Mixtures_2014_Header[3], set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 3 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round3Mixtures_2014[3], estimatesstats = Round3Mixtures80K_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round3Mixtures_2014_Header[3])

## Quick barplot
QuickBarplot(mixvec = Round3Mixtures_2014[3], estimatesstats = Round3Mixtures80K_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round3Mixtures_2014_Header[3])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round3Mixtures_2014[3], estimates = Round3Mixtures80K_2014_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round3Mixtures_2014_Header)
rm(Round3Mixtures80K_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round3Mixtures_2014_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures_2014_EstimatesStats.txt")

sapply(Round3Mixtures_2014[3], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round3Mixtures_2014_EstimatesStats[[Mix]][, "median"],
                                                 Round3Mixtures80K_2014_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round3Mixtures_2014_EstimatesStats[[Mix]][, "5%"],
                                               Round3Mixtures80K_2014_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round3Mixtures_2014_EstimatesStats[[Mix]][, "95%"],
                                               Round3Mixtures80K_2014_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round3Mixtures_2014_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA15GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 3 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SAYAKC14_3_Late, SKARLC14_3_Late, and SUYAKC14_3_Late as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from SUGANC14_3_Late
Round3Mixtures_2014_Estimates <- dget(file = "Estimates objects/Round3Mixtures_2014_Estimates.txt")

sapply(Round3Mixtures_2014_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC14_3_Late
sapply(Round3Mixtures80K_2014_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Resolved

str(Round3Mixtures_2014_Estimates)
str(Round3Mixtures80K_2014_Estimates)

Round3Mixtures_2014_Estimates_Final <- list(Stats = c(Round3Mixtures_2014_Estimates$Stats[1:2],
                                                      Round3Mixtures80K_2014_Estimates$Stats[1],
                                                      Round3Mixtures_2014_Estimates$Stats[4]),
                                            Output = c(Round3Mixtures_2014_Estimates$Output[1:2],
                                                       Round3Mixtures80K_2014_Estimates$Output[1],
                                                       Round3Mixtures_2014_Estimates$Output[4]))
str(Round3Mixtures_2014_Estimates_Final)
dput(x = Round3Mixtures_2014_Estimates_Final, file = "Estimates objects/Round3Mixtures_2014_Estimates_Final.txt")
dput(x = Round3Mixtures_2014_Estimates_Final$Stats, file = "Estimates objects/Round3Mixtures_2014_EstimatesStats_Final.txt")

Round3Mixtures_2014_Estimates_Final <- dget(file = "Estimates objects/Round3Mixtures_2014_Estimates_Final.txt")
Round3Mixtures_2014_EstimatesStats_Final <- dget(file = "Estimates objects/Round3Mixtures_2014_EstimatesStats_Final.txt")

sapply(Round3Mixtures_2014_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # No issues


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 3 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2014_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci/Output/40K Iterations", 
                                                                   mixvec = Round3Mixtures_2014, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round3Mixtures_2014_31RG_Estimates)
QuickBarplot(mixvec = Round3Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round3Mixtures_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round3Mixtures_2014_Header[c(3, 4, 2, 1)])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 4 MSA files for BAYES 2014 LateLate + Alitak Late Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Indetify Strata to Run
KMA2014Strata_4_LateLate <- grep(pattern = "4_LateLate", x = KMA2014Strata, value = TRUE)
Round4Mixtures_2014 <- c("SALITC14_3_Late", KMA2014Strata_4_LateLate)  # Add SALITC14_3_Late now that we've decided to keep 40K results
dput(x = Round4Mixtures_2014, file = "Objects/Round4Mixtures_2014.txt")

# Create rolling prior based on Round 1 and 3 estimates
str(Round1Mixtures_2014_EstimatesStats)
sapply(Round1Mixtures_2014_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})
str(Round3Mixtures_2014_EstimatesStats_Final)
sapply(Round3Mixtures_2014_EstimatesStats_Final, function(Mix) {table(Mix[, "GR"] > 1.2)})

str(c(Round1Mixtures_2014_EstimatesStats["SALITC14_2_Middle"], Round3Mixtures_2014_EstimatesStats_Final[c("SKARLC14_3_Late", "SUGANC14_3_Late", "SUYAKC14_3_Late")]))

Round4Mixtures_2014_Prior <- sapply(c(Round1Mixtures_2014_EstimatesStats["SALITC14_2_Middle"], Round3Mixtures_2014_EstimatesStats_Final[c("SKARLC14_3_Late", "SUGANC14_3_Late", "SUYAKC14_3_Late")]),
                                    function(Mix) {Prior.GCL(groupvec = KMA473PopsGroupVec15, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round4Mixtures_2014_Prior) <- gsub(pattern = "3_Late", replacement = "4_LateLate", x = names(Round4Mixtures_2014_Prior))  # This changes the names
names(Round4Mixtures_2014_Prior) <- gsub(pattern = "2_Middle", replacement = "3_Late", x = names(Round4Mixtures_2014_Prior))  # This changes the names
dput(x = Round4Mixtures_2014_Prior, file = "Objects/Round4Mixtures_2014_Prior.txt")
str(Round4Mixtures_2014_Prior)


## Dumping Mixture files
sapply(Round4Mixtures_2014, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round4Mixtures_2014, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round4Mixtures_2014_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round4Mixtures_2014, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 4 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round4Mixtures_2014_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
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
sapply(Round4Mixtures_2014[c(3, 4, 2, 1)], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round4Mixtures_2014_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round4Mixtures_2014_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round4Mixtures_2014_Estimates$Output)
Round4Mixtures_2014_Header <- setNames(object = c("Alitak Section Late July 26-August 29, 2014",
                                                  "Karluk Section LateLate August 25-29, 2014",
                                                  "Uganik Section LateLate August 25-29, 2014",
                                                  "Uyak Section LateLate August 25-29, 2014"), 
                                       nm = Round4Mixtures_2014)
dput(x = Round4Mixtures_2014_Header, file = "Objects/Round4Mixtures_2014_Header.txt")

PlotPosterior(mixvec = Round4Mixtures_2014, output = Round4Mixtures_2014_Estimates$Output, 
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round4Mixtures_2014_Header, set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 4 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round4Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round4Mixtures_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round4Mixtures_2014_Header)

## Quick barplot
QuickBarplot(mixvec = Round4Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round4Mixtures_2014_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round4Mixtures_2014_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round4Mixtures_2014[c(3, 4, 2, 1)], estimates = Round4Mixtures_2014_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round4Mixtures_2014_Header)
rm(Round4Mixtures_2014_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 4 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round4Mixtures_2014_31RG_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci/Output/40K Iterations", 
  mixvec = Round4Mixtures_2014, prior = "",  
  ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round4Mixtures_2014_31RG_Estimates)
QuickBarplot(mixvec = Round4Mixtures_2014[c(3, 4, 2, 1)], estimatesstats = Round4Mixtures_2014_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round4Mixtures_2014_Header[c(3, 4, 2, 1)])


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

max(HarvestByStrata2014, HarvestByStrata2015, na.rm = TRUE)


# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
# Need to test first to make sure it will work with mixtures of different length (Uganik 3 had to be run out to 80K)
sapply(rownames(HarvestByStrata2014)[!is.na(HarvestByStrata2014[, 4])], function(geomix) {
  assign(x = paste(geomix, "_3_Late_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
           mixvec = c(paste(geomix, "_3_Late", sep = ''), paste(geomix, "_4_LateLate", sep = '')), 
           catchvec = HarvestByStrata2014[geomix, 3:4], newname = paste(geomix, "_3_Late_Stratified", sep = ''), priorname = '',
           ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1), 
         pos = 1)
  dput(x = get(paste(geomix, "_3_Late_Stratified", sep = '')), file = paste("Estimates objects/", geomix, "_3_Late_Stratified.txt", sep = ''))
} ); beep(5)

dimnames(SKARLC14_3_Late_Stratified$Summary)
dimnames(Round4Mixtures_2014_Estimates$Stats$SKARLC14_4_LateLate)

# Look at Late and LateLate separately and then compare to stratified estimate
Round1Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")
Round2Mixtures_2014_EstimatesStats_Final <- dget("Estimates objects/Round2Mixtures_2014_EstimatesStats_Final.txt")
Round3Mixtures_2014_EstimatesStats_Final <- dget("Estimates objects/Round3Mixtures_2014_EstimatesStats_Final.txt")
Round4Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round4Mixtures_2014_EstimatesStats.txt")

KMA2014Strata_EstimatesStats <- c(Round1Mixtures_2014_EstimatesStats, 
                                  Round2Mixtures_2014_EstimatesStats_Final, 
                                  Round3Mixtures_2014_EstimatesStats_Final, 
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
Groups <- KMA15GroupsPC
Groups2Rows <- KMA15GroupsPC2Rows
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
  round(get(paste(geomix, "14_3_Late_Stratified", sep = ''))$Summary[, "MEDIAN"] * 100, 1)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Final Estimates Objects for Each Temporal Strata ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Round1Mixtures_2014_EstimatesStats <- dget("Estimates objects/Round1Mixtures_2014_EstimatesStats.txt")
Round2Mixtures_2014_EstimatesStats_Final <- dget("Estimates objects/Round2Mixtures_2014_EstimatesStats_Final.txt")
Round3Mixtures_2014_EstimatesStats_Final <- dget("Estimates objects/Round3Mixtures_2014_EstimatesStats_Final.txt")
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
str(Round3Mixtures_2014_EstimatesStats_Final)
str(Round4Mixtures_2014_EstimatesStats)
str(SKARLC14_3_Late_Stratified$Summary)

# dealing with different output from CustomCombineBAYESOutput and StratifiedEstimator
colnames(SKARLC14_3_Late_Stratified$Summary) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")
colnames(SUGANC14_3_Late_Stratified$Summary) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")
colnames(SUYAKC14_3_Late_Stratified$Summary) <- c("mean", "sd", "5%", "median", "95%", "P=0", "GR")

SKARLC14_3_Late_Stratified$Summary <- SKARLC14_3_Late_Stratified$Summary[, c("mean", "sd", "median", "5%", "95%", "P=0", "GR")]
SUGANC14_3_Late_Stratified$Summary <- SUGANC14_3_Late_Stratified$Summary[, c("mean", "sd", "median", "5%", "95%", "P=0", "GR")]
SUYAKC14_3_Late_Stratified$Summary <- SUYAKC14_3_Late_Stratified$Summary[, c("mean", "sd", "median", "5%", "95%", "P=0", "GR")]

KMA2014Strata_3_Late_EstimatesStats <- c(Round4Mixtures_2014_EstimatesStats[1],
                                         Round3Mixtures_2014_EstimatesStats_Final[1],
                                         list("SKARLC14_3_Late" = SKARLC14_3_Late_Stratified$Summary),
                                         list("SUGANC14_3_Late" = SUGANC14_3_Late_Stratified$Summary),
                                         list("SUYAKC14_3_Late" = SUYAKC14_3_Late_Stratified$Summary))

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
Groups <- KMA15GroupsPC
Groups2Rows <- KMA15GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5

sapply(names(TempMix14[c(4, 5, 3, 2, 1)]), function(geomix) {
  # emf(file = paste("Figures/2014/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "sans", bg = "white")
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
  # dev.off()
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
Groups <- KMA15GroupsPC
Groups2Rows <- KMA15GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5
ymax <- 140000

sapply(names(TempMix14[c(4, 5, 3, 2, 1)]), function(geomix) {
  # emf(file = paste("Figures/2014/", geomix, "Harvest.emf", sep = ''), width = 8.5, height = 6.5, family = "sans", bg = "white")
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
  # dev.off()
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
           groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
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
  Stats <- get(paste(geomix, "_Annual_Stratified", sep = ''))$Summary[, c("MEAN", "SD", "MEDIAN", "5%", "95%", "P0", "GR")]
  colnames(Stats) <- c("mean", "sd", "median", "5%", "95%", "P=0", "GR")
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
  text(x = colMeans(Barplot), y = -1, labels = KMA15GroupsPC2Rows, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)
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
  Prior.GCL(groupvec = KMA473PopsGroupVec15, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round1Mixtures_2015_Prior) <- gsub(pattern = "C14", replacement = "C15", x = names(Round1Mixtures_2015_Prior))  # This changes the names
# Use a flat prior for SALITC15_1_Early as there was no 2014 Early strata for Alitak
Round1Mixtures_2015_Prior <- c(list("SALITC15_1_Early" = KMA473Pops15FlatPrior), Round1Mixtures_2015_Prior)
Round1Mixtures_2015 %in% names(Round1Mixtures_2015_Prior)

dput(x = Round1Mixtures_2015_Prior, file = "Objects/Round1Mixtures_2015_Prior.txt")
str(Round1Mixtures_2015_Prior)


## Dumping Mixture files
sapply(Round1Mixtures_2015, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round1Mixtures_2015, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round1Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round1Mixtures_2015, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 2015 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC ,
  maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
  mixvec = Round1Mixtures_2015, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round1Mixtures_2015_Estimates, file = "Estimates objects/Round1Mixtures_2015_Estimates.txt")
dput(Round1Mixtures_2015_Estimates$Stats, file = "Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")

Round1Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2015_Estimates.txt")
Round1Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round1Mixtures_2015_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round1Mixtures_2015_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC14_3_Late
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round1Mixtures_2015[c(4, 5, 3, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round1Mixtures_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round1Mixtures_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
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
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round1Mixtures_2015_Header, set.mfrow = c(5, 3), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 2015 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round1Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round1Mixtures_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round1Mixtures_2015_Header)

## Quick barplot
QuickBarplot(mixvec = Round1Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round1Mixtures_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round1Mixtures_2015_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round1Mixtures_2015[c(4, 5, 3, 2, 1)], estimates = Round1Mixtures_2015_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round1Mixtures_2015_Header)
rm(Round1Mixtures_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 1 2015 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures_2015_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
                                                                   mixvec = Round1Mixtures_2015, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round1Mixtures_2015_31RG_Estimates)
QuickBarplot(mixvec = Round1Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round1Mixtures_2015_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round1Mixtures_2015_Header)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 1 2015 Repeated MSA files for BAYES 80K Iterations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping out a new control file to re-analyze SUGANC15_1_Early and SUYAKC15_1_Early Round 1 mixtures with a 80K posterior as opposed to our standard 40K
## Perhaps it will converge on one of the modes for South of Cape Suckling?

# Just use original mixture, no need to create a new one!

## Dumping a New Control file that goes out further!

sapply(Round1Mixtures_2015[4:5], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round1Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round1Mixtures_2015[4:5], function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})
# Will rename later, but leave as is for now for BayesCopyPaste.GCL


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 1 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round1Mixtures80K_2015_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
                                                                 maindir = "BAYES/2014-2015 Mixtures 46loci/Output/80K Iterations", 
                                                                 mixvec = Round1Mixtures_2015[4:5], prior = "",  
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
sapply(Round1Mixtures_2015[4:5], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round1Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round1Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round1Mixtures80K_2015_Estimates$Output)

PlotPosterior(mixvec = Round1Mixtures_2015[4:5], output = Round1Mixtures80K_2015_Estimates$Output, 
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round1Mixtures_2015_Header[4:5], set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 1 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round1Mixtures_2015[4:5], estimatesstats = Round1Mixtures80K_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round1Mixtures_2015_Header[3])

## Quick barplot
QuickBarplot(mixvec = Round1Mixtures_2015[4:5], estimatesstats = Round1Mixtures80K_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round1Mixtures_2015_Header[4:5])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round1Mixtures_2015[4:5], estimates = Round1Mixtures80K_2015_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round1Mixtures_2015_Header)
rm(Round1Mixtures80K_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round1Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round1Mixtures_2015_EstimatesStats.txt")

sapply(Round1Mixtures_2015[4:5], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round1Mixtures_2015_EstimatesStats[[Mix]][, "median"],
                                                 Round1Mixtures80K_2015_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round1Mixtures_2015_EstimatesStats[[Mix]][, "5%"],
                                               Round1Mixtures80K_2015_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round1Mixtures_2015_EstimatesStats[[Mix]][, "95%"],
                                               Round1Mixtures80K_2015_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round1Mixtures_2015_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA15GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 1 2015 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SALITC15_1_Early, SAYAKC15_1_Early, and SKARLC15_1_Early as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from SUGANC15_1_Early and SUYAKC15_1_Early, even though SUGANC15_1_Early did NOT converge with GR < 1.2 with 80K
Round1Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round1Mixtures_2015_Estimates.txt")

sapply(Round1Mixtures_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC15_1_Early and SUYAKC15_1_Early
sapply(Round1Mixtures80K_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUYAKC15_1_Early resolved, but SUGANC15_1_Early still has issues

str(Round1Mixtures_2015_Estimates)
str(Round1Mixtures80K_2015_Estimates)

Round1Mixtures_2015_Estimates_Final <- list(Stats = c(Round1Mixtures_2015_Estimates$Stats[Round1Mixtures_2015[1:3]],
                                                      Round1Mixtures80K_2015_Estimates$Stats[Round1Mixtures_2015[4:5]]),
                                            Output = c(Round1Mixtures_2015_Estimates$Output[Round1Mixtures_2015[1:3]],
                                                       Round1Mixtures80K_2015_Estimates$Output[Round1Mixtures_2015[4:5]]))
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
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
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
  Prior.GCL(groupvec = KMA473PopsGroupVec15, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round2Mixtures_2015_Prior) <- gsub(pattern = "1_Early", replacement = "2_Middle", 
                                         x = names(Round2Mixtures_2015_Prior))  # This changes the names
dput(x = Round2Mixtures_2015_Prior, file = "Objects/Round2Mixtures_2015_Prior.txt")
str(Round2Mixtures_2015_Prior)


## Dumping Mixture files
sapply(Round2Mixtures_2015, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round2Mixtures_2015, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round2Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2015, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 2015 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC ,
  maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
  mixvec = Round2Mixtures_2015, prior = "",  
  ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)

# Dput 1) estimates stats + posterior output & 2) estimates stats
dput(Round2Mixtures_2015_Estimates, file = "Estimates objects/Round2Mixtures_2015_Estimates.txt")
dput(Round2Mixtures_2015_Estimates$Stats, file = "Estimates objects/Round2Mixtures_2015_EstimatesStats.txt")

Round2Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2015_Estimates.txt")
Round2Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2015_EstimatesStats.txt")


# Verify that Gelman-Rubin < 1.2
sapply(Round2Mixtures_2015_Estimates$Stats, function(Mix) {Mix[, "GR"]})
sapply(Round2Mixtures_2015_Estimates$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUGANC14_3_Late
require(gplots)
par(mfrow = c(1, 1), mar = c(3.1, 4.6, 3.1, 1.1), oma = rep(0, 4))
sapply(Round2Mixtures_2015[c(4, 5, 3, 2, 1)], function(Mix) {
  BarPlot <- barplot2(Round2Mixtures_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '', cex.lab = 1.5, cex.main = 2)
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
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
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round2Mixtures_2015_Header, set.mfrow = c(5, 3), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 2015 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round2Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round2Mixtures_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round2Mixtures_2015_Header)

## Quick barplot
QuickBarplot(mixvec = Round2Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round2Mixtures_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round2Mixtures_2015_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round2Mixtures_2015[c(4, 5, 3, 2, 1)], estimates = Round2Mixtures_2015_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round2Mixtures_2015_Header)
rm(Round2Mixtures_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 2 2015 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures_2015_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
                                                                   mixvec = Round2Mixtures_2015, prior = "",  
                                                                   ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round2Mixtures_2015_31RG_Estimates)
QuickBarplot(mixvec = Round2Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round2Mixtures_2015_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round2Mixtures_2015_Header)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 2 2015 Repeated MSA files for BAYES 80K Iterations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping out a new control file to re-analyze SUGANC15_2_Late Round 2 mixtures with a 80K posterior as opposed to our standard 40K
## Perhaps it will converge on one of the modes for South of Cape Suckling?

# Just use original mixture, no need to create a new one!

## Dumping a New Control file that goes out further!

sapply(Round2Mixtures_2015[5], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round2Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round2Mixtures_2015[5], function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})
# Will rename later, but leave as is for now for BayesCopyPaste.GCL


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 2 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round2Mixtures80K_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
  mixvec = Round2Mixtures_2015[5], prior = "",  
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
sapply(Round2Mixtures_2015[5], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round2Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round2Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round2Mixtures80K_2015_Estimates$Output)

PlotPosterior(mixvec = Round2Mixtures_2015[5], output = Round2Mixtures80K_2015_Estimates$Output, 
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round2Mixtures_2015_Header[5], set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 2 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round2Mixtures_2015[5], estimatesstats = Round2Mixtures80K_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round2Mixtures_2015_Header[3])

## Quick barplot
QuickBarplot(mixvec = Round2Mixtures_2015[5], estimatesstats = Round2Mixtures80K_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round2Mixtures_2015_Header[5])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round2Mixtures_2015[5], estimates = Round2Mixtures80K_2015_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round2Mixtures_2015_Header)
rm(Round2Mixtures80K_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round2Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round2Mixtures_2015_EstimatesStats.txt")

sapply(Round2Mixtures_2015[5], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round2Mixtures_2015_EstimatesStats[[Mix]][, "median"],
                                                 Round2Mixtures80K_2015_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round2Mixtures_2015_EstimatesStats[[Mix]][, "5%"],
                                               Round2Mixtures80K_2015_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round2Mixtures_2015_EstimatesStats[[Mix]][, "95%"],
                                               Round2Mixtures80K_2015_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round2Mixtures_2015_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA15GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 2 2015 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SALITC15_2_Middle, SAYAKC15_2_Middle, SKARLC15_2_Middle, and SUGANC15_2_Middle as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from  and SUYAKC15_2_Middle
Round2Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round2Mixtures_2015_Estimates.txt")

sapply(Round2Mixtures_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SUYAKC15_2_Middle
sapply(Round2Mixtures80K_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # SUYAKC15_2_Middle resolved

str(Round2Mixtures_2015_Estimates)
str(Round2Mixtures80K_2015_Estimates)

Round2Mixtures_2015_Estimates_Final <- list(Stats = c(Round2Mixtures_2015_Estimates$Stats[Round2Mixtures_2015[1:4]],
                                                      Round2Mixtures80K_2015_Estimates$Stats[Round2Mixtures_2015[5]]),
                                            Output = c(Round2Mixtures_2015_Estimates$Output[Round2Mixtures_2015[1:4]],
                                                       Round2Mixtures80K_2015_Estimates$Output[Round2Mixtures_2015[5]]))
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
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
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
  Prior.GCL(groupvec = KMA473PopsGroupVec15, groupweights = Mix[, "mean"], minval = 0.01)}, simplify = FALSE)
names(Round3Mixtures_2015_Prior) <- gsub(pattern = "2_Middle", replacement = "3_Late", 
                                         x = names(Round3Mixtures_2015_Prior))  # This changes the names
dput(x = Round3Mixtures_2015_Prior, file = "Objects/Round3Mixtures_2015_Prior.txt")
str(Round3Mixtures_2015_Prior)


## Dumping Mixture files
sapply(Round3Mixtures_2015, function(Mix) {CreateMixture.GCL(sillys = Mix, loci = loci46, IDs = NULL, mixname = Mix, dir = "BAYES/2014-2015 Mixtures 46loci/Mixture", type = "BAYES", PT = FALSE)} )

## Dumping Control files
sapply(Round3Mixtures_2015, function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 40000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round3Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round3Mixtures_2015, function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Go run BAYES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 3 2015 Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC ,
  maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
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
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
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
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round3Mixtures_2015_Header, set.mfrow = c(5, 3), thin = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 3 2015 Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round3Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round3Mixtures_2015_Header)

## Quick barplot
QuickBarplot(mixvec = Round3Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round3Mixtures_2015_Header)

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round3Mixtures_2015[c(4, 5, 3, 2, 1)], estimates = Round3Mixtures_2015_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round3Mixtures_2015_Header)
rm(Round3Mixtures_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 3 2015 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2015_31RG_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
  mixvec = Round3Mixtures_2015, prior = "",  
  ext = "BOT", nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
str(Round3Mixtures_2015_31RG_Estimates)
QuickBarplot(mixvec = Round3Mixtures_2015[c(4, 5, 3, 2, 1)], estimatesstats = Round3Mixtures_2015_31RG_Estimates$Stats, groups = KMA31GroupsPC, header = Round3Mixtures_2015_Header)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Round 3 2015 Repeated MSA files for BAYES 80K Iterations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dumping out a new control file to re-analyze SUGANC15_2_Late Round 3 mixtures with a 80K posterior as opposed to our standard 40K
## Perhaps it will converge on one of the modes for South of Cape Suckling?

# Just use original mixture, no need to create a new one!

## Dumping a New Control file that goes out further!

sapply(Round3Mixtures_2015[c(1, 4)], function(Mix) {
  CreateControlFile.GCL(sillyvec = KMA473Pops, loci = loci46, mixname = Mix, basename = "KMA473Pops46Markers", suffix = "", nreps = 80000, nchains = 5,
                        groupvec = KMA473PopsGroupVec15, priorvec = Round3Mixtures_2015_Prior[[Mix]], initmat = KMA473PopsInits, dir = "BAYES/2014-2015 Mixtures 46loci/Control",
                        seeds = WASSIPSockeyeSeeds, thin = c(1, 1, 100), mixfortran = KMA47346MixtureFormat, basefortran = KMA47346Baseline, switches = "F T F T T T F")
})

## Create output directory
sapply(Round3Mixtures_2015[c(1, 4)], function(Mix) {dir.create(paste("BAYES/2014-2015 Mixtures 46loci/Output/", Mix, sep = ""))})
# Will rename later, but leave as is for now for BayesCopyPaste.GCL


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Round 3 Repeated 80K Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures80K_2015_Estimates <- CustomCombineBAYESOutput.GCL(
  groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, 
  maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
  mixvec = Round3Mixtures_2015[c(1, 4)], prior = "",  
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
sapply(Round3Mixtures_2015[c(1, 4)], function(Mix) {
  par(mfrow = c(1, 1))
  BarPlot <- barplot2(Round3Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"], col = "blue", ylim = c(1, pmax(1.5, max(Round3Mixtures80K_2015_EstimatesStats[[Mix]][, "GR"]))), ylab = "Gelman-Rubin", xpd = FALSE, main = Mix, names.arg = '')
  abline(h = 1.2, lwd = 3, xpd = FALSE)
  text(x = BarPlot, y = 1, labels = KMA15GroupsPC2Rows, srt = 0, pos = 1, xpd = TRUE, cex = 0.55)
})

# Quick look at raw posterior output
str(Round3Mixtures80K_2015_Estimates$Output)

PlotPosterior(mixvec = Round3Mixtures_2015[c(1, 4)], output = Round3Mixtures80K_2015_Estimates$Output, 
              groups = KMA15GroupsPC, colors = KMA15Colors, 
              header = Round3Mixtures_2015_Header[c(1, 4)], set.mfrow = c(5, 3), thin = 10)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Round 3 Repeated 80K Results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quick and dirty plot
QuickPlot(mixvec = Round3Mixtures_2015[c(1, 4)], estimatesstats = Round3Mixtures80K_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round3Mixtures_2015_Header[3])

## Quick barplot
QuickBarplot(mixvec = Round3Mixtures_2015[c(1, 4)], estimatesstats = Round3Mixtures80K_2015_Estimates, groups = KMA15GroupsPC, groups2rows = KMA15GroupsPC2Rows, header = Round3Mixtures_2015_Header[c(1, 4)])

## Make violin plots of posteriors with RGs sorted
ViolinPlot(mixvec = Round3Mixtures_2015[c(1, 4)], estimates = Round3Mixtures80K_2015_Estimates, groups = KMA15GroupsPC2Rows, colors = KMA15Colors, header = Round3Mixtures_2015_Header)
rm(Round3Mixtures80K_2015_Estimates)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compare Results from 40K to 80K
Round3Mixtures_2015_EstimatesStats <- dget(file = "Estimates objects/Round3Mixtures_2015_EstimatesStats.txt")

sapply(Round3Mixtures_2015[c(1, 4)], function(Mix) {
  par(mar = c(3.1, 5.1, 3.1, 2.1))
  Barplot <- barplot2(height = t(matrix(data = c(Round3Mixtures_2015_EstimatesStats[[Mix]][, "median"],
                                                 Round3Mixtures80K_2015_EstimatesStats[[Mix]][, "median"]) * 100, 
                                        ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      beside = TRUE, plot.ci = TRUE, ci.lwd = 2, 
                      ci.l = t(matrix(data = c(Round3Mixtures_2015_EstimatesStats[[Mix]][, "5%"],
                                               Round3Mixtures80K_2015_EstimatesStats[[Mix]][, "5%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ci.u = t(matrix(data = c(Round3Mixtures_2015_EstimatesStats[[Mix]][, "95%"],
                                               Round3Mixtures80K_2015_EstimatesStats[[Mix]][, "95%"]) * 100,
                                      ncol = 2, dimnames = list(KMA15GroupsPC, c("40K", "80K")))),
                      ylim = c(0, 100), col = c("white", "skyblue"), yaxt = 'n', xaxt = 'n', 
                      main = Round3Mixtures_2015_Header[Mix], ylab = "Percent of Mixture", cex.lab = 2, cex.main = 2)
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  mtext(text = KMA15GroupsPC2Rows, side = 1, line = 1, at = apply(Barplot, 2, mean), adj = 0.5, cex = 0.6)
  legend("topleft", legend = c("40K Iterations", "80K Iterations"), fill = c("White", "skyblue"), bty = 'n', cex = 2)
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dputting Final Round 3 2015 Estimates ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Keeping 40K iteration results from SAYAKC15_3_Late, SKARLC15_3_Late, and SUYAKC15_3_Late as GR was < 1.2 and thats "our business rule"
# Using 80K iteration results from  and SALITC15_3_Late and SUGANC15_3_Late
Round3Mixtures_2015_Estimates <- dget(file = "Estimates objects/Round3Mixtures_2015_Estimates.txt")

sapply(Round3Mixtures_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Issues with SALITC15_3_Late & SUGANC15_3_Late
sapply(Round3Mixtures80K_2015_EstimatesStats, function(Mix) {table(Mix[, "GR"] > 1.2)})  # Both resolved

str(Round3Mixtures_2015_Estimates)
str(Round3Mixtures80K_2015_Estimates)

Round3Mixtures_2015_Estimates_Final <- list(Stats = c(Round3Mixtures_2015_Estimates$Stats[Round3Mixtures_2015[c(2:3, 5)]],
                                                      Round3Mixtures80K_2015_Estimates$Stats[Round3Mixtures_2015[c(1, 4)]]),
                                            Output = c(Round3Mixtures_2015_Estimates$Output[Round3Mixtures_2015[c(2:3, 5)]],
                                                       Round3Mixtures80K_2015_Estimates$Output[Round3Mixtures_2015[c(1, 4)]]))
str(Round3Mixtures_2015_Estimates_Final)
dput(x = Round3Mixtures_2015_Estimates_Final, file = "Estimates objects/Round3Mixtures_2015_Estimates_Final.txt")
dput(x = Round3Mixtures_2015_Estimates_Final$Stats, file = "Estimates objects/Round3Mixtures_2015_EstimatesStats_Final.txt")

Round3Mixtures_2015_Estimates_Final <- dget(file = "Estimates objects/Round3Mixtures_2015_Estimates_Final.txt")
Round3Mixtures_2015_EstimatesStats_Final <- dget(file = "Estimates objects/Round3Mixtures_2015_EstimatesStats_Final.txt")

sapply(Round3Mixtures_2015_Estimates_Final$Stats, function(Mix) {table(Mix[, "GR"] > 1.2)})

# Dput final Middle Strata Estimates from 2015
dput(x = Round3Mixtures_2015_EstimatesStats_Final, file = "Estimates objects/Final/KMA2015Strata_3_Late_EstimatesStats.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Sneak Peak Round 3 Results with Finer Scale Reporting Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Round3Mixtures_2015_31RG_Estimates <- CustomCombineBAYESOutput.GCL(groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, 
                                                                   maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
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

ProportionColors <- colorpanel(n = 3, low = "blue", high = "white")
TempProportionColors15 <- sapply(KMA2015, function(geo) {
  ProportionColors[sapply(TempMix15[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
}, simplify = FALSE)

Estimates <- KMA2015Strata_EstimatesStats
Groups <- KMA15GroupsPC
Groups2Rows <- KMA15GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5

# dir.create("Figures/2015")

sapply(names(TempMix15[c(4, 5, 3, 2, 1)]), function(geomix) {
  emf(file = paste("Figures/2015/", geomix, ".emf", sep = ''), width = 8.5, height = 6.5, family = "sans", bg = "white")
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
Groups <- KMA15GroupsPC
Groups2Rows <- KMA15GroupsPC2Rows
cex.lab <- 1.5
cex.yaxis <- 1.3
cex.xaxis <- 0.6
cex.main <- 1.7
cex.leg <- 1.3
ci.lwd <- 2.5
ymax <- 200000

sapply(names(TempMix15[c(4, 5, 3, 2, 1)]), function(geomix) {
  emf(file = paste("Figures/2015/", geomix, "Harvest.emf", sep = ''), width = 8.5, height = 6.5, family = "sans", bg = "white")
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

# For this to work properly, the output mixture directories all need to be in the same place (can be moved back afterwards)
sapply(KMA2015, function(geomix) {
  assign(x = paste(geomix, "_Annual_Stratified", sep = ''), 
         value = StratifiedEstimator.GCL(
           groupvec = seq(KMA15GroupsPC), groupnames = KMA15GroupsPC, maindir = "BAYES/2014-2015 Mixtures 46loci/Output", 
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
  Stats <- get(paste(geomix, "_Annual_Stratified", sep = ''))$Summary[, c("MEAN", "SD", "MEDIAN", "5%", "95%", "P0", "GR")]
  colnames(Stats) <- c("mean", "sd", "median", "5%", "95%", "P=0", "GR")
  Stats
}, simplify = FALSE)
str(KMA2015_Annual_EstimatesStats)
dput(x = KMA2015_Annual_EstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_EstimatesStats.txt")




# Create a list object with all Stratified Annual Harvest
KMA2015_Annual_HarvestEstimatesStats <- sapply(KMA2015, function(strata) {
  cbind(KMA2015_Annual_EstimatesStats[[strata]][, c("mean", "sd", "median", "5%", "95%")] * sum(HarvestByStrata2015_Final[strata, ], na.rm = TRUE),
        KMA2015_Annual_EstimatesStats[[strata]][, c("P=0", "GR")])
}, simplify = FALSE )

dput(x = KMA2015_Annual_HarvestEstimatesStats, file = "Estimates objects/Final/KMA2015_Annual_HarvestEstimatesStats.txt")
str(KMA2015_Annual_HarvestEstimatesStats)


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

Annual2015_Stratified_Estimates <- sapply(KMA2015, function(geomix) {
  KMA2015_Annual_EstimatesStats[[geomix]][, "mean"]
})


require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2015_Stratified_Estimates[15:1, c(4,5,3,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2015 Proportion", at = seq(0, 1, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2014_Stratified_Estimates[15:1, c(4,5,3,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
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
levelplot(t(Annual2015_Stratified_HarvestEstimates[15:1, c(4,5,3,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2015 Harvest", at = seq(0, 250000, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

new.colors <- colorRampPalette(c("white", "black"))
levelplot(t(Annual2014_Stratified_HarvestEstimates[15:1, c(4,5,3,2,1)]), col.regions = new.colors, xlab = "Reporting Group", 
          ylab = "Fishery", main = "2014 Harvest", at = seq(0, 250000, length.out = 100), aspect = "fill", 
          scales = list(x = list(rot = 45)),
          panel = function(...) {
            panel.levelplot(...)
          })

sort(rowSums(Annual2015_Stratified_HarvestEstimates))
sort(rowSums(Annual2014_Stratified_HarvestEstimates))



# Bubble chart example
apply(col2rgb(KMA15Colors), 2, function(col) {rgb(red = col[1], green = col[2], blue = col[3], maxColorValue = 255)} )

# 2015
Annual2015_Stratified_HarvestEstimates_df <- melt(Annual2015_Stratified_HarvestEstimates)
names(Annual2015_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2015_Stratified_HarvestEstimates_df$RG <- factor(Annual2015_Stratified_HarvestEstimates_df$RG, levels = rev(KMA15GroupsPC))
Annual2015_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2015_Stratified_HarvestEstimates_df$Fishery, levels = KMA2015[c(4,5,3,2,1)])
Annual2015_Stratified_HarvestEstimates_df$Color <- rep(rev(KMA15Colors), 5)
str(Annual2015_Stratified_HarvestEstimates_df)


ggplot(data = Annual2015_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + geom_point() + scale_size_area(max_size = 25) + scale_color_manual(values = rev(rep(KMA15Colors, 5)))


ggplot(data = Annual2015_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, 250000), breaks = seq(50000, 250000, 50000), range = c(0, 25)) + 
  scale_color_manual(values = rev(rep(KMA15Colors, 5))) +
  ggtitle("2015 Harvest")


# 2014
Annual2014_Stratified_HarvestEstimates_df <- melt(Annual2014_Stratified_HarvestEstimates)
names(Annual2014_Stratified_HarvestEstimates_df) <- c("RG", "Fishery", "Harvest")
Annual2014_Stratified_HarvestEstimates_df$RG <- factor(Annual2014_Stratified_HarvestEstimates_df$RG, levels = rev(KMA15GroupsPC))
Annual2014_Stratified_HarvestEstimates_df$Fishery <- factor(Annual2014_Stratified_HarvestEstimates_df$Fishery, levels = KMA2014[c(4,5,3,2,1)])
Annual2014_Stratified_HarvestEstimates_df$Color <- rep(rev(KMA15Colors), 5)
str(Annual2014_Stratified_HarvestEstimates_df)

ggplot(data = Annual2014_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + geom_point() + scale_size_area(max_size = 25) + scale_color_manual(values = rev(rep(KMA15Colors, 5)))

ggplot(data = Annual2014_Stratified_HarvestEstimates_df, aes(x = Fishery, y = RG, size = Harvest, color = RG)) + 
  geom_point() + 
  scale_size_continuous(limits = c(0, 250000), breaks = seq(50000, 250000, 50000), range = c(0, 25)) + 
  scale_color_manual(values = rev(rep(KMA15Colors, 5))) +
  ggtitle("2014 Harvest")


# Create a matrix of early strata means
KMA2015Strata_1_Early_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_1_Early_EstimatesStats.txt")
EarlyStrata2015_Estimates <- sapply(KMA2015Strata_1_Early_EstimatesStats, function(geomix) {
  geomix[, "mean"]
})
colnames(EarlyStrata2015_Estimates) <- KMA2015

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Annual vs. Early Means
par(mar = c(4.1, 5.1, 3.1, 1.1))
sapply(KMA2015, function(geomix) {
  Barplot <- barplot2(height = rbind(Annual2015_Stratified_Estimates[, geomix], EarlyStrata2015_Estimates[, geomix]) * 100, 
                      beside = TRUE, col = c("blue", "white"), yaxt = 'n', xaxt = 'n', main = geomix, ylab = "Precent of Mixture",
                      cex.lab = 2, cex.main = 2, ylim = c(0, 100))
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = ",", digits = 0, format = "f"), cex.axis = 1.5)
  abline(h = 0, xpd = FALSE)
  text(x = colMeans(Barplot), y = -1, labels = KMA15GroupsPC2Rows, srt = 90, adj = 1, xpd = TRUE, cex = 0.6)
  legend("topleft", legend = c("Annual Means", "Early Strata Means"), fill = c("blue", "white"), bty = 'n', cex = 1.5)
})




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Percentages for KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_EstimatesStats.txt")
KMA2015Strata_EstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_EstimatesStats.txt")

str(KMA2014Strata_EstimatesStats)




# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             1, 3,
                             1, 4,
                             5, 6), nrow = 4, ncol = 2, byrow = TRUE)



GeoMix <- c("SUGANC", "SUYAKC", "SKARLC", "SAYAKC", "SALITC")

filenames <- setNames(object = c("Uganik Proportions 2014-2016", 
                                 "Uyak Proportions 2014-2016",
                                 "Karluk Proportions 2014-2016",
                                 "Ayakulik Proportions 2014-2016",
                                 "Alitak Proportions 2014-2016"), nm = GeoMix)

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
# TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)
# 
# Legend16 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
#                      nm = c("1_Early", "2_Middle", "3_Late"))
# TempLegend16 <- sapply(KMA2016, function(geo) {
#   Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
# }, simplify = FALSE)
# 
# 
# TempProportionColors16 <- sapply(KMA2016, function(geo) {
#   ProportionColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
# }, simplify = FALSE)
# 
# Estimates16 <- KMA2016Strata_EstimatesStats

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA15GroupsPC
Groups2Rows <- KMA15GroupsPC2Rows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.leg <- 1.3
ci.lwd <- 2.5


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)

sapply(GeoMix, function(geomix) {
  
  emf(file = paste("Figures/All Years/", filenames[geomix], ".emf", sep = ''), width = 7, height = 7, family = "sans", bg = "white")
  
  
  layout(mat = layoutmat, widths = c(0.075, 1, 1), heights = c(0.9, 0.9, 1, 0.1))
  par(mar = rep(0, 4))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Y-axis label
  plot.new()
  text(x = 0.25, y = 0.5, labels = "Percentage of Catch", srt = 90, cex = cex.lab)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2014 Barplot
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
  
  par(mar = c(2, 1, 1, 1))
  Barplot16 <- barplot2(height = matrix(data = rep(0, 45), nrow = 3, ncol = length(KMA15GroupsPC)),
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = matrix(data = rep(0, 45), nrow = 3, ncol = length(KMA15GroupsPC)),
                        ci.u = matrix(data = rep(0, 45), nrow = 3, ncol = length(KMA15GroupsPC)),
                        ylim = c(0, 100), col = "black", yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14[[geomix14]], x = "topleft", fill = TempProportionColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2016")
  abline(h = 0, xpd = FALSE)
  
  # geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
  # par(mar = c(2, 1, 1, 1))
  # Barplot16 <- barplot2(height = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][, "median"]})) * 100, 
  #                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
  #                       ci.l = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][, "5%"]})) * 100, 
  #                       ci.u = t(sapply(TempMix16[[geomix16]], function(tempmix) {Estimates16[[tempmix]][, "95%"]})) * 100, 
  #                       ylim = c(0, 100), col = TempProportionColors16[[geomix16]], cex.axis = cex.yaxis, yaxt = "n", xaxt = 'n')
  # axis(side = 2, at = seq(0, 100, 25), labels = formatC(x = seq(0, 100, 25), big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  # legend(legend = TempLegend16[[geomix16]], x = "topleft", fill = TempProportionColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")
  # abline(h = 0, xpd = FALSE)
  # 
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot16, 2, mean), adj = 0.5, cex = cex.xaxis)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Blank Corner
  par(mar = rep(0, 4))
  plot.new()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## x-axis label
  par(mar = rep(0, 4))
  plot.new()
  text(x = 0.5, y = 0.5, labels = "Reporting Group", cex = cex.lab)
  
  
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Plot Harvests for KMA Mixtures ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA2014Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2014Strata_HarvestEstimatesStats.txt")
KMA2015Strata_HarvestEstimatesStats <- dget(file = "Estimates objects/Final/KMA2015Strata_HarvestEstimatesStats.txt")

str(KMA2014Strata_EstimatesStats)




# Three barplot layout
layoutmat <- matrix(data=c(  1, 2,
                             1, 3,
                             1, 4,
                             5, 6), nrow = 4, ncol = 2, byrow = TRUE)



GeoMix <- c("SUGANC", "SUYAKC", "SKARLC", "SAYAKC", "SALITC")

filenames <- setNames(object = c("Uganik Harvest 2014-2016", 
                                 "Uyak Harvest 2014-2016",
                                 "Karluk Harvest 2014-2016",
                                 "Ayakulik Harvest 2014-2016",
                                 "Alitak Harvest 2014-2016"), nm = GeoMix)

# If showing proportions (percetages) use blue, otherwise green as "low"
HarvestColors <- colorpanel(n = 3, low = "green", high = "white")


#~~~~~~~~~~~~~~~~~~
# 2014
TempMix14 <- sapply(KMA2014, function(geo) {grep(pattern = geo, x = names(KMA2014Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

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
TempMix15 <- sapply(KMA2015, function(geo) {grep(pattern = geo, x = names(KMA2015Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)

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
# TempMix16 <- sapply(KMA2016, function(geo) {grep(pattern = geo, x = names(KMA2016Strata_EstimatesStats), value = TRUE)}, simplify = FALSE)
# 
# Legend16 <- setNames(object = c("June 1-27", "June 28-July 25", "July 26-August 29"), 
#                      nm = c("1_Early", "2_Middle", "3_Late"))
# TempLegend16 <- sapply(KMA2016, function(geo) {
#   Legend16[sapply(TempMix16[[geo]], function(strata) {unlist(strsplit(x = strata, split = paste(geo, "_", sep = '')))[2]} )]
# }, simplify = FALSE)
# 
# 
# TempHarvestColors16 <- sapply(KMA2016, function(geo) {
#   HarvestColors[sapply(TempMix16[[geo]], function(strata) {as.numeric(unlist(strsplit(x = strata, split = "_"))[2])} )]
# }, simplify = FALSE)
# 
# HarvestEstimates16 <- KMA2016Strata_HarvestEstimatesStats

#~~~~~~~~~~~~~~~~~~
# Size Parameters
Groups <- KMA15GroupsPC
Groups2Rows <- KMA15GroupsPC2Rows
cex.lab <- 1.5
cex.xaxis <- 0.5
cex.leg <- 1.3
ci.lwd <- 2.5
ymax <- max(sapply(c(HarvestEstimates14, HarvestEstimates15), function(strata) {strata[, "95%"]}))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make figures as .emf files

# dir.create("Figures/All Years")
require(devEMF)
require(gplots)

sapply(GeoMix, function(geomix) {
  
  emf(file = paste("Figures/All Years/", filenames[geomix], ".emf", sep = ''), width = 7, height = 7, family = "sans", bg = "white")
  
  
  layout(mat = layoutmat, widths = c(0.075, 1, 1), heights = c(0.9, 0.9, 1, 0.1))
  par(mar = rep(0, 4))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Y-axis label
  plot.new()
  text(x = 0.25, y = 0.5, labels = "Number of Fish Harvested (Thousands)", srt = 90, cex = cex.lab)
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2014 Barplot
  geomix14 <- grep(pattern = geomix, x = names(TempMix14), value = TRUE)
  par(mar = c(1, 1, 1, 1))
  Barplot14 <- barplot2(height = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][, "median"]})), 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][, "5%"]})), 
                        ci.u = t(sapply(TempMix14[[geomix14]], function(tempmix) {HarvestEstimates14[[tempmix]][, "95%"]})), 
                        ylim = c(0, ymax), col = TempHarvestColors14[[geomix14]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 40000), labels = formatC(x = seq(0, ymax, 40000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14[[geomix14]], x = "topleft", fill = TempHarvestColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2014")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot14, 2, mean), adj = 0.5, cex = cex.xaxis)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2015 Barplot
  geomix15 <- grep(pattern = geomix, x = names(TempMix15), value = TRUE)
  par(mar = c(1, 1, 1, 1))
  Barplot15 <- barplot2(height = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][, "median"]})), 
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][, "5%"]})), 
                        ci.u = t(sapply(TempMix15[[geomix15]], function(tempmix) {HarvestEstimates15[[tempmix]][, "95%"]})), 
                        ylim = c(0, ymax), col = TempHarvestColors15[[geomix15]], yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 40000), labels = formatC(x = seq(0, ymax, 40000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend15[[geomix15]], x = "topleft", fill = TempHarvestColors15[[geomix15]], border = "black", bty = "n", cex = cex.leg, title="2015")
  abline(h = 0, xpd = FALSE)
  
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot15, 2, mean), adj = 0.5, cex = cex.xaxis)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## 2016 Barplot
  
  par(mar = c(2, 1, 1, 1))
  Barplot16 <- barplot2(height = matrix(data = rep(0, 45), nrow = 3, ncol = length(KMA15GroupsPC)),
                        beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
                        ci.l = matrix(data = rep(0, 45), nrow = 3, ncol = length(KMA15GroupsPC)),
                        ci.u = matrix(data = rep(0, 45), nrow = 3, ncol = length(KMA15GroupsPC)),
                        ylim = c(0, ymax), col = "black", yaxt = "n", xaxt = 'n')
  axis(side = 2, at = seq(0, ymax, 40000), labels = formatC(x = seq(0, ymax, 40000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  legend(legend = TempLegend14[[geomix14]], x = "topleft", fill = TempHarvestColors14[[geomix14]], border = "black", bty = "n", cex = cex.leg, title="2016")

  
  # geomix16 <- grep(pattern = geomix, x = names(TempMix16), value = TRUE)
  # par(mar = c(2, 1, 1, 1))
  # Barplot16 <- barplot2(height = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][, "median"]})), 
  #                       beside = TRUE, plot.ci = TRUE, ci.lwd = ci.lwd,
  #                       ci.l = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][, "5%"]})), 
  #                       ci.u = t(sapply(TempMix16[[geomix16]], function(tempmix) {HarvestEstimates16[[tempmix]][, "95%"]})), 
  #                       ylim = c(0, ymax), col = TempHarvestColors16[[geomix16]], yaxt = "n", xaxt = 'n')
  # axis(side = 2, at = seq(0, ymax, 40000), labels = formatC(x = seq(0, ymax, 40000) / 1000, big.mark = "," , digits = 0, format = "f"), cex.axis = cex.yaxis)
  # legend(legend = TempLegend16[[geomix16]], x = "topleft", fill = TempHarvestColors16[[geomix16]], border = "black", bty = "n", cex = cex.leg, title="2016")
  # abline(h = 0, xpd = FALSE)
  # 
  mtext(text = Groups2Rows, side = 1, line = 1, at = apply(Barplot16, 2, mean), adj = 0.5, cex = cex.xaxis)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Blank Corner
  par(mar = rep(0, 4))
  plot.new()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## x-axis label
  par(mar = rep(0, 4))
  plot.new()
  text(x = 0.5, y = 0.5, labels = "Reporting Group", cex = cex.lab)
  
  
  dev.off()
})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table results ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Dates
DatesStrata2014 <- read.table(file = "Harvest/2014DatesByStrata.txt", header = TRUE, sep = "\t", as.is = TRUE)
DatesStrata2014.mat <- as.matrix(DatesStrata2014[-1])
dimnames(DatesStrata2014.mat) <- list(DatesStrata2014$location, c("1_Early", "2_Middle", "3_Late"))
dput(x = DatesStrata2014.mat, file = "Objects/DatesStrata2014_Final.txt"); rm(DatesStrata2014.mat)
DatesStrata2014_Final <- dget(file = "Objects/DatesStrata2014_Final.txt")

# Sample sizes
KMA2014_2015Strata_SampleSizes_Final <- KMA2014_2015Strata_SampleSizes[, "Final"]



# Get Final Estimates Objects
KMAfinalestimatesobjects <- list.files(path = "Estimates objects/Final", recursive = FALSE)
invisible(sapply(KMAfinalestimatesobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste("Estimates objects/Final", objct, sep = "/")), pos = 1) })); beep(2)
KMAfinalestimatesobjects; rm(KMAfinalestimatesobjects)

SheetNames <- sort(c(names(KMA2014Strata_EstimatesStats), names(KMA2014_Annual_EstimatesStats)))
names(SheetNames) <- SheetNames

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Defining caption variables

mixvec <- SheetNames



mix <- SheetNames[2]
harvest <- HarvestByStrata2014_Final
dates <- DatesStrata2014_Final
sampsize <- KMA2014_2015Strata_SampleSizes[, "Final"]



# String split by "_
SheetNames.split <- unlist(strsplit(x = mix, split = "_"))

# The first element is the geographic area + year
geomix <- SheetNames.split[1]

# If it is not an annual roll-up, then get the strata number + strata name
if(length(SheetNames.split) > 1) {
  tempmix <- paste(c(SheetNames.split[2], SheetNames.split[3]), collapse = "_")
  
  Stratum <- paste("stratum ", SheetNames.split[2], " (", dates[geomix, tempmix], "; Harvest=", formatC(x = harvest[geomix, tempmix], format = "f", digits = 0, big.mark = ","), "; n=", sampsize[mix], ")", sep = '')
} else {
    Stratum <-
  }







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
