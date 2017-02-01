## KMA Sockeye baseline
## Kyle Shedd Mon Oct 19 11:00:37 2015
date()

## The goal of this script is to create a coastwide sockeye baseline for the KMA sockecye MSA project that contains the following 14 RGs
## 1) West of Chignik 2) Black Lake 3) Chignik Lake 4-12) 9 intra-Kodiak RGs 13) Cook Inlet 14) PWS / Copper River and 15) SEAK / BC / WA


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Directory Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Create file structure based on OLD
# setwd("V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline")
# dirsOLD <- list.dirs(path = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/OLD", recursive = FALSE, full.names = FALSE)
# sapply(dirsOLD, dir.create)
# sapply(c("Baseline", "BAYES.baseline", "BAYES.control", "BAYES.mixture", "BAYES.output"), function(folder) {dir.create(path = paste(getwd(), "BAYES", folder, sep = "/"))})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline")

# This sources all of the new GCL functions to this workspace
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")
# This sources all of Kyle's non-GCL useful functions to this workspace (does not include "tempGCL" functions)
source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")

# username <- "krshedd"
# password <- "Ski10er!"

## save.image("V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/KMASockeyeBaseline.RData")
## load("V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/KMASockeyeBaseline.RData")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Get .gcl objects out of LOKI ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# KMA761CollectionsLOKI2R <- readClipboard()
# dput(x = KMA761CollectionsLOKI2R, file = "Objects/KMA761CollectionsLOKI2R.txt")
KMA761CollectionsLOKI2R <- dget(file = "Objects/KMA761CollectionsLOKI2R.txt")

# require(beepr)
# source("V:/WORK/Kyle/R Scripts/ReadLOKIv3.GCL.R")
# ReadLOKIv3.GCL(sillyvec = KMA745CollectionsReadLOKI, markersuite = "Sockeye2011_96SNPs", onlymarkersuitefish = TRUE) # Kyle's update to Jim Jasper's ReadLOKIv2.GCL to will only grab fish from the specified markersuite

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
LOKI2R.GCL(sillyvec = KMA761CollectionsLOKI2R, username = username, password = password)
objects(pattern = "\\.gcl")


## Save unaltered .gcl's as back-up:
dir.create("Raw genotypes/OriginalCollections")
invisible(sapply(KMA761CollectionsLOKI2R, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollections/" , silly, ".txt", sep = ''))} )); beep(8)

collection.size.original <- sapply(KMA761CollectionsLOKI2R, function(silly) get(paste(silly, ".gcl", sep = ""))$n)


#~~~~~~~~~~~~~~~~~~
## Remove NA Individuals
# Individuals not genotyped for the full suite of 96 SNPs
NA.indv.removed <- RemoveNAindv.GCL(sillyvec = KMA761CollectionsLOKI2R)
dput(x = NA.indv.removed, file = "Objects/NA.indv.removed.txt")

collection.size.postNA <- sapply(KMA761CollectionsLOKI2R, function(silly) get(paste(silly, ".gcl", sep = ""))$n)

# How many fish removed from each SILLY
hist(collection.size.original - collection.size.postNA, col = 8)


#~~~~~~~~~~~~~~~~~~
## Look for "weird" scores (i.e. "Unk" or "XXX")
table(sapply(KMA761CollectionsLOKI2R, function(silly) {length(grep(pattern = "Unk", x = get(paste(silly, ".gcl", sep = ""))$scores, , "Dose1"))} ))
table(sapply(KMA761CollectionsLOKI2R, function(silly) {length(grep(pattern = "XXX", x = get(paste(silly, ".gcl", sep = ""))$scores, , "Dose1"))} ))
# None found


#~~~~~~~~~~~~~~~~~~
## Make 0 -> NA
# NA fish (never genotyped for a locus) have been removed
# All fish remaining in .gcls have been genotyped for each locus in markersuite
# There are no "weird" scores
# Thus, to maintain compatability with older .GCL functions, I will push all 0
# scores (genotyped, but no-called) into NA.

invisible(sapply(KMA761CollectionsLOKI2R, function(silly) {
  my.gcl <- get(paste(silly, ".gcl", sep = ""))
  my.gcl$scores <- gsub(pattern = "0", replacement = NA, x = my.gcl$scores)
  my.gcl$scores <- gsub(pattern = "Unk", replacement = NA, x = my.gcl$scores)
  my.gcl$scores <- gsub(pattern = "XXX", replacement = NA, x = my.gcl$scores)
  assign(x = paste(silly, ".gcl", sep = ""), value = my.gcl, pos = 1)
}))


## Save clean .gcl's as back-up:
dir.create("Raw genotypes/OriginalCollectionsNAremoved")
invisible(sapply(KMA761CollectionsLOKI2R, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/OriginalCollectionsNAremoved/" , silly, ".txt", sep = ''))} )); beep(8)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pooling SAYAK11, SAYAK12, SAUKE13 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Need to pool collections for SAYAK11L and SAYAK12E and SAYAK12L
## Need to split SAUKE13
# Script comes from Kodiak 2014 baseline and SEAK 2015 baseline

InitialPooling <- c("SAYAK11", "SAYAK12", "SAUKE13")

# LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
# loci96 <- dget(file = "Objects/loci96.txt")
# 
# sapply(InitialPooling, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/", silly, ".txt", sep = "")), pos = 1)}); beep(8)
# objects(pattern = "\\.gcl")


## Ayakulik Weir samples (split into Early and Late run)
# First get late run fish from SAYAK11.gcl
table(SAYAK11.gcl$attributes$CAPTURE_DATE)
unique(SAYAK11.gcl$attributes$CAPTURE_DATE)[5]
SAYAK11LIDs <- AttributesToIDs.GCL(silly = "SAYAK11", attribute = "CAPTURE_DATE", matching = unique(SAYAK11.gcl$attributes$CAPTURE_DATE)[5])

SAYAK11LIDs <- list(as.numeric(na.omit(SAYAK11LIDs)))
names(SAYAK11LIDs) <- "SAYAK11"

PoolCollections.GCL("SAYAK11",loci = loci96, IDs = SAYAK11LIDs, newname = "SAYAK11L")
SAYAK11L.gcl$n ## Should be 70
dput(x = SAYAK11L.gcl, file = "Raw genotypes/OriginalCollectionsNAremoved/SAYAK11L.txt")


# Next get early run fish from SAYAK12.gcl
table(SAYAK12.gcl$attributes$CAPTURE_DATE)
unique(SAYAK12.gcl$attributes$CAPTURE_DATE)[1:4]
SAYAK12EIDs <- AttributesToIDs.GCL(silly = "SAYAK12", attribute = "CAPTURE_DATE", matching = unique(SAYAK12.gcl$attributes$CAPTURE_DATE)[1:4])

SAYAK12EIDs <- list(as.numeric(na.omit(SAYAK12EIDs)))
names(SAYAK12EIDs) <- "SAYAK12"

PoolCollections.GCL("SAYAK12", loci = loci96, IDs = SAYAK12EIDs, newname = "SAYAK12E")
SAYAK12E.gcl$n ## Should be 200
dput(x = SAYAK12E.gcl, file = "Raw genotypes/OriginalCollectionsNAremoved/SAYAK12E.txt")


# Then get late run fish from SAYAK12.gcl
table(SAYAK12.gcl$attributes$CAPTURE_DATE)
unique(SAYAK12.gcl$attributes$CAPTURE_DATE)[10:14]
SAYAK12LIDs <- AttributesToIDs.GCL(silly = "SAYAK12", attribute = "CAPTURE_DATE", matching = unique(SAYAK12.gcl$attributes$CAPTURE_DATE)[10:14])

SAYAK12LIDs <- list(as.numeric(na.omit(SAYAK12LIDs)))
names(SAYAK12LIDs) <- "SAYAK12"

PoolCollections.GCL("SAYAK12", loci = loci96, IDs = SAYAK12LIDs, newname = "SAYAK12L")
SAYAK12L.gcl$n ## Should be 250
dput(x = SAYAK12L.gcl, file = "Raw genotypes/OriginalCollectionsNAremoved/SAYAK12L.txt")


## Auke Creek weir for parentage analysis, just get 200 representative picked by Serena
SAUKE13.gcl$n

SAUKE13baselineIDs <- dget(file = "Objects/SAUKE13baseline.txt")

PoolCollections.GCL("SAUKE13", loci = loci96, IDs = SAUKE13baselineIDs, newname = "SAUKE13baseline")
SAUKE13baseline.gcl$n ## 200
dput(x = SAUKE13baseline.gcl, file = "Raw genotypes/OriginalCollectionsNAremoved/SAUKE13baseline.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all = TRUE))
# This sources all of the new GCL functions to this workspace
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")
# This sources all of Kyle's non-GCL useful functions to this workspace (does not include "tempGCL" functions)
source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")

## Locus Control
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
loci96 <- dget(file = "Objects/loci96.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get collections
# KMA762Collections <- readClipboard()
# dput(x = KMA762Collections, file = "Objects/KMA762Collections.txt")

KMA762Collections <- dget(file = "Objects/KMA762Collections.txt")
KMA762Collections

require(beepr)
invisible(sapply(KMA762Collections, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/OriginalCollectionsNAremoved/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Make sure all collections are in pops ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
KMA762Collections
# KMA473Pops <- readClipboard()
# dput(x = KMA473Pops, file = "Objects/KMA473Pops.txt")
KMA473Pops <- dget(file = "Objects/KMA473Pops.txt")

table(KMA762Collections %in% unlist(strsplit(x = KMA473Pops, split = "\\.")))
table(unlist(strsplit(x = KMA473Pops, split = "\\.")) %in% KMA762Collections)
# All looks good, moving forward

KMA762CollectionsPopNumber <- sapply(KMA762Collections, function(collection) grep(pattern = collection, x = KMA473Pops))
writeClipboard(as.character(KMA762CollectionsPopNumber))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Comparing sample sizes between these collections and other existing baselines ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Questions for Tyler

# There are a few subtle (1-2 fish) differences in the duplicate check. WASSIP
# appears to be more conservative (has more duplicates). Is this because of the 
# "quantile" argument?

# Why include SHOOD00? This collection was specifically not included in WASSIP
# due to unreliable metadata/outlier status
# Checked Tyler's BB Baseline and it looks like he has tracked down the origin of these samples
# They are indeed from the Nelson River, but were sampled into SHOOD00 vials. Go ahead and use them.
# Also decided to include the collection SORZI12 (collected for SEDM to get Orzinski up to >400 fish
# so it could be its own RG). It pools to an existing pop SORZI00


## Questions for Andy
# Why include SCRES941, SCRES942, and SCRESL09? These three collections were specifically
# not included in wassip due to unreliable metadata/outlier status OR mixture status
# Andy has investigated, these collections and they appear to pool well with known collections
# and perform well in fisheries mixtures (i.e. when we have samples near Crescent, most of the fish allocated there)


## Questions for Serena

# Why include SMAIN91? This Main Bay hatchery collection was specifically not included in WASSIP.
# Serena says that this collection was added to the SEAK baseline due to non-convergence
# of some NSEAK mixtures, so it was later added to the baseline.

# What is up with the 2000+ samples for SAUKE13?
# These samples were taken for parentage analysis and having all of them in the baseline
# Led to HWE issues, so Serena used HWLR to pick 200 "representative" fish.

# Why are there sample size discrepancies with SCRES03, SSPEE03, and SSPEE13?
# Accounting error, this was a non-issue.


# Looking into large number of duplicates found for KANAL
# KMA745Collections_DuplicateCheckReportSummary$SKANA07
# KMA745Collections_DuplicateCheckReportSummary$SKANA10
# KMA745Collections_DuplicateCheckReportSummary$SKANAL13
# Kanalku is SUPER inbred and there are lots of "duplicates"
# Thus these three collections ought to be QCed without looking for duplicates

## Pre-QC Massage

# SAYAK 11 and SAYAK 12 have already been dealt with.
# Remember not to inlude "SKANA07"  "SKANA10"  "SKANAL13" in the duplicate check.
# SAUKE13baseline is already pooled, don't use SAUKE13.
# SCRES03, SSPEE03, and SPEE13 check out. There must have been an accounting error.
# Check with Andy about inclusion of weird Crescent Lake collections ("SCRES941", "SCREE942", and "SCRESL09"). They are good.
# Check with Tyler about SHOOD00, include it, the metadata is reliable.
# Add SORZI12
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Data QC/Massage ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx); require(beepr)

KMA762Collections

KMA762Collections_SampleSizes <- matrix(data = NA, nrow = length(KMA762Collections), ncol = 5, 
                                         dimnames = list(KMA762Collections, c("Genotyped", "Alternate", "Missing", "Duplicate", "Final")))

#### Check loci
## Get sample size by locus
Original_KMA762Collections_SampleSizebyLocus <- SampSizeByLocus.GCL(KMA762Collections, loci96)
min(Original_KMA762Collections_SampleSizebyLocus)  # 2
sort(apply(Original_KMA762Collections_SampleSizebyLocus,1,min)/apply(Original_KMA762Collections_SampleSizebyLocus,1,max))  # Several under 0.8
table(apply(Original_KMA762Collections_SampleSizebyLocus,1,min)/apply(Original_KMA762Collections_SampleSizebyLocus,1,max) < 0.8)  # 28 SILLY's with at least one locus fail
str(Original_KMA762Collections_SampleSizebyLocus)

## Percent by locus
Original_KMA762Collections_PercentbyLocus <- apply(Original_KMA762Collections_SampleSizebyLocus, 1, function(row) {row / max(row)})
which(apply(Original_KMA762Collections_PercentbyLocus, 2, min) < 0.8)
# writeClipboard(as.character(apply(Original_KMA762Collections_PercentbyLocus, 2, min) < 0.8))

require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(Original_KMA762Collections_PercentbyLocus), col.regions = new.colors, xlab = "SILLY", ylab = "Locus", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares
## This looks pretty good, no holes are evident after tweaking ReadLOKIv2.GCL to ReadLOKIv3.GCL with "onlymarkersuitefish" switch turned to TRUE


#### Check individuals
### Initial
## Get number of individuals per silly before removing missing loci individuals
Original_KMA762Collections_ColSize <- sapply(paste(KMA762Collections, ".gcl", sep = ''), function(x) get(x)$n)
KMA762Collections_SampleSizes[, "Genotyped"] <- Original_KMA762Collections_ColSize


### Alternate
## Indentify alternate species individuals
ptm <- proc.time()
KMA762Collections_Alternate <- FindAlternateSpecies.GCL(sillyvec = KMA762Collections, species = "sockeye"); beep(8)
proc.time() - ptm

## Remove Alternate species individuals
RemoveAlternateSpecies.GCL(AlternateSpeciesReport = KMA762Collections_Alternate, AlternateCutOff = 0.5, FailedCutOff = 0.5); beep(2)

## Get number of individuals per silly after removing alternate species individuals
ColSize_KMA762Collections_PostAlternate <- sapply(paste(KMA762Collections, ".gcl", sep = ''), function(x) get(x)$n)
KMA762Collections_SampleSizes[, "Alternate"] <- Original_KMA762Collections_ColSize-ColSize_KMA762Collections_PostAlternate


### Missing
## Remove individuals with >20% missing data
KMA762Collections_MissLoci <- RemoveIndMissLoci.GCL(sillyvec = KMA762Collections, loci = loci96, proportion = 0.8); beep(8)

## Get number of individuals per silly after removing missing loci individuals
ColSize_KMA762Collections_PostMissLoci <- sapply(paste(KMA762Collections, ".gcl", sep = ''), function(x) get(x)$n)
KMA762Collections_SampleSizes[, "Missing"] <- ColSize_KMA762Collections_PostAlternate-ColSize_KMA762Collections_PostMissLoci


### Duplicate
## Check within collections for duplicate individuals.
KMA762Collections_DuplicateCheck95MinProportion <- CheckDupWithinSilly.GCL(sillyvec = KMA762Collections, loci = loci96, quantile = NULL, minproportion = 0.95); beep(8)
KMA762Collections_DuplicateCheckReportSummary <- sapply(KMA762Collections, function(x) KMA762Collections_DuplicateCheck95MinProportion[[x]]$report)

## Remove duplicate individuals
# Don't include "SKANA07"  "SKANA10"  "SKANAL13" (inbred/bottlenecked)
KMA762Collections_RemovedDups <- RemoveDups.GCL(KMA762Collections_DuplicateCheck95MinProportion[-grep(pattern = "SKANA", x = KMA762Collections)])

## Get number of individuals per silly after removing duplicate individuals
ColSize_KMA762Collections_PostDuplicate <- sapply(paste(KMA762Collections, ".gcl", sep = ''), function(x) get(x)$n)
KMA762Collections_SampleSizes[, "Duplicate"] <- ColSize_KMA762Collections_PostMissLoci-ColSize_KMA762Collections_PostDuplicate


### Final
KMA762Collections_SampleSizes[, "Final"] <- ColSize_KMA762Collections_PostDuplicate
KMA762Collections_SampleSizes

write.xlsx(KMA762Collections_SampleSizes, file = "Output/KMA762Collections_SampleSizes.xlsx")

## Save post-QC .gcl's as back-up:
dir.create(path = "Raw genotypes/PostQCCollections")
invisible(sapply(KMA762Collections, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostQCCollections/" , silly, ".txt", sep = ''))} )); beep(8)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Checking HWE in Collections ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Check HWE in collections (only use diploid loci)
mito.loci <- dget(file = "Objects/mito.loci.txt")

gcl2Genepop.GCL(sillyvec = KMA762Collections, loci = loci96[-mito.loci], path = "Genepop/KMA762Collections_93nuclearloci.txt", VialNums = TRUE); beep(2)

## Check HWE in Genepop
# Originally done with Exact test and subsequently re-done with MCMC
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")
HWE <- ReadGenepopHWE.GCL(file = "Genepop/KMA762Collections_93nuclearloci_MCMC.txt.P")
str(HWE)

ncollections <- length(KMA762Collections)
nloci <- length(loci96[-mito.loci])

# Summary p-values
HWE$SummaryPValues

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Overall loci
t(t(HWE$SummaryPValues[1:nloci, "Overall Pops"]))
sort(HWE$SummaryPValues[1:nloci, "Overall Pops"])  # Issues with One_c3-98 (known culprit in SEAK, cluster of red homozygotes often no-called)

# Markers that are out of HWE with p < 0.01 over all Pops (Fisher's method)
table(HWE$SummaryPValues[1:nloci, "Overall Pops"] < 0.01)  # Number of markers that are p < 0.01 over all pops (Fisher's method)
HWE$SummaryPValues[1:nloci, "Overall Pops"][HWE$SummaryPValues[1:nloci, "Overall Pops"] < 0.01]  # One_c3-98 & One_ACBP-79

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WASSIP approach was to toss out any locus that is p < 0.01 over all Populations (Fisher's method)
KMA762Collections_HWE_locusfail <- loci96[-mito.loci][HWE$SummaryPValues[1:nloci, "Overall Pops"] < 0.01]
dput(x = KMA762Collections_HWE_locusfail, file = "Objects/KMA762Collections_HWE_locusfail.txt")
# KMA762Collections_HWE_locusfail <- loci96[-mito.loci][HWE[[2]][seq(loci96[-mito.loci]), length(colnames(HWE[[2]]))] < 0.01] # 11 markers if done by Exact test
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# How many collections is a marker p < 0.05
sort(apply(HWE$SummaryPValues[1:nloci, 1:ncollections], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} ))  # One_c3-98 out of HWE in 93 / 762 Pops

# Histogram of # pops out of HWP p < 0.05 per loci
par(mar = c(5.1, 5.1, 4.1, 2.1))
hist(apply(HWE$SummaryPValues[1:nloci, 1:ncollections], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} ), col = 8, xlab = "# Collections p < 0.05", ylab = "# Loci", main = "HWE for 762 Collections", breaks = seq(0, 110, 1), cex.main = 2, cex.lab = 1.5)  # How many populations < 0.05 per marker
abline(v = 0.05 * length(KMA762Collections), col = "red", lwd = 5)  # Expectation that it should be p < 0.05 if 5% of Collections
text(x = 0.05 * length(KMA762Collections), y = 8.2, labels = "5% of Collections", col = "red", pos = 4, cex = 1.5)
text(x = unique(sort(apply(HWE$SummaryPValues[1:nloci, 1:ncollections], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} )))[29:35] - 0.5,
     y = table(sort(apply(HWE$SummaryPValues[1:nloci, 1:ncollections], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} )))[29:35],
     labels = sapply(unique(sort(apply(HWE$SummaryPValues[1:nloci, 1:ncollections], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} )))[29:35], function(ncollectionsout) {paste(names(which(apply(HWE$SummaryPValues[1:nloci, 1:ncollections], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} ) == ncollectionsout)), collapse = ", ")} ),
     srt = 38, pos = 4, cex = 1.2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Overall collections
t(t(HWE$SummaryPValues["Overall Loci", 1:ncollections]))
sort(HWE$SummaryPValues["Overall Loci", 1:ncollections])

# Collections that are out of HWE with p < 0.01 over all Pops (Fisher's method)
table(HWE$SummaryPValues["Overall Loci", 1:ncollections] < 0.01)  # 7 collections that are p < 0.01 over all pops (Fisher's method)
HWE$SummaryPValues["Overall Loci", 1:ncollections][HWE$SummaryPValues["Overall Loci", 1:ncollections] < 0.01]  # SLYNXLK09_95    SCHAU01_96   SBATL04B_96    SSKIL92_96     SEEK04_32    SBAKE96_97   SCEDAR94_95


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WASSIP approach was to toss out any population that is p < 0.05 for > 10% loci
sort(apply(HWE$SummaryPValues[1:nloci, 1:ncollections], 2, function(pop) {sum(pop < 0.05, na.rm = TRUE)} ))
KMA762Collections_HWE_popsfail <- KMA762Collections[apply(HWE$SummaryPValues[1:nloci, 1:ncollections], 2, function(pop) {sum(pop < 0.05, na.rm = TRUE)} ) > (0.1 * nloci)]
dput(x = KMA762Collections_HWE_popsfail, file = "Objects/KMA762Collections_HWE_popsfail.txt")
# NOTE: WE ARE NOT GOING TO TOSS OUT ANY POPULATIONS HERE, DESPITE "RULE BREAKING"
# The decision was made to keep ALL pops from regional baselines, presuming that
# project leaders have already gone through this process and there is no need to
# re-invent the wheel.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


write.xlsx(x = HWE$SummaryPValues, file = "Output/KMA762Collections_HWEmatrix_93nuclearloci_MCMC_Genepop.xlsx")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## View patterns of HWE
# NOTE: white is NA
require(lattice)
# KMA762CollectionsGroupVec15 <- as.numeric(readClipboard())
# dput(x = KMA762CollectionsGroupVec15, file = "Objects/KMA762CollectionsGroupVec15.txt")
KMA762CollectionsGroupVec15 <- dget(file = "Objects/KMA762CollectionsGroupVec15.txt")
Colors15 <- dget(file = "Objects/Colors15.txt")

# 762 Collections 93 nuclear loci + Overall Collections/Loci
levelplot(t(HWE$SummaryPValues), 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA762CollectionsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.75, x1 = cumsum(table(KMA762CollectionsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.75, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA762CollectionsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 94.25, x1 = cumsum(table(KMA762CollectionsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 94.25, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          col.regions = c("red", "darkorange", "yellow", "green3", "grey90"), xlab = "SILLY", ylab = "Locus", 
          at = c(0, 0.001, 0.01, 0.05, 1), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares

# 762 Collections 2 HWE FAIL markers (for pops, because only One_c3-98 failed for collections)
levelplot(t(HWE$SummaryPValues)[, c(loci96[-mito.loci] %in% c(KMA762Collections_HWE_locusfail, "One_ACBP-79"), FALSE)], 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA762CollectionsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.45, x1 = cumsum(table(KMA762CollectionsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.45, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA762CollectionsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 2.55, x1 = cumsum(table(KMA762CollectionsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 2.55, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          col.regions = c("red", "darkorange", "yellow", "green3", "grey90"), xlab = "SILLY", ylab = "Locus", 
          at = c(0, 0.001, 0.01, 0.05, 1), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares



## View distribution of HWE p-values
# By Locus
maxylim <- max(sapply(loci96[-mito.loci], function(loci) {sum(HWE$SummaryPValues[loci, 1:ncollections] > 0.95, na.rm = TRUE)}))  # Determine ylim
pdf(file = "HWE/HWEPvalueHistograms_byLocus_KMA762Collections_93nuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)   
invisible(sapply(loci96[-mito.loci], function(loci) {hist(HWE$SummaryPValues[loci, 1:ncollections], col = 8, breaks = seq(from = 0, to = 1, by = 0.05), main = loci, xlab = "P-value", ylim = c(0, maxylim))
                                                     abline(h = ncollections / 20, col = "red", lwd = 3)} ))
dev.off()

# Only 2 problem markers
maxylim <- max(sapply(c(KMA762Collections_HWE_locusfail, "One_ACBP-79"), function(loci) {sum(HWE[[2]][loci, seq(dim(HWE[[2]])[2] - 1)] > 0.95, na.rm = TRUE)}))  # Determine ylim
pdf(file = "HWE/HWEPvalueHistograms_byLocus_KMA762Collections_2failnuclearloci_MCMC_Genepop.pdf",width=11, height=8.5, family="Times",pointsize=20)   
invisible(sapply(c(KMA762Collections_HWE_locusfail, "One_ACBP-79"), function(loci) {hist(HWE$SummaryPValues[loci, 1:ncollections], col = 8, breaks = seq(from = 0, to = 1, by = 0.05), main = loci, xlab = "P-value", ylim = c(0, maxylim))
                                                                                    abline(h = ncollections / 20, col = "red", lwd = 3)} ))
dev.off()


# By Collections
maxylim <- max(sapply(seq(ncollections), function(silly) {sum(HWE$SummaryPValues[1:nloci, silly] > 0.95, na.rm = TRUE)}))  # Determine ylim
pdf(file = "HWE/HWEPvalueHistograms_byPop_KMA762Collections_93nuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
invisible(sapply(seq(ncollections), function(silly) {hist(HWE$SummaryPValues[1:nloci, silly], col = 8, breaks = seq(from = 0, to = 1, by = 0.05), main = KMA762Collections[silly], xlab = "P-value", ylim = c(0, maxylim))
                                                     abline(h = nloci / 20, col = "red", lwd = 3)} ))
dev.off()

# Only 11 failed Collections
KMA762Collections_HWE_popsfail_number <- sapply(KMA762Collections_HWE_popsfail, function(pop) {which(KMA762Collections == pop)} )

maxylim <- max(sapply(KMA762Collections_HWE_popsfail_number, function(silly) {sum(HWE$SummaryPValues[1:nloci, silly] > 0.95, na.rm = TRUE)}))  # Determine ylim
pdf(file = "HWE/HWEPvalueHistograms_byPop_KMA11CollectionssFail_93nuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
invisible(sapply(KMA762Collections_HWE_popsfail_number, function(silly) {hist(HWE$SummaryPValues[1:nloci, silly], col = 8, breaks = seq(from = 0, to = 1, by = 0.05), main = KMA762Collections[silly], xlab = "P-value", ylim = c(0, maxylim))
                                                                         abline(h = nloci / 20, col = "red", lwd = 3)} ))
dev.off()
# There are a few populations that look suspicious, but I will lean on previous baseline work and keep them in
# No need to re-invent the wheel

# Where are these 11 Collections from, just out of curiousity
Groups15 <- dget(file = "Objects/Groups15.txt")
setNames(object = Groups15[KMA762CollectionsGroupVec15[KMA762Collections_HWE_popsfail_number]], nm = KMA762Collections_HWE_popsfail)



## View Fis vs. Pvalue for all markers
# All 93 nuclear markers
Groups15Short <- dget(file = "Objects/Groups15Short.txt")

pdf(file = "HWE/HWEPvalueVsFis_byPop_KMA762Collections_93nuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
invisible(sapply(loci96[-mito.loci], function(locus) {plot(x = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"], y = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"], pch = 16, xlim = c(-1, 1), ylim = c(0, 1), xlab = "WC Fis", ylab = "P-value", main = locus, col = Colors15[KMA762CollectionsGroupVec15], cex = 0.7)
                                                      legend("topright", legend = Groups15Short, fill = Colors15, xpd = TRUE, bty = "n", cex = 0.7)
                                                      abline(h = 0.05, col = "red", lwd = 3)
                                                      text(x = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"][which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05)], y = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"][which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05)],
                                                           labels = which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05), cex = 0.4)} ))
dev.off()


# Only 2 problem markers
pdf(file = "HWE/HWEPvalueVsFis_byPop_KMA762Collections_2failnuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
invisible(sapply(c(KMA762Collections_HWE_locusfail, "One_ACBP-79"), function(locus) {plot(x = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"], y = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"], pch = 16, xlim = c(-1, 1), ylim = c(0, 1), xlab = "WC Fis", ylab = "P-value", main = locus, col = Colors15[KMA762CollectionsGroupVec15], cex = 0.7)
                                                      legend("topright", legend = Groups15Short, fill = Colors15, xpd = TRUE, bty = "n", cex = 0.7)
                                                      abline(h = 0.05, col = "red", lwd = 3)
                                                      text(x = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"][which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05)], y = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"][which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05)],
                                                           labels = which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05), cex = 0.4)} ))
dev.off()
# There are a few collections that look suspicious, but I will lean on previous baseline work and keep them in

# Histograms of Fis by marker for 762 Collections, out of curiousity
pdf(file = "HWE/HistogramofFis_byLocus_KMA762Collections_93nuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
invisible(sapply(loci96[-mito.loci], function(locus) {hist(x = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"], breaks = seq(-1, 1, 0.05), col = 8, xlab = "WC Fis", main = locus)
                                                       abline(v = median(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"], na.rm = TRUE), lwd = 4, col = "red")
                                                       abline(v = mean(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"], na.rm = TRUE), lwd = 4, col = "blue")
                                                       legend("topright", legend = c("Median", "Mean"), col = c("red", "blue"), lwd = 5, cex = 1.5, bty = "n")} ))
dev.off()

# Mean Fis for each marker over 762 Collections
hist(sapply(loci96[-mito.loci], function(locus) {mean(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"], na.rm = TRUE)} ), col = 8, breaks = seq(-0.13, 0.13, 0.005), xlab = "WC Fis", main = "Histogram of Mean Fis over 762 Sockeye Collections\nfor 93 nuclear SNPs")
abline(v = 0, col = "black", lwd = 4)
text(x = min(sapply(loci96[-mito.loci], function(locus) {mean(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"], na.rm = TRUE)} )), y = 1, labels = "One_c3-98", cex = 1.5, pos = 4, srt = 15)


# Digging into "One_ACBP-79"
HWE$DataByPop$Silly <- rep(x = KMA762Collections, each = length(loci96[-mito.loci]))  # Add column for silly, as "Pop" is the genepop bastardization
HWEfail_ACBP_sillys <- sapply(HWE$DataByPop[which(HWE$DataByPop$Locus == "One_ACBP-79" & HWE$DataByPop$PValue < 0.05), "Silly"], function(silly) {unlist(strsplit(x = silly, split = "_"))[1]} )  # all sillys with p < 0.05 for ACBP
HWEfail_ACBP_sillys_FreqFis <- FreqFisPlot4SNPs.GCL(sillyvec = HWEfail_ACBP_sillys, loci = loci96[-mito.loci], groupvec = KMA762CollectionsGroupVec15[KMA762Collections %in% HWEfail_ACBP_sillys], groupcol = Colors15, file = "HWE/FreqFisPlot_HWE_collectionsfail_ACBP.pdf")
cbind(HWEfail_ACBP_sillys,Groups15[KMA762CollectionsGroupVec15[KMA762Collections %in% HWEfail_ACBP_sillys]])
# It looks like most of the time when One_ACBP-79 is out of HWE, Fis is + (missing hets), especially in SW Kodiak
# I looked at some of the Fluidigm plots for S041, where these silly's were run and there does not appear to be any systematic no-calling of hets
# This marker is legitmately not conforming to HWE, and should be dropped

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Allele Frequency Plots Collections ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### All 762 Collections
# Counts
KMA762CollectionsPostRemovalsAlleleCounts <- FreqPop.GCL(sillyvec = KMA762Collections, loci = loci96)
str(KMA762CollectionsPostRemovalsAlleleCounts)

# Frequencies
KMA762CollectionsPostRemovalsFreqs <- KMA762CollectionsPostRemovalsAlleleCounts[,,"Allele 1"] / (KMA762CollectionsPostRemovalsAlleleCounts[,,"Allele 2"] + KMA762CollectionsPostRemovalsAlleleCounts[,,"Allele 1"])
str(KMA762CollectionsPostRemovalsFreqs)
write.xlsx(KMA762CollectionsPostRemovalsFreqs, file = "Output/KMA762CollectionsPostRemovalsFreqs.xlsx")

#### Create plots
## Create groupvec
# KMA762CollectionsGroups15 <- readClipboard()
# 
# Groups15 <- unique(KMA762CollectionsGroups15)
# dput(x = Groups15, file = "Objects/Groups15.txt")
# 
# Groups15Short <- c("W of Chig", "Black L", "Chig L", "UStat/Ak", "Frazer", "Ayak", "Karluk", "Uganik", "NW Minor", "Afognak", "E Minor", "Saltery", "CI", "PWS/Cop", "SEAK")
# dput(x = Groups15Short, file = "Objects/Groups15Short.txt")
# 
# KMA762CollectionsGroupVec15 <- rep(x = seq(Groups15), times = table(KMA762CollectionsGroups15)[Groups15])
# dput(x = KMA762CollectionsGroupVec15, file = "Objects/KMA762CollectionsGroupVec15.txt")
# writeClipboard(as.character(KMA762CollectionsGroupVec15))
# 
# Colors15 <- c("grey60", "black", "brown", "red", "cyan", "blue", "purple", "darkorange2", "magenta", "yellow3", "green", "darkgreen", "pink2", "darkgoldenrod4", "turquoise4")
# plot(seq(Colors15), col = Colors15, pch = 16, cex = 5)
# dput(x = Colors15, file = "Objects/Colors15.txt")

Groups15 <- dget(file = "Objects/Groups15.txt")
Groups15Short <- dget(file = "Objects/Groups15Short.txt")
KMA762CollectionsGroupVec15 <- dget(file = "Objects/KMA762CollectionsGroupVec15.txt")
Colors15 <- dget(file = "Objects/Colors15.txt")

## Plot points and smoothing line
getwd()
pdf(file = "FreqPlots/KMA762Collections_96SNPs_GroupVec15_FreqPlots.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
par(mar = c(2.1, 4.1, 4.1, 4.6))
for(locus in loci96){
  plot(KMA762CollectionsPostRemovalsFreqs[, locus], main = locus, col = Colors15[KMA762CollectionsGroupVec15], pch = 19, ylim = c(0, 1), ylab = "Frequency", xlab = "Collections", xaxt = 'n', cex = 0.7)
  lines(supsmu(seq(length(KMA762Collections)), KMA762CollectionsPostRemovalsFreqs[, locus]), lwd = 2)
  mtext(text = "Collections", side = 1, line = 0.7)
  text(x = seq(KMA762Collections), y = KMA762CollectionsPostRemovalsFreqs[, locus], labels = seq(KMA762Collections), cex = 0.2)
  legend(x = 780, y = 0.9, legend = Groups15Short, fill = Colors15, xpd = TRUE, bty = "n", cex = 0.7)
}; rm(locus)
dev.off()



## View levelplot of frequencies
require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(KMA762CollectionsPostRemovalsFreqs), col.regions = new.colors, xlab = "Locus", ylab = "SILLY", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares
levelplot(KMA762CollectionsPostRemovalsFreqs, col.regions = new.colors, xlab = "SILLY", ylab = "Locus", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares

#~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!!
#### Allele 1/2 are different between LOKI2R.GCL and ReadLOKI.GCL!!!! ####
#~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!~~~!!!!
# 10/29/215, this was fixed
LocusControl$alleles$One_U1101
LocusControl$alleles


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pooling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save first
## save.image("V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/KMASockeyeBaseline.RData")

## All pooling is from previous baselines, there is NO new pooling
KMA473pops <- dget(file = "Objects/KMA473Pops.txt")

## Confirm that all Collections are present in the Populations
table(KMA762Collections %in% unlist(strsplit(x = paste(KMA473Pops, collapse = "."), split = "\\.")))

## List of Pops to pool
KMA473PopstoCollections <- strsplit(x = KMA473Pops, split = "\\.")
KMAPool <- KMA473PopstoCollections[sapply(KMA473PopstoCollections, function(pop) {length(pop) > 1} )]

## Pooling
invisible(lapply(KMAPool, function(pop) {PoolCollections.GCL(collections = pop, loci = loci96, IDs = NULL, newname = paste(pop, collapse = "."))} )); beep(8)

## Confirm that all Populations exist as .gcl objects
table(KMA473Pops %in% unlist(strsplit(x = objects(pattern = "\\.gcl"), split = ".gcl")))

## Save post-QC Pop .gcl's as back-up:
dir.create(path = "Raw genotypes/PostQCPops")
invisible(sapply(KMA473Pops, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostQCPops/" , silly, ".txt", sep = ''))} )); beep(8)

#~~~~~~~~~~~~~~~~~~
## Investigation of temporal English Bay 1992 collections...
EnglishBay.Collections <- c("SENG92E", "SENG92L")
EnglishBay96SNPFreq <- FreqPop.GCL(sillyvec = EnglishBay.Collections, loci = loci96)  # Quick look

TemporalEnglishBayTest <- list(EnglishBay.Collections)
TemporalEnglishBayTestResult <- FishersTest.GCL(freq = EnglishBay96SNPFreq, loci = loci96, tests = TemporalEnglishBayTest)
TemporalEnglishBayTestResult$OverallResults
hist(TemporalEnglishBayTestResult$ResultsByLocus$SENG92E.SENG92L$pval, breaks = seq(0, 1, 0.05), col = 8, main = "Histogram of Fisher P-values for 96 SNPs")  # histogram of p-values, more than expected
hist(sort(apply(sapply(unlist(TemporalEnglishBayTest), function(silly) {apply(EnglishBay96SNPFreq[silly, , ], 1, function(locus) {locus["Allele 1"] / sum(locus)} )} ), 1, function(locus) {locus[1] - locus[2]} )), col = 8, main = "Histogram of allele frequency differences 96 SNPs", xlab = "Allele Freq Differences")
sapply(EnglishBay.Collections, function(silly) {get(paste(silly, ".gcl", sep = ''))$n} )

# Also ran fullsniplings since these are hatchery smolt collections!
# SENG92E table was 1 - 78; 2 - 4 (i.e. most individuals are completely unrelated)
# SENG92L table was 1 - 40; 2 - 14; 3 - 6; 4 - 2 (i.e. more relationships, but nothing too crazy, however, the bigger fullsib groupsing had the most support in MCMC reps)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pooling for new colletions ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Need to figure out where "new" Kodiak and SEAK collections pool
## To get a final, ordered KMA762CollectionsOrdered

KMA762Collections96SNPsFreq <- FreqPop.GCL(sillyvec = KMA762Collections, loci = loci96)  # Quick look
mito.loci <- dget(file = "Objects/mito.loci.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Kodiak Baseline Collections

#~~~~~~~~~~~~~~~~~~
## SCAIDA14 Regional
Frazer.Collections <- c("SPINNM08", "SSTUM08", "SCOUR08", "SMIDWS08", "SHOLFS08", "SMIDWM08", "SLINDM08", "SVALA08", "SOUTS08", "SDOGSC08", "SCAIDA14")
FrazerTest <- combn(x = Frazer.Collections, m = 2, simplify = FALSE)
FrazerTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = FrazerTest)

sapply(Frazer.Collecionts, function(silly) {FrazerTestResults$OverallResults[grep(pattern = silly, x = rownames(FrazerTestResults$OverallResults)), ]}, USE.NAMES = TRUE, simplify = FALSE)

# New pooling DO NOT USE
FrazerPoolTest <- list(c("SPINNM08", "SSTUM08"), c("SMIDWM08", "SLINDM08", "SHOLFS08", "SDOGSC08", "SVALA08"), c("SCOUR08", "SMIDWS08"))
FrazerPoolTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = FrazerPoolTest)
FrazerPoolTestResults$OverallResults  # SOUTS08 & SCAIDA14 will remain alone in pops by themselves.

# WASSIP pooling, done by ecotype USE THIS
FrazerWASSIPPoolTest <- list(c("SPINNM08", "SSTUM08"), c("SCOUR08", "SMIDWS08", "SHOLFS08"), c("SMIDWM08", "SLINDM08"))
FrazerWASSIPPoolTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = FrazerWASSIPPoolTest)
FrazerWASSIPPoolTestResults$OverallResults  # 0.9819; 0.0462; 0.9995
# Leave as per WASSIP and add SCAIDA14 as its own pop.

source("V:/DATA/R_GEN/JJs GCL/FreqFisPlot4SNPs.GCL.R")
FrazerFreqFis <- FreqFisPlot4SNPs.GCL(sillyvec = Frazer.Collections, loci = loci96[-mito.loci], groupvec = seq(Frazer.Collections), file = "FreqFisPlots for Pooling/FrazerFreqFisPlot.pdf")
str(FrazerFreqFis)
sapply(rownames(FrazerFreqFis$HWE.pval), function(silly) {hist(FrazerFreqFis$HWE.pval[silly, ], breaks = seq(0, 1, 0.05), col = 8, ylim = c(0, 20), main = silly)
                                                          abline(h = 5, col = 2, lwd = 2)} )

#~~~~~~~~~~~~~~~~~~
## SUGAN15 Temporal
TemporalUganikTest <- list(c("SUGAN97", "SUGAN15"))
TemporalUganikTestReulst <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalUganikTest)
TemporalUganikTestReulst$OverallResults  # 0.0102
hist(TemporalUganikTestReulst$ResultsByLocus$SUGAN15.SUGAN97$pval, breaks = seq(0, 1, 0.05), col = 8)  # histogram of p-values, more than expected
# Sorted differences in allele frequency
sort(apply(sapply(unlist(TemporalUganikTest), function(silly) {apply(KMA762Collections96SNPsFreq[silly, , ], 1, function(locus) {locus["Allele 1"] / sum(locus)} )} ), 1, function(locus) {locus[1] - locus[2]} ))
# One_U1004-183 is pretty different (0.2)
gcl2Genepop.GCL(sillyvec = unlist(TemporalUganikTest), loci = loci96[-mito.loci], path = "Genepop/UganikCollections_93nuclearloci.txt")
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")
HWE <- ReadGenepopHWEKS.GCL(file = "Genepop/UganikCollections_93nuclearloci.txt.P")
HWE$SummaryPValues["Overall Loci", ]
sapply(colnames(HWE$SummaryPValues), function(silly) {hist(HWE$SummaryPValues[, silly], breaks = seq(0, 1, by = 0.05), col = 8, main = silly, ylim = c(0,40))
                                                      abline(h = 5, col = 2, lwd = 2)} )
sort(HWE$SummaryPValues[, "SUGAN15_190 "])
# not out of HWE, just separate populations?, do not pool.
# Spoke with Birch Foster on 10/21/15 and he indicated that multiple spawning pops might be possible

# Pool collections anyways to check and see how they look pooled
PoolCollections.GCL(collections = unlist(TemporalUganikTest), loci = loci96, IDs = NULL, newname = paste(unlist(TemporalUganikTest), collapse = "."))
source("V:/DATA/R_GEN/JJs GCL/FreqFisPlot4SNPs.GCL.R")
UganikFreqFis <- FreqFisPlot4SNPs.GCL(sillyvec = c("SUGAN97", "SUGAN15", "SUGAN97.SUGAN15"), loci = loci96[-mito.loci], groupvec = 1:3, file = "FreqFisPlots for Pooling/UganikFreqFisPlot.pdf")
str(UganikFreqFis)
sapply(rownames(UganikFreqFis$HWE.pval), function(silly) {hist(UganikFreqFis$HWE.pval[silly, ], breaks = seq(0, 1, 0.05), col = 8, ylim = c(0, 20), main = silly)
                                                          abline(h = 5, col = 2, lwd = 2)} )

#~~~~~~~~~~~~~~~~~~
## SBARAB12 & 15 Temporal
TemporalBarabraTest <- list(c("SBARAB12", "SBARAB15"))
TemporalBarabraTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalBarabraTest)
TemporalBarabraTestResults$OverallResults  # 0.9961
hist(sort(TemporalBarabraTestResults$ResultsByLocus$SBARAB12.SBARAB15$pval), breaks = seq(0, 1, 0.05), col = 8)
# Pool all temporal collections


#~~~~~~~~~~~~~~~~~~
## SBUSKL15 Temporal
TemporalBuskinTest <- list(c("SBUSK05", "SBUSKL10", "SBUSKL15"))
TemporalBuskinTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalBuskinTest)
TemporalBuskinTestResults$OverallResults  # 0.6156
hist(sort(TemporalBuskinTestResults$ResultsByLocus$SBUSK05.SBUSKL10.SBUSKL15$pval), breaks = seq(0, 1, 0.05), col = 8)

PoolCollections.GCL(collections = unlist(TemporalBuskinTest), loci = loci96, IDs = NULL, newname = paste(unlist(TemporalBuskinTest), collapse = "."))
source("V:/DATA/R_GEN/JJs GCL/FreqFisPlot4SNPs.GCL.R")
BuskinFreqFis <- FreqFisPlot4SNPs.GCL(sillyvec = c("SBUSK05", "SBUSKL10", "SBUSKL15", "SBUSK05.SBUSKL10.SBUSKL15"), loci = loci96[-mito.loci], groupvec = 1:4, file = "FreqFisPlots for Pooling/BuskinFreqFisPlot.pdf")
str(BuskinFreqFis)
sapply(rownames(BuskinFreqFis$HWE.pval), function(silly) {hist(BuskinFreqFis$HWE.pval[silly, ], breaks = seq(0, 1, 0.05), col = 8, ylim = c(0, 20), main = silly)
                                                          abline(h = 5, col = 2, lwd = 2)} )
# Pool all temporal collections


#~~~~~~~~~~~~~~~~~~
## SLKLOU14 Temporal
TemporalLakeLouiseTest <- list(c("SLKLOU05", "SLKLOU10", "SLKLOU14"))
TemporalLakeLouiseTestReults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalLakeLouiseTest)
TemporalLakeLouiseTestReults$OverallResults  # 0.0351; 10 and 14 are 0.7311; 05 and 10 are 0.9353; but 05 and 14 are 0.0372
TemporalLakeLouiseTestReults$ResultsByLocus$SLKLOU05.SLKLOU10.SLKLOU14[order(TemporalLakeLouiseTestReults$ResultsByLocus$SLKLOU05.SLKLOU10.SLKLOU14$pval), ]
hist(sort(TemporalLakeLouiseTestReults$ResultsByLocus$SLKLOU05.SLKLOU10.SLKLOU14$pval), breaks = seq(0, 1, 0.05), col = 8)
# Should we pool? Stick with 05/10 pooled as per WASSIP, but leave 14 separate.

# Check HWE
gcl2Genepop.GCL(sillyvec = unlist(TemporalLakeLouiseTest), loci = loci96[-mito.loci], path = "Genepop/LakeLouiseCollections_93nuclearloci.txt")
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")
HWE <- ReadGenepopHWEKS.GCL(file = "Genepop/LakeLouiseCollections_93nuclearloci.txt.P")
HWE$SummaryPValues["Overall Loci", ] # all pops appear to be conform to HWE
sapply(colnames(HWE$SummaryPValues), function(silly) {hist(HWE$SummaryPValues[, silly], breaks = seq(0, 1, by = 0.05), col = 8, main = silly, ylim = c(0,40))
                                                      abline(h = 5, col = 2, lwd = 2)} )  # looks good, no issues

PoolCollections.GCL(collections = unlist(TemporalLakeLouiseTest), loci = loci96, IDs = NULL, newname = paste(unlist(TemporalLakeLouiseTest), collapse = "."))
source("V:/DATA/R_GEN/JJs GCL/FreqFisPlot4SNPs.GCL.R")
LakeLouiseFreqFis <- FreqFisPlot4SNPs.GCL(sillyvec = c("SLKLOU05", "SLKLOU10", "SLKLOU14", "SLKLOU05.SLKLOU10.SLKLOU14"), loci = loci96[-mito.loci], groupvec = 1:4, file = "FreqFisPlots for Pooling/LakeLouiseFreqFisPlot.pdf")
str(LakeLouiseFreqFis)
sapply(rownames(LakeLouiseFreqFis$HWE.pval), function(silly) {hist(LakeLouiseFreqFis$HWE.pval[silly, ], breaks = seq(0, 1, 0.05), col = 8, ylim = c(0, 20), main = silly)
                                                              abline(h = 5, col = 2, lwd = 2)} )
sapply(rownames(LakeLouiseFreqFis$HWE.pval), function(silly) {pchisq(q = -2 * sum(log(LakeLouiseFreqFis$HWE.pval[silly, ]), na.rm = TRUE), df = 2 * sum(!is.na(LakeLouiseFreqFis$HWE.pval[silly, ])), lower.tail = FALSE)} )
HWE$SummaryPValues["Overall Loci", ]

cbind(t(round(LakeLouiseFreqFis$HWE.pval, 4))[, 1:3], HWE$SummaryPValues[1:93, 1:3])  # Not exactly the same because Genepop assigns NA and this program assigns 1.000
t(LakeLouiseFreqFis$allele.freqs)


#~~~~~~~~~~~~~~~~~~
## SSALT14 Temporal
TemporalSalteryTest <- list(c("SSALT94", "SSALT99", "SSALT14"))
TemporalSalteryTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalSalteryTest)
TemporalSalteryTestResults$OverallResults  # 0.0058; 94 and 99 0.6766; 99 and 14 0.0256; 94 and 14 0.5511

TemporalSalteryTestResults$ResultsByLocus$SSALT94.SSALT99.SSALT14[order(TemporalSalteryTestResults$ResultsByLocus$SSALT94.SSALT99.SSALT14$pval), ]
hist(sort(TemporalSalteryTestResults$ResultsByLocus$SSALT94.SSALT99.SSALT14$pval), breaks = seq(0, 1, 0.05), col = 8); abline(h = 5, col = 2, lwd = 2)
# Should we pool? Stick with 94/99 pooled as per WASSIP, but leave 2014 separate as it was collected much earlier over the entire month of July
# Birch also noted that since these were weir fish, they may have been headed towards the N end of Saltery as "shoal" spawners as opposed to older
# Saltery collections that were collected as part of the egg take near Saltery Creek ("creek" spawners)

source("V:/DATA/R_GEN/JJs GCL/FreqFisPlot4SNPs.GCL.R")
SalteryFreqFis <- FreqFisPlot4SNPs.GCL(sillyvec = c("SSALT94", "SSALT99", "SSALT14", "SLKIT15"), loci = loci96[-mito.loci], groupvec = 1:4, file = "FreqFisPlots for Pooling/SalteryFreqFisPlot.pdf")
str(SalteryFreqFis)
sapply(rownames(SalteryFreqFis$HWE.pval), function(silly) {hist(SalteryFreqFis$HWE.pval[silly, ], breaks = seq(0, 1, 0.05), col = 8, ylim = c(0, 20), main = silly)
                                                           abline(h = 5, col = 2, lwd = 2)} )
sapply(rownames(SalteryFreqFis$HWE.pval), function(silly) {pchisq(q = -2 * sum(log(SalteryFreqFis$HWE.pval[silly, ]), na.rm = TRUE), df = 2 * sum(!is.na(SalteryFreqFis$HWE.pval[silly, ])), lower.tail = FALSE)} )
# Appear to be conform to HWE
table(SSALT14.gcl$attributes$CAPTURE_DATE)



#~~~~~~~~~~~~~~~~~~
## SLKIT15 Regional
KitoiPoolTest <- list(c("SSALT94", "SSALT99", "SSALT14", "SLKIT15"))
KitoiPoolTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = KitoiPoolTest)
KitoiPoolTestResults$OverallResults  # with all SSALT 0; with SSALT94 0.5748; with SSALT99 0.0499; with SSALT14 0

KitoiPoolTestResults$ResultsByLocus$SSALT94.SSALT99.SSALT14[order(TemporalSalteryTestResults$ResultsByLocus$SSALT94.SSALT99.SSALT14$pval), ]

# Do not pool SLKIT15 with Saltery, but still in Saltery RG
sapply(rownames(SalteryFreqFis$HWE.pval), function(silly) {hist(SalteryFreqFis$HWE.pval[silly, ], breaks = seq(0, 1, 0.05), col = 8, ylim = c(0, 20), main = silly)
                                                           abline(h = 5, col = 2, lwd = 2)} )
sapply(rownames(SalteryFreqFis$HWE.pval), function(silly) {pchisq(q = -2 * sum(log(SalteryFreqFis$HWE.pval[silly, ]), na.rm = TRUE), df = 2 * sum(!is.na(SalteryFreqFis$HWE.pval[silly, ])), lower.tail = FALSE)} )
# Check with Birch if these came from Kitoi or they keep getting brood from Saltery. Looks like Saltery late SSALT94 (September) vs. SSALT 15 (end of June).



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### SEAK Baseline Collections
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~
## SDATLAS12 Regional
DatlaskaPoolTest <- list(c("SDATLAS12", "SKWAT11"))
DatlaskaPoolTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = DatlaskaPoolTest)
DatlaskaPoolTestResults$OverallResults  # 0
# Do NOT pool

#~~~~~~~~~~~~~~~~~~
## SNAHL12 Temporal
TemporalNahlinTest <- list(c("SNAHL03", "SNAHL07", "SNAHL12"))
TemporalNahlinTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalNahlinTest)
TemporalNahlinTestResults$OverallResults  # 0.109
# Pool all temporal collections

#~~~~~~~~~~~~~~~~~~
## SSHAKS11 Temporal
TemporalShakesTest <- list(c("SSHAKS06", "SSHAKES07", "SSHAKS09", "SSHAKS11"))
TemporalShakesTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalShakesTest)
TemporalShakesTestResults$OverallResults  # 0; 06, 07, and 09 0.3087; 2011 looks different?
# Should we pool? YES, force it according to SRO 10/28/15, small n

#~~~~~~~~~~~~~~~~~~
## SBIGLA14 Temporal
TemporalBigLakeTest <- list(c("SBIGLK10", "SBIGLA11", "SBIGLA14"))
TemporalBigLakeTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalBigLakeTest)
TemporalBigLakeTestResults$OverallResults  # 0; 10 and 14 0.5849; 10 and 11 0; 11 and 14 0
# Should we pool? YES, force it according to SRO 10/28/15, small n

#~~~~~~~~~~~~~~~~~~
## STHRE10 Temporal
TemporalKlawockTest <- list(c("STHRE04", "SHALF08", "STHRE10"))
TemporalKlawockTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalKlawockTest)
TemporalKlawockTestResults$OverallResults  # 0.0153; 04 and 10 temporal 0.3095
# Should we pool? Put STHRE04 & STHRE10 together, but not SHALF08

source("V:/DATA/R_GEN/JJs GCL/FreqFisPlot4SNPs.GCL.R")
KlawockFreqFis <- FreqFisPlot4SNPs.GCL(sillyvec = c("STHRE04", "SHALF08", "STHRE10"), loci = loci96[-mito.loci], groupvec = 1:3, file = "FreqFisPlots for Pooling/KlawockFreqFisPlot.pdf")
str(KlawockFreqFis)
sapply(rownames(KlawockFreqFis$HWE.pval), function(silly) {hist(KlawockFreqFis$HWE.pval[silly, ], breaks = seq(0, 1, 0.05), col = 8, ylim = c(0, 20), main = silly)
                                                           abline(h = 5, col = 2, lwd = 2)} )
sapply(rownames(KlawockFreqFis$HWE.pval), function(silly) {pchisq(q = -2 * sum(log(KlawockFreqFis$HWE.pval[silly, ]), na.rm = TRUE), df = 2 * sum(!is.na(KlawockFreqFis$HWE.pval[silly, ])), lower.tail = FALSE)} )

#~~~~~~~~~~~~~~~~~~
## SINCK03 & 08
TemporalKlawockCreekTest <- list(c("SINCK03", "SINCK08"))
TemporalKlawockCreekTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalKlawockCreekTest)
TemporalKlawockCreekTestResults$OverallResults  # 0.9977
# Pool all temporal collections

Test <- list(c("SHALF08", "SINCK03", "SINCK08"))
TestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = Test)
TestResults$OverallResults
# Poolt SHALF08 with Inlet Creek collections


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Loose End Baseline Collections
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TemporalNelsonRiverTest <- list(c("SNELSR07", "SHOOD00"))
TemporalNelsonRiverTestResults <- FishersTest.GCL(freq = KMA762Collections96SNPsFreq, loci = loci96, tests = TemporalNelsonRiverTest)
TemporalNelsonRiverTestResults$OverallResults
hist(sort(TemporalNelsonRiverTestResults$ResultsByLocus$SNELSR07.SHOOD00$pval), breaks = seq(0, 1, 0.05), col = 8); abline(h = 5, col = 2, lwd = 2)

# Okay, good thing we don't pool



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Checking HWE in Populations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Check HWE in collections (only use diploid loci)
mito.loci <- dget(file = "Objects/mito.loci.txt")
gcl2Genepop.GCL(sillyvec = KMA473Pops, loci = loci96[-mito.loci], path = "Genepop/KMA473Pops_93nuclearloci.txt", VialNums = TRUE); beep(2)

## Check HWE in Genepop
# Originally done with Exact test and subsequently re-done with MCMC
source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")
HWE <- ReadGenepopHWE.GCL(file = "Genepop/KMA473Pops_93nuclearloci_MCMC.txt.P")
str(HWE)

npops <- length(KMA473Pops)
nloci <- length(loci96[-mito.loci])

# Summary p-values
HWE$SummaryPValues



# Looking into Ayakulik weir samples HWE and Fis
hist(HWE$DataByPop[as.numeric(HWE$DataByPop$Pop) == 187, "WC Fis"], col=8, xlim = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.01))

sapply(which(KMA473PopsGroupVec17 %in% c(5,6,7)), function(pop) {
  hist(HWE$DataByPop[as.numeric(HWE$DataByPop$Pop) == pop, "WC Fis"], 
       col=8, xlim = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.01), main = KMA473Pops[pop])
})




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Overall loci
t(t(HWE$SummaryPValues[1:nloci, "Overall Pops"]))
sort(HWE$SummaryPValues[1:nloci, "Overall Pops"])  # Issues with One_c3-98 (known culprit in SEAK) & One_ACBP-79 (known culprit in Kodiak)

# Markers that are out of HWE with p < 0.01 over all Pops (Fisher's method)
table(HWE$SummaryPValues[1:nloci, "Overall Pops"] < 0.01)  # Number of markers that are p < 0.01 over all pops (Fisher's method)
HWE$SummaryPValues[1:nloci, "Overall Pops"][HWE$SummaryPValues[1:nloci, "Overall Pops"] < 0.01]  # One_c3-98 & One_ACBP-79

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WASSIP approach was to toss out any locus that is p < 0.01 over all Populations (Fisher's method)
KMA473Pops_HWE_locusfail <- loci96[-mito.loci][HWE$SummaryPValues[1:nloci, "Overall Pops"] < 0.01]
dput(x = KMA473Pops_HWE_locusfail, file = "Objects/KMA473Pops_HWE_locusfail.txt")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# How many pops is a marker p < 0.05
sort(apply(HWE$SummaryPValues[1:nloci, 1:npops], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} )) # One_c3-98 out of HWE in 84 / 473 Pops

# Histogram of # pops out of HWP p < 0.05 per loci
par(mar = c(5.1, 5.1, 4.1, 2.1))
hist(apply(HWE$SummaryPValues[1:nloci, 1:npops], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} ), col = 8, xlab = "# Populations p < 0.05", ylab = "# Loci", main = "HWE for 473 Populations", breaks = seq(0, 100, 1), cex.main = 2, cex.lab = 1.5)  # How many populations < 0.05 per marker
abline(v = 0.05 * length(KMA473Pops), col = "red", lwd = 5) # Expectation that it should be p < 0.05 if 5% of Pops
text(x = 0.05 * length(KMA473Pops), y = 8.2, labels = "5% of Populations", col = "red", pos = 4, cex = 1.5)
text(x = unique(sort(apply(HWE$SummaryPValues[1:nloci, 1:npops], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} )))[23:30] - 0.5,
     y = table(sort(apply(HWE$SummaryPValues[1:nloci, 1:npops], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} )))[23:30] + c(-0.1, 0, 0, 0.66, 0.33, 0, 0 , 0),
     labels = sapply(unique(sort(apply(HWE$SummaryPValues[1:nloci, 1:npops], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} )))[23:30], function(npopsout) {paste(names(which(apply(HWE$SummaryPValues[1:nloci, 1:npops], 1, function(locus) {sum(locus < 0.05, na.rm = TRUE)} ) == npopsout)), collapse = ", ")} ),
     srt = 10, pos = 4, cex = 1.2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Overall populations
t(t(HWE$SummaryPValues["Overall Loci", 1:npops]))
sort(HWE$SummaryPValues["Overall Loci", 1:npops])

# Populations that are out of HWE with p < 0.01 over all Pops (Fisher's method)
table(HWE$SummaryPValues["Overall Loci", 1:npops] < 0.01)  # Number of populations that are p < 0.01 over all pops (Fisher's method)
HWE$SummaryPValues["Overall Loci", 1:npops][HWE$SummaryPValues["Overall Loci", 1:npops] < 0.01]  # SLYNXLK09_95    SCHAU01_96   STAKWA11_41   SBIGLA14_95    SMCDO13_70   SCOBB07_101    SBAKE96_97   SCEDAR94_95

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WASSIP approach was to toss out any population that is p < 0.05 for > 10% loci
sort(apply(HWE$SummaryPValues[1:nloci, 1:npops], 2, function(pop) {sum(pop < 0.05, na.rm = TRUE)} ))
KMA473Pops_HWE_popsfail <- KMA473Pops[apply(HWE$SummaryPValues[1:nloci, 1:npops], 2, function(pop) {sum(pop < 0.05, na.rm = TRUE)} ) > (0.1 * nloci)]
dput(x = KMA473Pops_HWE_popsfail, file = "Objects/KMA473Pops_HWE_popsfail.txt")

# KMA473PopsGroupVec15 <- as.numeric(readClipboard())
# dput(x = KMA473PopsGroupVec15, file = "Objects/KMA473PopsGroupVec15.txt")

# Where are these pops?
cbind(KMA473Pops_HWE_popsfail, Groups15[KMA473PopsGroupVec15[KMA473Pops %in% KMA473Pops_HWE_popsfail]])

# NOTE: WE ARE NOT GOING TO TOSS OUT ANY POPULATIONS HERE, DESPITE "RULE BREAKING"
# The decision was made to keep ALL pops from regional baselines, presuming that
# project leaders have already gone through this process and there is no need to
# re-invent the wheel.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


require(xlsx)
write.xlsx(HWE$SummaryPValues, file = "Output/KMA473Pops_HWEmatrix_93nuclearloci_MCMC_Genepop.xlsx")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## View patterns of HWE
# NOTE: white is NA
require(lattice)

# 473 Pops 93 nuclear loci + Overall Pops/Loci
levelplot(t(HWE$SummaryPValues), 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.75, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.75, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 94.25, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 94.25, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          col.regions = c("red", "darkorange", "yellow", "green3", "grey90"), xlab = "SILLY", ylab = "Locus", 
          at = c(0, 0.001, 0.01, 0.05, 1), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares

# 473 Pops 2 HWE FAIL markers
levelplot(t(HWE$SummaryPValues)[, c(loci96[-mito.loci] %in% KMA473Pops_HWE_locusfail, FALSE)], 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.45, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.45, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 2.55, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 2.55, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          col.regions = c("red", "darkorange", "yellow", "green3", "grey90"), xlab = "SILLY", ylab = "Locus", 
          at = c(0, 0.001, 0.01, 0.05, 1), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares



## View distribution of HWE p-values
# By Locus
maxylim <- max(sapply(loci96[-mito.loci], function(loci) {sum(HWE$SummaryPValues[loci, 1:npops] > 0.95, na.rm = TRUE)}))  # Determine ylim
pdf(file = "HWE/HWEPvalueHistograms_byLocus_KMA473Pops_93nuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)   
invisible(sapply(loci96[-mito.loci], function(loci) {hist(HWE$SummaryPValues[loci, 1:npops], col = 8, breaks = seq(from = 0, to = 1, by = 0.05), main = loci, xlab = "P-value", ylim = c(0, maxylim))
                                                     abline(h = npops / 20, col = "red", lwd = 3)} ))
dev.off()

# Only 2 problem markers
maxylim <- max(sapply(loci96[-mito.loci], function(loci) {sum(HWE[[2]][loci, seq(dim(HWE[[2]])[2] - 1)] > 0.95, na.rm = TRUE)}))  # Determine ylim
pdf(file = "HWE/HWEPvalueHistograms_byLocus_KMA473Pops_2failnuclearloci_MCMC_Genepop.pdf",width=11, height=8.5, family="Times",pointsize=20)   
invisible(sapply(KMA473Pops_HWE_locusfail, function(loci) {hist(HWE$SummaryPValues[loci, 1:npops], col = 8, breaks = seq(from = 0, to = 1, by = 0.05), main = loci, xlab = "P-value", ylim = c(0, maxylim))
                                                           abline(h = npops / 20, col = "red", lwd = 3)} ))
dev.off()


# By Pops
maxylim <- max(sapply(seq(npops), function(silly) {sum(HWE$SummaryPValues[1:nloci, silly] > 0.95, na.rm = TRUE)}))  # Determine ylim
pdf(file = "HWE/HWEPvalueHistograms_byPop_KMA473Pops_93nuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
invisible(sapply(seq(npops), function(silly) {hist(HWE$SummaryPValues[1:nloci, silly], col = 8, breaks = seq(from = 0, to = 1, by = 0.05), main = KMA473Pops[silly], xlab = "P-value", ylim = c(0, maxylim))
                                              abline(h = nloci * 0.05, col = "red", lwd = 3)} ))
dev.off()

# Only 17 failed Pops
KMA473Pops_HWE_popsfail_number <- sapply(KMA473Pops_HWE_popsfail, function(pop) {which(KMA473Pops == pop)} )

maxylim <- max(sapply(KMA473Pops_HWE_popsfail_number, function(silly) {sum(HWE$SummaryPValues[1:nloci, silly] > 0.95, na.rm = TRUE)}))  # Determine ylim
pdf(file = "HWE/HWEPvalueHistograms_byPop_KMA17PopsFail_93nuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
invisible(sapply(KMA473Pops_HWE_popsfail_number, function(silly) {hist(HWE$SummaryPValues[1:nloci, silly], col = 8, breaks = seq(from = 0, to = 1, by = 0.05), main = KMA473Pops[silly], xlab = "P-value", ylim = c(0, maxylim))
                                                                  abline(h = nloci * 0.05, col = "red", lwd = 3)} ))
dev.off()
# There are a few populations that look suspicious, but I will lean on previous baseline work and keep them in
# No need to re-invent the wheel

# Where are these 16 pops from, just out of curiousity
setNames(object = Groups15[KMA473PopsGroupVec15[KMA473Pops_HWE_popsfail_number]], nm = KMA473Pops_HWE_popsfail)



## View Fis vs. Pvalue for all markers
# All 93 nuclear markers
pdf(file = "HWE/HWEPvalueVsFis_byPop_KMA473Pops_93nuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
invisible(sapply(loci96[-mito.loci], function(locus) {plot(x = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"], y = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"], pch = 16, xlim = c(-1, 1), ylim = c(0, 1), xlab = "WC Fis", ylab = "P-value", main = locus, col = Colors15[KMA473PopsGroupVec15], cex = 0.7)
                                                      legend("topright", legend = Groups15Short, fill = Colors15, xpd = TRUE, bty = "n", cex = 0.7)
                                                      abline(h = 0.05, col = "red", lwd = 3)
                                                      text(x = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"][which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05)], y = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"][which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05)],
                                                           labels = which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05), cex = 0.4)} ))
dev.off()


# Only 2 problem markers
pdf(file = "HWE/HWEPvalueVsFis_byPop_KMA473Pops_2failnuclearloci_MCMC_Genepop.pdf", width = 11, height = 8.5, family = "Times", pointsize = 20)
invisible(sapply(KMA473Pops_HWE_locusfail, function(locus) {plot(x = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"], y = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"], pch = 16, xlim = c(-1, 1), ylim = c(0, 1), xlab = "WC Fis", ylab = "P-value", main = locus, col = Colors15[KMA473PopsGroupVec15], cex = 0.7)
                                                            legend("topright", legend = Groups15Short, fill = Colors15, xpd = TRUE, bty = "n", cex = 0.7)
                                                            abline(h = 0.05, col = "red", lwd = 3)
                                                            text(x = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "WC Fis"][which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05)], y = HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"][which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05)], 
                                                                 labels = which(HWE$DataByPop[which(HWE$DataByPop$Locus == locus), "PValue"] < 0.05), cex = 0.4)} ))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Deciding what Loci to drop for HWE ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### NOTE: THIS SECTION BELOW IS OLD AND PRESERVED HERE FOR ARCHIVAL PURPOSES, THIS WAS WRITTEN USING EXACT TEST RESULTS AS OPPOSED TO MCMC RESULTS FOR TESTING HWE IN GENEPOP
### ONLY 2 MARKERS ARE DROPPED DUE TO NON-CONFORMANCE TO HWP ("One_ACBP-79" "One_c3-98")


#??????????????????????????????????????????????????????????????????????????????
# Decided against the WASSIP approach of Fisher's p < 0.01 since there were 11 loci out, inluding 2 of the most variable (important for MSA) markers.
# The big problem with Fisher's method is that if there are hard 0's, then Fisher's Chisqr goes to infinity (p-value = 0).
# Used a combination of looking at the geographic patterns of HWP p-values, histograms of those p-values, and plots of HWP p-values vs. Fis values.
# Several of these markers had 1 or a few pops with very high + Fis values, indicating that these pops are likely almost fixed and a heterozygote got mis-scored as an alternate homozygote.
# i.e. all blues and one red (no greens)
# On 7/27/15 Tyler suggested that I was on the right track and to come up with something defensible
# I went through all 11 markers one at a time with Chris and we made the following decisions.


# One_ACBP-79:    DROP. Weirdness in Kodiak, specifically with Ayakulik/Frazer, which will be one of the toughest RGs to pull apart. Fis was +.
# One_aldB-152:   KEEP. SHARRIET12 [220] has high Fis, likely due to reason listed above (mostly red homo, one het, and one blue homo; chip looks good, just a weird result).  # table(SHARRIET12.gcl$counts[, "One_aldB-152", 1])
# One_c3-98:      DROP. Fis is highly skewed -. LOTS of pops are out of HWE for PWS-SEAK. No red homozygotes. Lots of overcalling heterozygotes!!! Likely a duplicated locus???  # table(SWIND03.SWIND07.gcl$counts[, "One_c3-98", 1]); table(SKSLK10.SKSLK11.gcl$counts[, "One_c3-98", 1])
# One_E2:         KEEP. SESHAR08.SESHA91 [296] has extremely high Fis, likely due to reason listed above.  # table(SESHAR08.SESHA91.gcl$counts[, "One_E2", 1])
# One_IL8r-362:   KEEP. SEEK04.SEEK07 [431] has extremely high Fis, likely due to reason listed above.  # table(SEEK04.SEEK07.gcl$counts[, "One_lpp1-44", 1]); low sample size...; table(SKLAG09.gcl$counts[, "One_lpp1-44", 1])
# One_lpp1-44:    KEEP. SMAHO03.SMAHO07 [424] has extremely high Fis, likely due to reason listed above.  # table(SMAHO03.SMAHO07.gcl$counts[, "One_lpp1-44", 1])
# One_MHC2_190:   KEEP. SNADE95 [460] and SMISS98 [63] have extremely high Fis, likely due to reason listed above. No heterozygotes!  # table(SNADE95.gcl$counts[, "One_MHC2_190", 1]); table(SMISS98.gcl$counts[, "One_MHC2_190", 1])
# One_sast-211:   DROP. Fis is highly skewed +. Also, not a very informative marker so it will likely get dropped anyways. High Fis due to reason listed above.  # table(SONGI06.gcl$counts[, "One_sast-211", 1]); table(SMCDON06.gcl$counts[, "One_sast-211", 1])
# One_txnip-401:  KEEP. SMUDA05 [123] has extremely high Fis, likely due to reason listed above.  # table(SMUDA05.gcl$counts[, "One_txnip-401", 1])
# One_U1004-183:  KEEP. Fis is normally distributed with some + and some -. This pop may contain an Early/Late run SNEVA09.SNEVA13  # table(SNEVA09.SNEVA13.gcl$counts[, "One_U1004-183", 1])
# One_U1012-68:   KEEP. SANVI06 [46] has high Fis, likely due to reason listed above.  # table(SANVI06.gcl$counts[, "One_U1012-68", 1])

# KMA473Pops_HWE_locusfail[c(1, 3, 8)]  # "One_ACBP-79"  "One_c3-98"    "One_sast-211"
# KMA473Pops_HWE_locusfail_remove <- KMA473Pops_HWE_locusfail[c(1, 3, 8)]



# # Worth re-running HWE in Genepop with MCMC as opposed to Exact test. Then convert those hard 0's to "repzero <- format(1 / (batches * iterations), scientific = FALSE, digits = 6)"
# source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopHWEKS.GCL.R")
# HWEv2 <- ReadGenepopHWEKS.GCL(file = "Genepop/KMA473Pops_93nuclearloci_MCMC.txt.P")
# # source("V:/WORK/Kyle/R Scripts/OLD/ReadGenepopHWEKS_OLD.GCL.R")
# # HWEcomp <- ReadGenepopHWEKS.GCL(file = "Genepop/KMA473Pops_93nuclearloci.txt.P")
# 
# # New function works for both Exact and MCMC
# # USE MCMC for actual project!!!!!
# rownames(HWEcomp[[2]])[which(HWEcomp[[2]][, length(colnames(HWEcomp[[2]]))] < 0.01)]  # "One_ACBP-79"   "One_aldB-152"  "One_c3-98"     "One_E2"        "One_IL8r-362"  "One_lpp1-44"   "One_MHC2_190"  "One_sast-211"  "One_txnip-401" "One_U1004-183" "One_U1012-68"
# rownames(HWE[[2]])[which(HWE[[2]][, length(colnames(HWE[[2]]))] < 0.01)]  # "One_ACBP-79"   "One_aldB-152"  "One_c3-98"     "One_E2"        "One_IL8r-362"  "One_lpp1-44"   "One_MHC2_190"  "One_sast-211"  "One_txnip-401" "One_U1004-183" "One_U1012-68"
# rownames(HWEv2[[2]])[which(HWEv2[[2]][, length(colnames(HWEv2[[2]]))] < 0.01)]  # "One_ACBP-79"  "One_c3-98"
#??????????????????????????????????????????????????????????????????????????????

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Toss out markers not conforming to HWP ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# See above for detailed explanation
# Briefly, One_ACBP-79 has +Fis in Kodiak (Ayakulik/Frazer) & One_c3-98 has very skewed -Fis especially in SEAK
KMA473Pops_HWE_locusfail
dput(x = KMA473Pops_HWE_locusfail, file = "Objects/KMA473Pops_HWE_locusfail_remove.txt")

loci94 <- loci96[!loci96 %in% KMA473Pops_HWE_locusfail]
dput(x = loci94, file = "Objects/loci94.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace; dget .gcl objects and Locus Control ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline")

rm(list = ls(all = TRUE))
# This sources all of the new GCL functions to this workspace
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")
# This sources all of Kyle's non-GCL useful functions to this workspace (does not include "tempGCL" functions)
source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")


## Get objects
baselineobjects <- sapply(list.files(path = "Objects"), function(file) {unlist(strsplit(x = file, split = ".txt"))} )
invisible(sapply(baselineobjects, function(file) {assign(x = file, value = dget(file = paste("Objects/", file, ".txt", sep = "")), pos = 1)} ))
LocusControl <- OriginalLocusControl; rm(OriginalLocusControl)

## Get Populations

require(beepr)
invisible(sapply(KMA473Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPops/", silly, ".txt", sep = "")), pos = 1)} )); beep(8)
objects(pattern = "\\.gcl")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Allele Frequency Plots Populations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### All 473 Pops
# Counts
KMA473PopsPostRemovalsAlleleCounts <- FreqPop.GCL(sillyvec = KMA473Pops, loci = loci96)
str(KMA473PopsPostRemovalsAlleleCounts)

# Frequencies
KMA473PopsPostRemovalsFreqs <- KMA473PopsPostRemovalsAlleleCounts[,,"Allele 1"] / (KMA473PopsPostRemovalsAlleleCounts[,,"Allele 2"] + KMA473PopsPostRemovalsAlleleCounts[,,"Allele 1"])
dput(x = KMA473PopsPostRemovalsFreqs, file = "Objects/KMA473PopsPostRemovalsFreqs.txt")
str(KMA473PopsPostRemovalsFreqs)
write.xlsx(KMA473PopsPostRemovalsFreqs, file = "Output/KMA473PopsPostRemovalsFreqs.xlsx")

#### Create plots
## Create groupvec
KMA473PopsGroupVec15


## Plot points and smoothing line
getwd()
pdf(file = "FreqPlots/KMA473Pops_96SNPs_GroupVec15_FreqPlots.pdf",width=11, height=8.5, family="Times",pointsize=20)   
par(mar = c(2.1, 4.1, 4.1, 4.6))
for(locus in loci96){
  plot(KMA473PopsPostRemovalsFreqs[, locus], main = locus, col = Colors15[KMA473PopsGroupVec15], pch = 19, ylim = c(0, 1), ylab = "Frequency", xlab = "Collections", xaxt = 'n', cex = 0.7)
  lines(supsmu(seq(length(KMA473Pops)), KMA473PopsPostRemovalsFreqs[, locus]), lwd = 2)
  mtext(text = "Collections", side = 1, line = 0.7)
  text(x = seq(KMA473Pops), y = KMA473PopsPostRemovalsFreqs[, locus], labels = seq(KMA473Pops), cex = 0.2)
  legend(x = 490, y = 0.9, legend = Groups15Short, fill = Colors15, xpd = TRUE, bty = "n", cex = 0.7)
}; rm(locus)
dev.off()

## View levelplot of frequencies
require(lattice)
new.colors <- colorRampPalette(c("black", "white"))
levelplot(t(KMA473PopsPostRemovalsFreqs), col.regions = new.colors, xlab = "Locus", ylab = "SILLY", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares
levelplot(KMA473PopsPostRemovalsFreqs, col.regions = new.colors, xlab = "SILLY", ylab = "Locus", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MDS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dump a GENEPOP file for genind
require(adegenet)
gcl2Genepop.GCL(sillyvec = KMA473Pops, loci = loci94, path = "Genepop/KMA473Pops_94loci.gen")
genind <- read.genepop(file = "Genepop/KMA473Pops_94loci.gen")

genpop <- genind2genpop(genind)

AdegenetNei473Pop94loci <- dist.genpop(genpop, method = 1, diag = TRUE, upper = TRUE)
dput(x = AdegenetNei473Pop94loci, file = "Trees/AdegenetNei473Pop94loci.txt")
str(AdegenetNei473Pop94loci)

require(ape)
Nei473NJtree <- nj(AdegenetNei473Pop94loci)
str(Nei473NJtree)
plot.phylo(x = Nei473NJtree, cex = 0.5, no.margin = TRUE)


## For FigTree
# CommonNames473 <- readClipboard()
# dput(x = CommonNames473, file = "Objects/CommonNames473.txt")

CommonNames473 <- dget(file = "Objects/CommonNames473.txt")
Nei473NJtree$tip.label <- CommonNames473
plot(Nei473NJtree, font = 2, cex = 0.5)
axis(1)  # Adds scale to bottom of plot 

write.tree(Nei473NJtree, "GDA/NJofNei473Popstree.nex", tree.names = CommonNames473)

# TahltanLakeTemporalResults <- FishersTest.GCL(freq = KMA473PopsPostRemovalsAlleleCounts, loci = loci96, tests = list(KMA473Pops[grep(pattern = "Tahltan Lake", x = CommonNames473)]))
# TahltanLakeTemporalResults
# TahltanLakeTemporalResults$OverallResults
# hist(TahltanLakeTemporalResults[[2]][[1]][seq(loci96), "pval"]) 
# # Correct to not Pool, just need to change CommonNames473

## MDS
library('rgl')

# MDS <- cmdscale(as.matrix(AdegenetNei473Pop94loci), k = 3)  # Did in base R as it kept crashing in RStudio...dunno why.
# dput(x = MDS, file = "Objects/MDSAdegenetNei473Pop94loci.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei473Pop94loci.txt")

x <- as.vector(MDS[, 1])   
y <- as.vector(MDS[, 2])
z <- as.vector(MDS[, 3])

KMA473PopsGroupVec15 <- dget(file = "Objects/KMA473PopsGroupVec15.txt")
Colors15 <- dget(file = "Objects/Colors15.txt")

par(family = "serif")
plot3d(x, y, z + abs(range(z)[1]), xlab = '', ylab = '', zlab = '', aspect = FALSE, col = Colors15[KMA473PopsGroupVec15], size = 0.5, type = 's', axes = TRUE, box = TRUE, top = TRUE, cex = 1)
box3d()
# plot3d(x, y, z + abs(range(z)[1]), aspect = FALSE, col = "black", size = 0.5, type = 'h', box = TRUE, axes = FALSE, top = FALSE, add = TRUE)  # adds pins to spheres
texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.2, 0), text = seq(KMA473Pops), font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)

dir.create("MDS")
rgl.snapshot("MDS/MDSAdegenetNei473Pop94loci.png", fmt="png", top=TRUE )


# Plot individually
plot(x, col = Colors15[KMA473PopsGroupVec15], pch = 16)
plot(y, col = Colors15[KMA473PopsGroupVec15], pch = 16)
CommonNames473[which.min(y)]; CommonNames473[which.max(y)]
plot(z, col = Colors15[KMA473PopsGroupVec15], pch = 16)
CommonNames473[which.min(z)]; CommonNames473[which.max(z)]

CommonNames473[302:320]  # Part of Copper is different

# Outliers? Kanalku Lake, King Salmon Lake, Kaflia, Upper Tlikakila River (Lake Clark), Mendeltna Creek (Copper), and Kitlope Lake (BC/WA)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Linkage Disequilibrium ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## It took two weeks to to LD on all 473 Pops
# AND it only did the first 200 or so pops!
# New strategy is to make several Genepop files (.gen) and run LD in batches of ~ 50 Pops or so
# It should be simple to combine them on the back end


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Go by 50 Pops at a time...
# Creating input genepop files
loci93nuclear <- loci96[-mito.loci]
invisible(sapply(1:9, function(i) {
  gcl2Genepop.GCL(sillyvec = KMA473Pops[(((i - 1) * 50) + 1):(i * 50)], loci = loci93nuclear, path = paste("Genepop/", i, "KMA473Pops_", (((i - 1) * 50) + 1), "_to_", (i * 50), "_93nuclearloci.txt", sep = ""), VialNums = TRUE)
} )); beep(8)

# Last set of 14 Pops
gcl2Genepop.GCL(sillyvec = KMA473Pops[451:473], loci = loci93nuclear, path ="Genepop/10KMA473Pops_451_to_473_93nuclearloci.txt", VialNums = TRUE)

# Number of pairwise comparisons
sum(seq(length(loci93nuclear) - 1)) * length(KMA473Pops)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("V:/DATA/R_GEN/tempGCL Source Scripts/ReadGenepopDisKS.GCL.R")

# Get file paths
LD.files <- grep(pattern = "KMA473Pops", x = list.files(path = "Genepop", pattern = ".DIS", full.names = TRUE), value = TRUE)
LD.files

# Read in .DIS files
LD <- lapply(1:10, function(i) {ReadGenepopDisKS.GCL(file = grep(pattern = paste(i, "KMA473Pops", sep = ""), x = LD.files, value = TRUE))} ); beep(8)
str(LD, max.level = 1)
str(LD[[1]])
table(LD[[1]]$Locus1 == LD[[2]]$Locus1)

# Combine mutliple dataframes in list into a single dataframe
LD.all <- cbind(Locus1 = LD[[1]]$Locus1, Locus2 = LD[[1]]$Locus2, do.call(what = "cbind", lapply(LD, function(ld) {ld[, 3:(dim(ld)[2] - 1)]} )))
str(LD.all); dim(LD.all)  # dim should be = npops + 2 (Locus1, Locus2)
dput(x = LD.all, file = "Objects/KMA473Pops_93nuclearlociLD.txt")

LD.all <- dget(file = "Objects/KMA473Pops_93nuclearlociLD.txt")

str(LD.all, max.level = 1)
dimnames(LD.all)
head(LD.all)
tail(LD.all)

# Number of pops p < 0.05 for all pairs
LD.all$nPopsFail <- apply(LD.all[, 3:(length(KMA473Pops) + 2)] < 0.05, 1, function(locuspair) {sum(locuspair)} )

# Which pairs have 1/2 of all pops at p < 0.05
LD.all$Fail <- LD.all$nPopsFail > (length(KMA473Pops) / 2)

LD.all[which(LD.all$Fail == TRUE), ]

# Only two possible pairs of loci exhibit LD in > half of pops (> half being 26 or more of 50 pops), these three are known culprits.
# Marker pairs that fail are: 1) One_GPDH2 & One_GPDH & 2)One_MHC2_251 & One_MHC2_190. Both of these were linked for WASSIP.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create LD Histogram Plot for MS
LD.all[which(LD.all$nPopsFail > 60), c(1:2, 476:477)]

library('gplots')
library('devEMF')
emf(file = "LD/Histogram of proportion of pops in LD.emf", width = 7.5, height = 6, family = "Times")
par(mar = c(5.1, 4.1, 1.1, 1.1))
hist(LD.all$nPopsFail / length(KMA473Pops), breaks = 50, col = 8, xlab = expression(paste("Proportion of populations ", italic("P"), " < 0.05", sep = '')), main = "", cex.lab = 1.2, cex.axis = 1.2); table(LD.all$nPopsFail)
text(x = LD.all[435, "nPopsFail"] / length(KMA473Pops) - 0.01, y = 200, labels = expression(paste(italic("One_MHC2_190_251"))), srt = 90, pos = 4, col = "blue", cex = 1.2)
text(x = LD.all[120, "nPopsFail"] / length(KMA473Pops) - 0.01, y = 200, labels = expression(paste(italic("One_GPDH-201_GPDH2-187"))), srt = 90, pos = 4, col = "darkgreen", cex = 1.2)
text(x = LD.all[1540, "nPopsFail"] / length(KMA473Pops) - 0.01, y = 200, labels = expression(paste(italic("One_Tf_ex11-750_in3-182"))), srt = 90, pos = 4, col = "black", cex = 1.2)
arrows(x0 = LD.all[c(435, 120, 1540), "nPopsFail"] / length(KMA473Pops), y0 = rep(190, 3), x1 = LD.all[c(435, 120, 1540), "nPopsFail"] / length(KMA473Pops), y1 = rep(10, 3), lwd = 3, length = 0.15)
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hist(LD.all$nPopsFail / length(KMA473Pops), breaks = 50, col = 8, xlab = expression(paste("Proportion of populations ", italic(P), " < 0.05", sep = '')), ylim = c(0, 10))
hist(LD.all$nPopsFail / length(KMA473Pops), breaks = seq(0, 0.7, 0.005), col = 8, xlab = expression(paste("Proportion of populations ", italic(P), " < 0.05", sep = '')), main = "", ylim = c(0, 10))

hist(LD.all$nPopsFail, breaks = 50, col = 8, ylim = c(0,10))
hist(LD.all$nPopsFail, breaks = 300, col = 8, ylim = c(0,10))

LD.all[which(LD.all$nPopsFail > 60), c(1:2, 476:477)]

# 3) One_Tf_ex3-182 & One_Tf_ex10-750 is also far enough outside the distribution to be considered linked

# No other marker pairs appear to be problematic ("substaintially more" linked than other pairs as defined in WASSIP, i.e. not outside of distribution/histogram)
# Even though there is a cluster of 6 pairs of loci that are linked for ~15% of pops (outside the tail of the distribution), we did not consider these linked in WASSIP and won't here

## Looking at all LD pairs outside of the distribution (i.e. any gap at all)
# See if they are in any of my marker sets
LD.all[which(LD.all$nPopsFail > 60), c(1:2, 476:477)]

# LD.all[which(LD.all$nPopsFail > 60), 1:2][apply(apply(LD.all[which(LD.all$nPopsFail > 60), 1:2], 1, function(LDpair) {as.character(LDpair) %in% loci22} ), 2, function(pr) sum(pr) == 2), ]
# LD.all[which(LD.all$nPopsFail > 60), 1:2][apply(apply(LD.all[which(LD.all$nPopsFail > 60), 1:2], 1, function(LDpair) {as.character(LDpair) %in% loci24} ), 2, function(pr) sum(pr) == 2), ]
# LD.all[which(LD.all$nPopsFail > 60), 1:2][apply(apply(LD.all[which(LD.all$nPopsFail > 60), 1:2], 1, function(LDpair) {as.character(LDpair) %in% loci24Kodiak} ), 2, function(pr) sum(pr) == 2), ]
# LD.all[which(LD.all$nPopsFail > 60), 1:2][apply(apply(LD.all[which(LD.all$nPopsFail > 60), 1:2], 1, function(LDpair) {as.character(LDpair) %in% loci46} ), 2, function(pr) sum(pr) == 2), ]


# Use levelplot to see geographic patterns of LD
require(lattice)
new.colors <- colorRampPalette(c("white", "black"))

levelplot(t(data.matrix(LD.all[which(LD.all$nPopsFail > 60), 3:(2 + length(KMA473Pops))])), 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.45, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.45, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 9.55, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 9.55, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          col.regions = c("red", "darkorange", "yellow", "green3", "grey90"), xlab = "SILLY", ylab = "Locus Pair", 
          at = c(0, 0.001, 0.01, 0.05, 1), scales = list(x = list(rot = 90)), aspect = "fill")  # aspect = "iso" will make squares


#### FINAL DECISION: only consider GDPH, MHC, and Tf_ex pairs as linked

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### fORCA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## fORCA, run in separate R script 500 locus pairs at a time, so I can move on with other things
# KMA473PopsFullfORCA <- RandomPairfORCA.GCL(sillyvec = KMA473Pops, loci = loci96[-mito.loci], groupvec = KMA473PopsGroupVec15, nblocks = 4278, NREP = 1000) # Run over weekend?

# ## dget fORCA 1-9 and put together into one object
# # Get file paths
# fORCA.files <- list.files(path = "fORCA", pattern = ".txt", full.names = TRUE)
# fORCA.files
# 
# # Read in .txt files
# KMA473PopsFullfORCA <- lapply(fORCA.files, function(i) {dget(file = paste(getwd(), i, sep = "/"))} ); beep(2)
# str(KMA473PopsFullfORCA)
# 
# KMA473PopsFullfORCA <- do.call(what = "rbind", args = KMA473PopsFullfORCA)
# str(KMA473PopsFullfORCA)  # n obs. = choose(n = nloci, k = 2) where nloci = 93
# head(KMA473PopsFullfORCA)
# dput(x = KMA473PopsFullfORCA, file = "Objects/KMA473PopsFullfORCA.txt")

KMA473PopsFullfORCA <- dget(file = "Objects/KMA473PopsFullfORCA.txt")
str(KMA473PopsFullfORCA)

# Calculate delta fORCA
deltafORCA <- KMA473PopsFullfORCA$Locus1.Locus2 - pmax(KMA473PopsFullfORCA$Locus1, KMA473PopsFullfORCA$Locus2)
names(deltafORCA) <- rownames(KMA473PopsFullfORCA)
dput(x = deltafORCA, file = "Objects/deltafORCA.txt")
deltafORCA <- dget(file = "Objects/deltafORCA.txt")

# View distribution of delta fORCA
hist(deltafORCA, breaks = 100, col = 8, xlab = expression(paste(Delta, "fORCA", sep = "")), main = expression(paste("Full ", Delta, "fORCA", sep = "")), xlim = c(-0.05, 0.07), cex.main = 2, cex.axis = 1.5, cex.lab = 1.5)
mtext(text = "nblocks = all 4278 possible combinations", side = 3, line = 0, cex = 1.5)

# Determine and plot delta90
delta90 <- quantile(x = deltafORCA, probs = 0.9)
abline(v = delta90, col = "red", lwd = 5)  # delta90 = 0.022

delta90fORCA <- names(deltafORCA[which(deltafORCA >= delta90)])
delta90fORCA[intersect(grep("MHC2_251", delta90fORCA), grep("MHC2_190", delta90fORCA))]
delta90fORCA[intersect(grep("GPDH", delta90fORCA), grep("GPDH2", delta90fORCA))]
delta90fORCA[intersect(grep("Tf_ex10", delta90fORCA), grep("Tf_ex3", delta90fORCA))]


## Add lines for three identified linked pairs
abline(v = deltafORCA[intersect(grep("MHC2_251", names(deltafORCA)), grep("MHC2_190", names(deltafORCA)))], col = 4, lwd = 5)  # 0.007
abline(v = deltafORCA[grep("GPDH.One_GPDH", names(deltafORCA))], col = "darkgreen", lwd = 5)  # 0.005
abline(v = deltafORCA[intersect(grep("Tf_ex10", names(deltafORCA)), grep("Tf_ex3", names(deltafORCA)))], lwd = 5)  # 0.061

legend("topleft", legend = c(expression(Delta[90]), "MHC2_190_251", "GPDH_GDPH2", "Tf_ex10_ex3"), text.col = c("red", "blue", "darkgreen", "black"), bty = 'n', cex = 1.5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create fORCA Histogram Plot for MS
library('gplots')
library('devEMF')
emf(file = "LD/Histogram of deltafORCA.emf", width = 7.5, height = 6, family = "Times")
par(mar = c(5.1, 4.1, 1.1, 1.1))
hist(deltafORCA, breaks = 100, col = 8, xlab = expression(paste("Pair ", italic("f")["ORCA"], " - max(single ", italic("f")["ORCA"], ")", sep = "")), main = '', xlim = c(-0.07, 0.07), cex.main = 2, cex.axis = 1.4, cex.lab = 1.4)
legend("topleft", legend = c(expression(Delta[90]), expression(italic("One_MHC2_190_251")), expression(italic("One_GPDH-201_GPDH2-187")), expression(italic("One_Tf_ex11-750_in3-182"))), text.col = c("red", "blue", "darkgreen", "black"), bty = 'n', cex = 1.2)

abline(h = 0)

segments(x0 = delta90, y0 = 130, x1 = delta90, y1 = 0, col = "red", lwd = 4)
segments(x0 = deltafORCA[intersect(grep("MHC2_251", names(deltafORCA)), grep("MHC2_190", names(deltafORCA)))], y0 = 130, x1 = deltafORCA[intersect(grep("MHC2_251", names(deltafORCA)), grep("MHC2_190", names(deltafORCA)))], y1 = 0, col = "blue", lwd = 4)
segments(x0 = deltafORCA[grep("GPDH.One_GPDH", names(deltafORCA))], y0 = 130, x1 = deltafORCA[grep("GPDH.One_GPDH", names(deltafORCA))], y1 = 0, col = "darkgreen", lwd = 4)
segments(x0 = deltafORCA[intersect(grep("Tf_ex10", names(deltafORCA)), grep("Tf_ex3", names(deltafORCA)))], y0 = 130, x1 = deltafORCA[intersect(grep("Tf_ex10", names(deltafORCA)), grep("Tf_ex3", names(deltafORCA)))], y1 = 0, col = "black", lwd = 4)

# abline(v = delta90, col = "red", lwd = 4)  # delta90 = 0.022
# abline(v = deltafORCA[intersect(grep("MHC2_251", names(deltafORCA)), grep("MHC2_190", names(deltafORCA)))], col = 4, lwd = 4)  # 0.007
# abline(v = deltafORCA[grep("GPDH.One_GPDH", names(deltafORCA))], col = "darkgreen", lwd = 4)  # 0.005
# abline(v = deltafORCA[intersect(grep("Tf_ex10", names(deltafORCA)), grep("Tf_ex3", names(deltafORCA)))], lwd = 4)  # 0.061

dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






## Need to check with Jim regarding fORCA vs. LinkageCorrelationsJAGS
delta90
KMA473PopsFullfORCA[intersect(grep("MHC2_251", names(deltafORCA)), grep("MHC2_190", names(deltafORCA))), ] # Combine > MHC2_190 > MHC2_251 (but not enough to combine)
KMA473PopsFullfORCA[grep("GPDH.One_GPDH", names(deltafORCA)), ]  # Combine > GPDH > GPDH2 (but not enough to combine)
KMA473PopsFullfORCA[intersect(grep("Tf_ex10", names(deltafORCA)), grep("Tf_ex3", names(deltafORCA))), ]  # Combine > Tf_ex3 > Tf_ex10 (enough to combine)


# Which single loci are the best for fORCA?
t(t(sort(setNames(object = c(KMA473PopsFullfORCA[1, 1], KMA473PopsFullfORCA[1:92, 2]), nm = loci96[-mito.loci]))))


t(t(sort(setNames(object = c(KMA473PopsFullfORCA[1849, 1], KMA473PopsFullfORCA[1849:1940, 2]), nm = loci96[-mito.loci]))))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### LinkageCorrelation.GCL ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dir.create("BRugs")
# require(BRugs)
# 
# KMA473PopsLinkageCorrelationMHC <- LinkageCorrelation.GCL(sillyvec = KMA473Pops, markerset = grep(pattern = "MHC", x = loci96, value = TRUE), groupvec = KMA473PopsGroupVec15, groupnames = Groups15Short, dir = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/BRugs"); beep(8)
# dput(x = KMA473PopsLinkageCorrelationMHC, file = "Objects/KMA473PopsLinkageCorrelationMHC.txt")
# 
# rm(list = objects(pattern = "temp.gcl"))  # remove temp.gcl objects
# 
# KMA473PopsLinkageCorrelationGPDH <- LinkageCorrelation.GCL(sillyvec = KMA473Pops, markerset = grep(pattern = "GPDH", x = loci96, value = TRUE), groupvec = KMA473PopsGroupVec15, groupnames = Groups15Short, dir = "BRugs")
# dput(x = KMA473PopsLinkageCorrelationGPDH, file = "Objects/KMA473PopsLinkageCorrelationGPDH.txt")
# 
# KMA473PopsLinkageCorrelationTf_ex <- LinkageCorrelation.GCL(sillyvec = KMA473Pops, markerset = grep(pattern = "Tf_ex", x = loci96, value = TRUE), groupvec = KMA473PopsGroupVec15, groupnames = Groups15Short, dir = "BRugs")
# dput(x = KMA473PopsLinkageCorrelationTf_ex, file = "Objects/KMA473PopsLinkageCorrelationTf_ex.txt")
# 
# # These didn't work!!! Going to use the "newer" LinkageCorrelationsJAGS.GCL and modify to get population level correlation coefficients

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### LinkageCorrelationsJAGS.GCL ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Jim's "new" function to check and see if should combine or not
# Testing for known culprits in Sockeye2011_96SNPs, that were subsequently confirmed with Genepop.DIS, see "LD", only these three pairs failed at > half of pops

## MHC
KMA473PopsLinkageCorMHC <- LinkageCorrelationJAGS.GCL(sillyvec = KMA473Pops, markerset = c("One_MHC2_190", "One_MHC2_251"), groupvec = KMA473PopsGroupVec15, groupnames = Groups15Short); beep(8)
KMA473PopsLinkageCorMHC  # MHC2_190 > Combine > MHC2_251; Drop MHC2_251
dput(x = KMA473PopsLinkageCorMHC, file = "Objects/KMA473PopsLinkageCorMHC.txt")

## GDPH
KMA473PopsLinkageCorGPDH <- LinkageCorrelationJAGS.GCL(sillyvec = KMA473Pops, markerset = c("One_GPDH", "One_GPDH2"), groupvec = KMA473PopsGroupVec15, groupnames = Groups15Short); beep(8)
KMA473PopsLinkageCorGPDH  # GPDH > GPDH2 > Combine; Drop GDPH2, barely, both are very similar in performance (GDPH is better in Chignik and Kodiak for both 11 and 9 RGs and for WASSIP)
dput(x = KMA473PopsLinkageCorGPDH, file = "Objects/KMA473PopsLinkageCorGPDH.txt")

## Tf_ex
KMA473PopsLinkageCorTf_ex <- LinkageCorrelationJAGS.GCL(sillyvec = KMA473Pops, markerset = c("One_Tf_ex10-750", "One_Tf_ex3-182"), groupvec = KMA473PopsGroupVec15, groupnames = Groups15Short); beep(8)  # 
KMA473PopsLinkageCorTf_ex  # Tf_ex3 > Combine > Tf_ex10; Drop Tf_ex10-750 (Same as WASSIP)
dput(x = KMA473PopsLinkageCorTf_ex, file = "Objects/KMA473PopsLinkageCorTf_ex.txt")

## This method says to drop and not combine any of the three linkage groups. The choice of dropping is similar for both Kodiak and Chignik baselines
## This differs from the results for Bristol Bay, which combined MHC, but dropped GDPH2 and Tf_ex10-750
## In SEAK, Serena ended up dropping MHC2_251 (same as here) but dropping GDPH instead of GPDH2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## View Linkage Correlations
require(gplots)
# MHC
barplot2(height = KMA473PopsLinkageCorMHC[5:dim(KMA473PopsLinkageCorMHC)[1], "mean"], plot.ci = TRUE,
         ci.l = KMA473PopsLinkageCorMHC[5:dim(KMA473PopsLinkageCorMHC)[1], "2.5%"],
         ci.u = KMA473PopsLinkageCorMHC[5:dim(KMA473PopsLinkageCorMHC)[1], "97.5%"],
         ylab = expression(paste("Correlation Coefficient (   ",italic(r)," )",sep='')),
         xlab = "Reporting Group", ylim = c(-1, 1), col = Colors15, names.arg = NA, main = "MHC")
legend("topright", legend = Groups15Short, fill = Colors15, bty = "n", cex = 0.8)

# GPDH
barplot2(height = KMA473PopsLinkageCorGPDH[5:dim(KMA473PopsLinkageCorGPDH)[1], "mean"], plot.ci = TRUE,
         ci.l = KMA473PopsLinkageCorGPDH[5:dim(KMA473PopsLinkageCorGPDH)[1], "2.5%"],
         ci.u = KMA473PopsLinkageCorGPDH[5:dim(KMA473PopsLinkageCorGPDH)[1], "97.5%"],
         ylab = expression(paste("Correlation Coefficient (   ",italic(r)," )",sep='')),
         xlab = "Reporting Group", ylim = c(-1, 1), col = Colors15, names.arg = NA, main = "GPDH")
legend("topright", legend = Groups15Short, fill = Colors15, bty = "n", cex = 0.8)

# Tf_ex
barplot2(height = KMA473PopsLinkageCorTf_ex[5:dim(KMA473PopsLinkageCorTf_ex)[1], "mean"], plot.ci = TRUE,
         ci.l = KMA473PopsLinkageCorTf_ex[5:dim(KMA473PopsLinkageCorTf_ex)[1], "2.5%"],
         ci.u = KMA473PopsLinkageCorTf_ex[5:dim(KMA473PopsLinkageCorTf_ex)[1], "97.5%"],
         ylab = expression(paste("Correlation Coefficient (   ",italic(r)," )",sep='')),
         xlab = "Reporting Group", ylim = c(-1, 1), col = Colors15, names.arg = NA, main = "Tf_ex")
legend("topright", legend = Groups15Short, fill = Colors15, bty = "n", cex = 0.8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Redo to get Population level correclation coefficients
source("V:/WORK/Kyle/R Scripts/LinkageCorrelationJAGSvKyle.GCL.R")  # Modified from Jim's function to provide population correlation coefficients

# MHC
KMA473PopsLinkageCor.rvalues.MHC <- LinkageCorrelationJAGSvKyle.GCL(sillyvec = KMA473Pops, markerset = c("One_MHC2_190", "One_MHC2_251"), groupvec = KMA473PopsGroupVec15, groupnames = Groups15Short); beep(8)
dput(x = KMA473PopsLinkageCor.rvalues.MHC, "Objects/KMA473PopsLinkageCor.rvalues.MHC.txt")
str(KMA473PopsLinkageCor.rvalues.MHC)
KMA473PopsLinkageCor.rvalues.MHC[1:2]
rownames(KMA473PopsLinkageCor.rvalues.MHC$Gst) <- c("Gst", "One_MHC2_190.One_MHC2_251.Gst", "One_MHC2_190.Gst", "One_MHC2_251.Gst")

# GPDH
KMA473PopsLinkageCor.rvalues.GPDH <- LinkageCorrelationJAGSvKyle.GCL(sillyvec = KMA473Pops, markerset = c("One_GPDH", "One_GPDH2"), groupvec = KMA473PopsGroupVec15, groupnames = Groups15Short); beep(8)
dput(x = KMA473PopsLinkageCor.rvalues.GPDH, "Objects/KMA473PopsLinkageCor.rvalues.GPDH.txt")
KMA473PopsLinkageCor.rvalues.GPDH[1:2]
rownames(KMA473PopsLinkageCor.rvalues.GPDH$Gst) <- c("Gst", "One_GPDH.One_GPDH2.Gst", "One_GPDH.Gst", "One_GPDH2.Gst")

# Tf_ex
KMA473PopsLinkageCor.rvalues.Tf_ex <- LinkageCorrelationJAGSvKyle.GCL(sillyvec = KMA473Pops, markerset = c("One_Tf_ex10-750", "One_Tf_ex3-182"), groupvec = KMA473PopsGroupVec15, groupnames = Groups15Short); beep(8)  # 
dput(x = KMA473PopsLinkageCor.rvalues.Tf_ex, "Objects/KMA473PopsLinkageCor.rvalues.Tf_ex.txt")
KMA473PopsLinkageCor.rvalues.Tf_ex[1:2]
rownames(KMA473PopsLinkageCor.rvalues.Tf_ex$Gst) <- c("Gst", "One_Tf_ex10-750.One_Tf_ex3-182.Gst", "One_Tf_ex10-750.Gst", "One_Tf_ex3-182.Gst")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## View population specific plots
# MHC
par(mar = c(5.1, 5.1, 4.1, 3.1))
barplot2(height = KMA473PopsLinkageCor.rvalues.MHC$Pop[, "mean"], width = 1, space = 0, plot.ci = TRUE,
         ci.l = KMA473PopsLinkageCor.rvalues.MHC$Pop[, "2.5%"],
         ci.u = KMA473PopsLinkageCor.rvalues.MHC$Pop[, "97.5%"],
         ylab = expression(paste("Correlation Coefficient (   ",italic(r)," )",sep='')),
         xlab = "Populations", ylim = c(-1, 1), col = Colors15[KMA473PopsGroupVec15],
         border = Colors15[KMA473PopsGroupVec15], names.arg = NA, main = "MHC")
legend(x = 474, y = 1, legend = Groups15Short, fill = Colors15, bty = "n", cex = 1, xpd = TRUE)

# GPDH
par(mar = c(5.1, 5.1, 4.1, 3.1))
barplot2(height = KMA473PopsLinkageCor.rvalues.GPDH$Pop[, "mean"], width = 1, space = 0, plot.ci = TRUE,
         ci.l = KMA473PopsLinkageCor.rvalues.GPDH$Pop[, "2.5%"],
         ci.u = KMA473PopsLinkageCor.rvalues.GPDH$Pop[, "97.5%"],
         ylab = expression(paste("Correlation Coefficient (   ",italic(r)," )",sep='')),
         xlab = "Populations", ylim = c(-1, 1), col = Colors15[KMA473PopsGroupVec15], 
         border = Colors15[KMA473PopsGroupVec15], names.arg = NA, main = "GPDH")
legend(x = 474, y = 1, legend = Groups15Short, fill = Colors15, bty = "n", cex = 1, xpd = TRUE)

# Tf_ex
par(mar = c(5.1, 5.1, 4.1, 3.1))
barplot2(height = KMA473PopsLinkageCor.rvalues.Tf_ex$Pop[, "mean"], width = 1, space = 0, plot.ci = TRUE,
         ci.l = KMA473PopsLinkageCor.rvalues.Tf_ex$Pop[, "2.5%"],
         ci.u = KMA473PopsLinkageCor.rvalues.Tf_ex$Pop[, "97.5%"],
         ylab = expression(paste("Correlation Coefficient (   ",italic(r)," )",sep='')),
         xlab = "Populations", ylim = c(-1, 1), col = Colors15[KMA473PopsGroupVec15], 
         border = Colors15[KMA473PopsGroupVec15], names.arg = NA, main = "Tf_ex")
legend(x = 474, y = 1, legend = Groups15Short, fill = Colors15, bty = "n", cex = 1, xpd = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Combine Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Contrasting results from fORCA and Gst ratios from LinakageCorellationsJAGC.GCL
## Also incorporating results from Kodiak 49 Pops (note this is with 11 Kodiak RGs)
Kodiak49PopsRandomfORCA <- dget("V:/WORK/Sockeye/Kodiak/Kodiak Afognak Baseline Look 2014/Kyle/Objects/Kodiak49PopsRandomfORCA.txt")


# MHC
KMA473PopsFullfORCA[intersect(grep("MHC2_251", names(deltafORCA)), grep("MHC2_190", names(deltafORCA))), ] # MHC2_190 > MHC2_251 > Combine
Kodiak49PopsRandomfORCA[intersect(grep("MHC2_251",rownames(Kodiak49PopsRandomfORCA)),grep("MHC2_190",rownames(Kodiak49PopsRandomfORCA))), ]  # Keep MHC2_190
KMA473PopsLinkageCor.rvalues.MHC$Gst  # Keep MHC2_190

# GPDH
KMA473PopsFullfORCA[grep("GPDH.One_GPDH", names(deltafORCA)), ]  # GPDH2 > GPDH > Combine
Kodiak49PopsRandomfORCA[grep("GPDH2.One_GPDH",rownames(Kodiak49PopsRandomfORCA)), ]  # Keep GPDH
KMA473PopsLinkageCor.rvalues.GPDH$Gst  # Keep GPDH

# Tf_ex
KMA473PopsFullfORCA[intersect(grep("Tf_ex10", names(deltafORCA)), grep("Tf_ex3", names(deltafORCA))), ]  # Combine > Tf_ex3 > Tf_ex10
Kodiak49PopsRandomfORCA[intersect(grep("Tf_ex10",rownames(Kodiak49PopsRandomfORCA)),grep("Tf_ex3",rownames(Kodiak49PopsRandomfORCA))), ]  # Keep Tf_ex3
KMA473PopsLinkageCor.rvalues.Tf_ex$Gst  # Keep Tf_ex3

## Jim suggests that we want to stick with fORCA as it is a better metric for MSA
## The Gst ratios were just supposed to save time, not muddy the water.
## It is worth visualizing the correlation coefficient plots by pop to verify
## And also looking to visualize haplotype frequencies for combined loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CombineLoci.GCL(sillyvec = KMA473Pops, markerset = c("One_MHC2_190", "One_MHC2_251"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = KMA473Pops, markerset = c("One_GPDH", "One_GPDH2"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = KMA473Pops, markerset = c("One_Tf_ex10-750", "One_Tf_ex3-182"), update = TRUE); beep(2)

# Get allele counts, sample sizes, frequencies
y <- FreqPop.GCL(sillyvec = KMA473Pops, loci = LocusControl$locusnames[which(LocusControl$nalleles == 9)])
str(y)

n <- SampSizeByLocus.GCL(sillyvec = KMA473Pops, loci = LocusControl$locusnames[which(LocusControl$nalleles == 9)])
n <- as.matrix(n)
str(n)

q <- array(data = NA, dim = dim(y), dimnames = dimnames(y))
for(allele in seq(dim(y)[3])) {
  q[, , allele] <- y[, , allele] / n
}
str(q)

table(apply(q, c(1,2), sum))
dimnames(q)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Haplotype frequency plot
require(lattice)
new.colors <- colorRampPalette(c("white", "black"))

# MHC
levelplot(q[, "One_MHC2_190.One_MHC2_251", ], 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.45, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.45, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 9.55, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 9.55, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          , col.regions = new.colors, xlab = "Pop", ylab = "Allele", aspect = "fill", main = "MHC")

# GDPH
levelplot(q[, "One_GPDH.One_GPDH2", ], 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.45, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.45, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 9.55, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 9.55, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          , col.regions = new.colors, xlab = "Pop", ylab = "Allele", aspect = "fill", main = "GPDH")

# Tf_ex
levelplot(q[, "One_Tf_ex10-750.One_Tf_ex3-182", ], 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.45, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.45, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 9.55, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 9.55, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          , col.regions = new.colors, xlab = "Pop", ylab = "Allele", aspect = "fill", main = "Tf_ex")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Linkage Decisions ####
# Made with Chris/Tyler/Bill on Thu Jul 30 13:12:27 2015
# 
## MHC
# DON'T COMBINE: We decided not to combine based on fORCA for all 473 Pops and 49 Kodiak Populations
# MHC2_190: We decided to keep MHC2_190, as it is more informative for the 49 Kodiak Populations, and for SEAK
#
## GPDH
# DON'T COMBINE: We decided not to combine based on fORCA for all 473 Pops and 49 Kodiak Populations
# GPDH: We decided to keep GPDH, as it is more informative for the 49 Kodiak Populations, Chignik, BB, and WASSIP
#
## Tf_ex
# COMBINE: We decided to combine based on fORCA for all 473 Pops, even though 49 Kodiak Populations suggests to NOT combine
# Tf_ex3: We decided to keep Tf_ex3, as it is more informative for the 49 Kodiak Populations, Chignik, BB, and WASSIP
# Try combined and alone and see!
# 
# Confirmed with Chris on Fri Aug 14 11:07:52 2015 that these are the only 3 linked pairs of loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Remove linked markers
#### CONFIRM THAT THESE ARE INDEED LINKED IN THE WHOLE BASELINE ####

KMA473Pops_LD_locusfail_remove <- c("One_MHC2_251", "One_GPDH2", "One_Tf_ex10-750", "One_Tf_ex3-182", names(mito.loci))
dput(x = KMA473Pops_LD_locusfail_remove, file = "Objects/KMA473Pops_LD_locusfail_remove.txt")

KMA473Pops_LD_locuscombine_add <- c("One_Tf_ex10-750.One_Tf_ex3-182", "One_CO1.One_Cytb_17.One_Cytb_26")
dput(x = KMA473Pops_LD_locuscombine_add, file = "Objects/KMA473Pops_LD_locuscombine_add.txt")


loci89 <- c(loci94[!loci94 %in% KMA473Pops_LD_locusfail_remove], KMA473Pops_LD_locuscombine_add)
dput(x = loci89, file = "Objects/loci89.txt")


## Need to combine the three mitochondrial markers for Jim's GstRatiosForThreeHaploidSNPs.GCL function
# This will determine which mitochondrial markers to keep/combine
CombineLoci.GCL(sillyvec = KMA473Pops, markerset = c("One_CO1", "One_Cytb_17", "One_Cytb_26"), update = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Before determining final markerset, we want to evaluate if it is worth keeping all 3 mitochondrial markers ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Need to combine the three mitochondrial markers for Jim's GstRatiosForThreeHaploidSNPs.GCL function
# This will determine which mitochondrial markers to keep/combine
CombineLoci.GCL(sillyvec = KMA473Pops, markerset = c("One_CO1", "One_Cytb_17", "One_Cytb_26"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = KMA473Pops, markerset = c("One_CO1", "One_Cytb_17"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = KMA473Pops, markerset = c("One_CO1", "One_Cytb_26"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = KMA473Pops, markerset = c("One_Cytb_17", "One_Cytb_26"), update = TRUE); beep(2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save post-QC Pop .gcl's as back-up:
LocusControl$locusnames
dput(x = LocusControl, file = "Objects/LocusControl103.txt")

loci103 <- LocusControl$locusnames
dput(x = loci103, file = "Objects/loci103.txt")

dir.create(path = "Raw genotypes/PostQCPopsloci103")
invisible(sapply(KMA473Pops, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostQCPopsloci103/" , silly, ".txt", sep = ''))} )); beep(8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## All three mitocondrial markers combined
# Get allele counts, sample sizes, frequencies
y <- FreqPop.GCL(sillyvec = KMA473Pops, loci = LocusControl$locusnames[which(LocusControl$nalleles == 8)])
str(y)

n <- SampSizeByLocus.GCL(sillyvec = KMA473Pops, loci = LocusControl$locusnames[which(LocusControl$nalleles == 8)])
n <- as.matrix(n)
str(n)

q <- array(data = NA, dim = dim(y), dimnames = dimnames(y))
for(allele in seq(dim(y)[3])) {
  q[, , allele] <- y[, , allele] / n
}
str(q)

table(apply(q, c(1,2), sum))
dimnames(q)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Haplotype frequency plot
require(lattice)
new.colors <- colorRampPalette(c("white", "black"))

# One_CO1.One_Cytb_17.One_Cytb_26
levelplot(q[, "One_CO1.One_Cytb_17.One_Cytb_26", ], 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.45, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.45, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 8.55, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 8.55, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          , col.regions = new.colors, xlab = "Pop", ylab = "Allele", aspect = "fill", main = "One_CO1.One_Cytb_17.One_Cytb_26")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Combinations of 2 mitochondrial markers

y <- FreqPop.GCL(sillyvec = KMA473Pops, loci = LocusControl$locusnames[which(LocusControl$nalleles == 4)])
str(y)

n <- SampSizeByLocus.GCL(sillyvec = KMA473Pops, loci = LocusControl$locusnames[which(LocusControl$nalleles == 4)])
n <- as.matrix(n)
str(n)

q <- array(data = NA, dim = dim(y), dimnames = dimnames(y))
for(allele in seq(dim(y)[3])) {
  q[, , allele] <- y[, , allele] / n
}
str(q)

table(apply(q, c(1,2), sum))
dimnames(q)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Haplotype frequency plot
# One_CO1.One_Cytb_17
levelplot(q[, "One_CO1.One_Cytb_17", ], 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.45, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.45, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 4.55, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 4.55, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          , col.regions = new.colors, xlab = "Pop", ylab = "Allele", aspect = "fill", main = "One_CO1.One_Cytb_17")

# One_CO1.One_Cytb_26
levelplot(q[, "One_CO1.One_Cytb_26", ], 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.45, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.45, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 4.55, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 4.55, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          , col.regions = new.colors, xlab = "Pop", ylab = "Allele", aspect = "fill", main = "One_CO1.One_Cytb_26")

# One_Cytb_17.One_Cytb_26
levelplot(q[, "One_Cytb_17.One_Cytb_26", ], 
          panel = function(...){
            panel.levelplot(...)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 0.45, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 0.45, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
            panel.arrows(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + c(0.75, rep(0.5, 14)), y0 = 4.55, x1 = cumsum(table(KMA473PopsGroupVec15)) + c(rep(0.5, 14), 1.25), y1 = 4.55, col = Colors15, lwd = 4, length = 10, code = 3, angle = 90)
          },
          , col.regions = new.colors, xlab = "Pop", ylab = "Allele", aspect = "fill", main = "One_Cytb_17.One_Cytb_26")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Jim's new GstRatio function
GstRatios_3mitochondrialSNPs <- GstRatiosForThreeHaploidSNPs.GCL(sillyvec = KMA473Pops, combinedhaploidmarker = "One_CO1.One_Cytb_17.One_Cytb_26", groupvec = KMA473PopsGroupVec15)
dput(x = GstRatios_3mitochondrialSNPs, file = "Objects/GstRatios_3mitochondrialSNPs.txt")
GstRatios_3mitochondrialSNPs  # One_CO1 is clearly the best choice here

## The Gst ratio function uses "haploid gene diversity" as heterozygosity (see Hedrick pg 99-100)
## And takes the ratio of Gst for RGs for different markers
## While this initially seems intuitive, the maximum value of Gst diminishes with increasing numbers of alleles,
## such that Gst no longer ranges from 0 -> 1, but 0 -> GstMAX < 1

## Jim is working on an extension of fORCA to handle haploid data so that we can determine whether it is "worth it"
## To keep/combine multiple haploid markers.

# KMA473Pops_mitochondrial_locusfail_remove <- names(mito.loci)[-1]
# dput(x = KMA473Pops_mitochondrial_locusfail_remove, file = "Objects/KMA473Pops_mitochondrial_locusfail_remove.txt")
# 
# loci89 <- loci91[!loci91 %in% KMA473Pops_mitochondrial_locusfail_remove]
# dput(x = loci89, file = "Objects/loci89.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Jim's update to fORCA.GCL to handle haploid markers
# No distribution to compare it to, but still get a handle on what is helpful
allmito.loci <- LocusControl$locusnames[LocusControl$ploidy == 1][-c(4:6)]

allmito.473Pops.15RG.fORCA <- fORCA.GCL(sillyvec = KMA473Pops, loci = allmito.loci, groupvec = KMA473PopsGroupVec15, NREP = 100000); beep(8)
dput(x = allmito.473Pops.15RG.fORCA, file = "Objects/allmito.473Pops.15RG.fORCA.txt")


Kodiak57PopsGroupVec9 <- KMA473PopsGroupVec15[which(KMA473PopsGroupVec15 >= 4 & KMA473PopsGroupVec15 <= 12)] - 3
dput(x = Kodiak57PopsGroupVec9, file = "Objects/Kodiak57PopsGroupVec9.txt")
Kodiak57Pops <- KMA473Pops[which(KMA473PopsGroupVec15 >= 4 & KMA473PopsGroupVec15 <= 12)]
dput(x = Kodiak57Pops, file = "Objects/Kodiak57Pops.txt")
Groups9.nospace.Kodiak <- Groups15.nospace[4:12]
dput(x = Groups9.nospace.Kodiak, file = "Objects/Groups9.nospace.Kodiak.txt")


allmito.50Pops.9RG.fORCA <- fORCA.GCL(sillyvec = Kodiak57Pops, loci = allmito.loci, groupvec = Kodiak57PopsGroupVec9, NREP = 100000); beep(8)
dput(x = allmito.50Pops.9RG.fORCA, file = "Objects/allmito.50Pops.9RG.fORCA.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Clean workspace and dget Population .gcls, LocusControl, and necessary objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")
# source("V:\\Analysis\\R files\\Scripts\\PROD\\Functions.GCL.r")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 96 Loci
## Locus Control
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
loci96 <- dget(file = "Objects/loci96.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get Populations
KMA473Pops <- dget(file = "Objects/KMA473Pops.txt")

require(beepr)
invisible(sapply(KMA473Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPops/", silly, ".txt", sep = "")), pos = 1)})); beep(2)
length(objects(pattern = "\\.gcl"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
# ## Dput all necessary items for loci103
# dput(x = LocusControl, file = "Objects/LocusControl103.txt")
# 
# ## Save post-QC Pop 103loci .gcl's as back-up:
# dir.create(path = "Raw genotypes/PostQCPopsloci103")
# sapply(KMA473Pops, function(silly) {dput(x = get(paste(silly, ".gcl", sep = '')), file = paste("Raw genotypes/PostQCPopsloci103/" , silly, ".txt", sep = ''))}); beep(8)

# RG size for hypothetical 31 reporting groups
cbind(setNames(object = sapply(1:31, function(i) {sum(sapply(KMA473Pops[KMA473PopsGroupVec31 == i], function(silly) {get(paste0(silly, ".gcl"))$n}))}), nm = KMA31GroupsPC))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summary Statistics ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # All 96 loci
# SumStats_KMA473Pops_96loci <- HoFisFstTable.GCL(sillyvec = KMA473Pops, loci = loci96, dir = "Tables")
# dput(x = SumStats_KMA473Pops_96loci, file = "Objects/SumStats_KMA473Pops_96loci.txt"); beep(8)
# str(SumStats_KMA473Pops_96loci)
# 
# # Write table for baseline .xlsx file
# require(xlsx)
# write.xlsx(x = SumStats_KMA473Pops_96loci, file = "Tables/HoFisFstTable_KMA473Pops_loci96.xlsx")
# 
# # Order by Fst
# SumStats_KMA473Pops_96loci[order(SumStats_KMA473Pops_96loci$Fst[seq(loci96)]), ]
# 
# # Histogram of Fst
# hist(SumStats_KMA473Pops_96loci[seq(loci96), "Fst"], col = 8, xlab = "Fst", main = "Histogram of Fst")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# All possible 103 loci (linked diploid and haploid are combined)
SumStats_KMA473Pops_103loci <- HoFisFstTable.GCL(sillyvec = KMA473Pops, loci = loci103, dir = "Tables"); beep(8)
dput(x = SumStats_KMA473Pops_103loci, file = "Objects/SumStats_KMA473Pops_103loci.txt"); beep(2)
write.xlsx(x = SumStats_KMA473Pops_103loci, file = "Output/HoFisFstTable_KMA473Pops_loci103.xlsx")

# Order by Fst
SumStats_KMA473Pops_103loci[order(SumStats_KMA473Pops_103loci$Fst[seq(loci96)]), ]

# Histogram of Fst
hist(SumStats_KMA473Pops_103loci[order(SumStats_KMA473Pops_103loci$Fst[seq(loci96)]), "Fst"], col = 8, xlab = "Fst", main = "Histogram of Fst")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Final set of 89 loci (linked diploid and haploid are combined) to get overall values
SumStats_KMA473Pops_89loci <- HoFisFstTable.GCL(sillyvec = KMA473Pops, loci = loci89, dir = "Tables"); beep(8)
dput(x = SumStats_KMA473Pops_89loci, file = "Objects/SumStats_KMA473Pops_89loci.txt"); beep(2)
write.xlsx(x = SumStats_KMA473Pops_89loci, file = "Output/HoFisFstTable_KMA473Pops_loci89.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Final set of 46 loci (linked diploid and haploid are combined) to get overall values
SumStats_KMA473Pops_46loci <- HoFisFstTable.GCL(sillyvec = KMA473Pops, loci = loci46, dir = "Tables"); beep(8)
dput(x = SumStats_KMA473Pops_46loci, file = "Objects/SumStats_KMA473Pops_46loci.txt"); beep(2)
write.xlsx(x = SumStats_KMA473Pops_46loci, file = "Output/HoFisFstTable_KMA473Pops_loci46.xlsx")


#~~~~~~~~~~~~~~~~~~
# Get variance components to speed things up
require(package = hierfstat)

# KMA473Pops103loci4levelVarCompByLocus <- dget(file = "Variance Components/KMA473Pops103loci4levelVarCompByLocus.txt")
# str(KMA473Pops103loci4levelVarCompByLocus)
# 
# KMA473Pops15Group103lociVarComps <- dget(file = "Variance Components/KMA473Pops15Group103lociVarComps.txt")
# str(KMA473Pops15Group103lociVarComps)

KMA473Pops0Group103lociVarComps <- dget(file = "Variance Components/KMA473Pops0Group103lociVarComps.txt")
str(KMA473Pops0Group103lociVarComps)

loci <- loci46
nloci <- length(loci)
ploidy <- LocusControl$ploidy[loci]

MyVC <- KMA473Pops0Group103lociVarComps
dimnames(MyVC)[[2]] <- c("P", "I", "G")

MyTable=MyVC[loci, c("P", "I")]/apply(MyVC[loci, ], 1, sum)

MyTable[loci[ploidy==1], "I"] <- 0

MyTable=rbind(MyTable,Overall=apply(MyVC[loci,c("P","I")],2,sum)/sum(apply(MyVC[loci,],1,sum)))

Ho=rep("-",nloci+1)

names(Ho)=c(loci,"Overall")

Ho[loci[ploidy==2]] <- as.numeric(as.character(SumStats_KMA473Pops_89loci[loci[ploidy==2], "Ho"]))

# Ho[loci[ploidy==2]]=apply(basic.stats(dat[,c("Pop",loci[ploidy==2])])$Ho,1,mean)

Ho["Overall"]=mean(as.numeric(Ho[loci[ploidy==2]]))

MyTable=data.frame(Ho=Ho,Fis=MyTable[,"I"],Fst=MyTable[,"P"])

SumStats_KMA473Pops_46loci <- MyTable
dput(x = SumStats_KMA473Pops_46loci, file = "Objects/SumStats_KMA473Pops_46loci.txt")
write.xlsx(x = SumStats_KMA473Pops_46loci, file = "Output/HoFisFstTable_KMA473Pops_loci46.xlsx")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Pairwise Fst Tree ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run with "final" markerset out of 96 (remove HWE, LD, mitochondrial)
# KMA473Pops89lociFstTree <- PairwiseFstTree.GCL(sillyvec = KMA473Pops, loci = loci89, dir = "Trees" , nboots = 1000)
# dput(x = KMA473Pops89lociFstTree, file = "Objects/KMA473Pops89lociFstTree.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Takes too long, need to do in pieces of pairs
dir.create(path = "Trees/FstTreeObjects")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## dget Tree output
# Do not use as the "bootstrapFst" objects are MASSIVE
# invisible(sapply(list.files(path = "Trees/FstTreeObjects"), function(fle) {assign(x = unlist(strsplit(x = fle, split = ".txt")), value =  dget(file = paste(getwd(), "/Trees/FstTreeObjects/", fle, sep = '')), pos = 1)} )); beep(8)
# str(Fst_1); str(vc_1, max.level = 0)
# seq(objects(pattern = "Fst_"))

rm(list = ls(all = TRUE))
setwd("V:/WORK/Sockeye/Kodiak/KMA Commercial Harvest 2014-2016/Baseline")
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")
source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")


# Combine Fst objects
Fst_files <- list.files(path = "Trees/FstTreeObjects", pattern = "Fst_")
Fst_files <- Fst_files[-grep(pattern = "bootstrap", x = Fst_files)]
Fst <- lapply(seq(Fst_files), function(i) {dget(file = paste(getwd(), "/Trees/FstTreeObjects/Fst_", i, ".txt", sep = ""))} ); #beep(2)
Fst <- Reduce(f = "+", x = Fst)
str(Fst)
# Fst <- Reduce(f = "+", x = sapply(1:22, function(i) {get(paste("Fst", i, sep = "_"))}, simplify = FALSE))
# write.table(x = Fst, file = "Output/KMA473Pops89lociPairwiseFst.txt")

# require(lattice)
# new.colors <- colorRampPalette(c("white", "black"))
# levelplot(Fst, col.regions = new.colors, xlab = "SILLY", ylab = "SILLY", main = "Fst 473Pops 89loci", scales = list(draw = FALSE), aspect = "fill", Rowv = FALSE)

# Combine vc objects into one list
vc <- lapply(seq(list.files(path = "Trees/FstTreeObjects", pattern = "vc_")), function(i) {dget(file = paste(getwd(), "/Trees/FstTreeObjects/vc_", i, ".txt", sep = ""))} ); #beep(8)
vc <- lapply(vc, function(vc_i) {vc_i[which(sapply(vc_i, function(pair) {!is.null(pair)} ) == TRUE)]} )
vc <- do.call(what = "c", args = vc)
str(vc, max.level = 0); choose(n = length(KMA473Pops), k = 2)

require(ape); require(hierfstat)
tree <- nj(Fst)

# # Define objects
# nboots = 1000
# sillyvec <- KMA473Pops
# loci <- loci89
# 
# nloci <- length(loci)
# nsillys <- length(sillyvec)
# pairs <- combn(sillyvec, 2)
# pairnames <- apply(pairs, 2, function(col) {paste(col, collapse = ".")} )
# dimnames(pairs)[[2]] <- pairnames
# 
# 
# # Finish tree
# require(ape); require(hierfstat)
# 
# tree <- nj(Fst)
# 
# trees <- bootstrapFst <- vector("list", nboots)
# 
# cat("2 Main function tasks: \n1) Bootstrap\n2) each pair of SILLYs\n")
# if (.Platform$OS.type == "windows") flush.console()
# 
# cat("\nBootstrap over nboots (1000)\n", sep = '')
# if (.Platform$OS.type == "windows") flush.console()
# pb <- txtProgressBar(min = 0, max = nboots, style = 3)
# 
# 
# # Did this in piecies (100 boots at time) on the server
# for(boot in 1:nboots) {
#   setTxtProgressBar(pb = pb, value = boot)
#   
#   # Sample loci with replacement
#   temploci <- sample(loci, nloci, replace = TRUE)
#   
#   # Determine pairwise Fst
#   tempFstVec <- sapply(pairnames, function(pair) {sum(vc[[pair]][temploci, 1]) / sum(vc[[pair]][temploci, 1:3])} )
#   
#   # Translate vector into Fst matrix
#   tempFst <- array(0, c(nsillys, nsillys), dimnames = list(sillyvec, sillyvec))
#   tempFst[lower.tri(tempFst, diag = FALSE)] <- tempFstVec
#   tempFst <- tempFst + t(tempFst)
#   
#   # Save this bootstrap of the Fst matrix
#   bootstrapFst[[boot]] <- tempFst
#   
#   # Create a tree for this bootstrap
#   trees[[boot]] <- nj(tempFst)   
# }

## dget "bootstrapFst" and "trees" objects
# bootstrap_files <- list.files(path = "Trees/FstTreeObjects", pattern = "bootstrap")
# bootstrapFst <- do.call(what = "c", args = lapply(seq(bootstrap_files), function(i) {dget(file = paste(getwd(), "/Trees/FstTreeObjects/bootstrapFst_", i, ".txt", sep = ""))} ))
# str(bootstrapFst, max.level = 1)

trees_files <- list.files(path = "Trees/FstTreeObjects", pattern = "trees")
trees <- lapply(seq(trees_files), function(i) {dget(file = paste(getwd(), "/Trees/FstTreeObjects/trees_", i, ".txt", sep = ""))} ); #beep(2)
trees <- do.call(what = "c", args = trees)
str(trees, max.level = 0)

date()
ptm <- proc.time()
bootstrap <- prop.clades(tree,trees); #beep(8)
proc.time() - ptm

# Save tree
KMA473Pops89lociFstTree <- list(tree = tree, bootstrap = bootstrap, PairwiseFst = Fst, vc = vc)  # , BootstrapFst = bootstrapFst # Too big
dput(x = KMA473Pops89lociFstTree, file = "Trees/KMA473Pops89lociFstTree.txt")
dput(x = KMA473Pops89lociFstTree[1:3], file = "Trees/KMA473Pops89lociFstTreeNoVC.txt")
str(KMA473Pops89lociFstTree, max.level = 1)

# Cleanup
rm(list = c(objects(pattern = "Fst_"), 
            objects(pattern = "vc_"),
            objects(pattern = "trees_"),
            vc, Fst, nboots, sillyvec, loci, nloci, nsillys, pairs, pairnames, tree, trees, boot, tempFst, temploci, pair, sillys, bootstrapFst)
   )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA473Pops89lociFstTree <- dget(file = "Trees/KMA473Pops89lociFstTree.txt")

str(KMA473Pops89lociFstTree)
CommonNames473

# Check out bootstrap support
KMA473Pops89lociFstTree$bootstrap

# Check out edge.length
KMA473Pops89lociFstTree$tree$edge.length <- pmax(0,KMA473Pops89lociFstTree$tree$edge.length) # Get rid of negative branches
KMA473Pops89lociFstTree$tree$bootstrap <- KMA473Pops89lociFstTree$bootstrap

# Add common names to tip labels
KMA473Pops89lociFstTree$tree$tip.label <- CommonNames473


# Write a tree
require(ape)
write.tree(KMA473Pops89lociFstTree$tree, "Trees/KMA473Pops89lociFstTree.nex")

## Plot in R
# Need to run colortree line by line to get colors to work
colortree <- treeColor.GCL(tree = KMA473Pops89lociFstTree$tree, currentnames = CommonNames473, treenames = CommonNames473, groupvec = KMA473PopsGroupVec15, regioncol = sapply(Colors15, function(color) {which(colors() == color)}), regionpch = NULL)
plot.phylo(colortree$tree, edge.color = colortree$color, edge.width = 6, use.edge.length = TRUE, show.tip.label = TRUE, adj = 0.015, cex = 1, font = 1, label.offset = 0.001)

# Add bootstrap values at the nodes in order to place on FigTree plot

par(xpd = TRUE)
boots <- ifelse(KMA473Pops89lociFstTree$bootstrap > 500, KMA473Pops89lociFstTree$bootstrap, "")
nodelabels(boots, frame = "none", adj = c(0, 0), col = "black", cex = 1)  # shows consensus node labels on tree
par(xpd = FALSE)
axis(1)  #Adds scale to bottom of plot
mtext(text = expression(italic(F)[ST]), side = 1, cex = 1.5, outer = FALSE, line = 3)

# Get RGB colors
names(Colors15) <- Groups15Short
Colors15RGB <- t(col2rgb(col = Colors15))
dput(x = Colors15RGB, file = "Objects/Colors15RGB.txt")

# Get RG sample sizes
RGSampleSizes <- setNames(object = sapply(seq(Groups15), function(i) {sum(sapply(KMA473Pops[which(KMA473PopsGroupVec15 == i)], function(silly) {get(paste(silly, ".gcl", sep = ""))$n} ))} ), nm = Groups15Short)
getwd()
write.table(RGSampleSizes, file = "sampsize.txt", sep = ",")
# Saltery and Uganik need more fish!






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Tree for manuscript

KMA473Pops89lociFstTree <- dget(file = "Trees/KMA473Pops89lociFstTree.txt")
str(KMA473Pops89lociFstTree, max.level = 1)

# Create a colored Fst tree for final figures
require(ape)
par(mar = c(5, 5, 5, 5), oma = c(3, 0, 0, 0))
plot.phylo(x = KMA473Pops89lociFstTree$tree, cex = 0.5, no.margin = TRUE, type = "p")
axisPhylo(1, las = 1, backward = FALSE)

str(KMA473Pops89lociFstTree$tree, max.level = 1)

Colors14 <- c(Colors15[1:4], "dodgerblue", Colors15[7:15])
dput(x = Colors14, file = "Objects/Colors14.txt")

Colors14RGB <- Colors15RGB[c(1:4, 6:15), ]
Colors14RGB[5, ] <- c(30, 144, 255)
rownames(Colors14RGB)[5] <- "Ayak/Fzr"
dput(x = Colors14RGB, file = "Objects/Colors14RGB.txt")
Colors14RGB[4:11, ]



ColorTree <- treeColor.GCL(tree = KMA473Pops89lociFstTree$tree, currentnames = KMA473Pops, treenames = CommonNames473, groupvec = KMA473PopsGroupVec14, regioncol = match(Colors14,colors()), regionpch = NULL)
dput(x = ColorTree, file = "Trees/ColorTreeFstKMA473Pops14RG89loci.txt")
str(ColorTree)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Big Tree, no legend
par(family="serif")
require(ape)
require(devEMF)

# 14 RGs
emf(file = "Trees/KMA473Pops14RG89LociFstColorTree.emf", width = 4, height = 6)
par(oma = c(4, 0, 0, 0))
plot.phylo(x = ColorTree$tree, edge.color = ColorTree$color, edge.width = 3, use.edge.length = TRUE, 
           show.tip.label = FALSE, adj = 0.05, font = 0.2, cex = 0.6, no.margin = TRUE, type = "p")  # cex = 0.8
axisPhylo(1, las = 1, backward = FALSE)
mtext(text=expression(italic(F)[ST]), side = 1, cex = 1.5, outer = FALSE, line = 3, adj = 0.4)
dev.off()



# 14 RGs
emf(file = "Trees/KMA473Pops14RG89LociFstColorTreeBig.emf", width = 8, height = 50)
par(oma = c(4, 0, 0, 0))
plot.phylo(x = ColorTree$tree, edge.color = ColorTree$color, edge.width = 3, use.edge.length = TRUE, 
           show.tip.label = TRUE, adj = 0.05, font = 0.2, cex = 0.6, no.margin = TRUE, type = "p")  # cex = 0.8
axisPhylo(1, las = 1, backward = FALSE)
mtext(text=expression(italic(F)[ST]), side = 1, cex = 1.5, outer = FALSE, line = 3, adj = 0.4)
dev.off()


# 6 RGs
Colors6 <- c(Colors15[1], rep("brown4", 2), rep("green4", 8), Colors15[13:15])

ColorTree <- treeColor.GCL(tree = KMA473Pops89lociFstTree$tree, currentnames = KMA473Pops, treenames = CommonNames473, groupvec = KMA473PopsGroupVec14, regioncol = match(Colors6,colors()), regionpch = NULL)
dput(x = ColorTree, file = "Trees/ColorTreeFstKMA473Pops6RG89loci.txt")
emf(file = "Trees/KMA473Pops6RG89LociFstColorTree.emf", width = 4, height = 6)
par(oma = c(4, 0, 0, 0))
plot.phylo(x = ColorTree$tree, edge.color = ColorTree$color, edge.width = 3, use.edge.length = TRUE, 
           show.tip.label = FALSE, adj = 0.05, font = 0.2, cex = 0.6, no.margin = TRUE, type = "p")  # cex = 0.8
axisPhylo(1, las = 1, backward = FALSE)
mtext(text=expression(italic(F)[ST]), side = 1, cex = 1.5, outer = FALSE, line = 3, adj = 0.4)
dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Kodiak57Pops89lociFstTree <- PairwiseFstTree.GCL(sillyvec = Kodiak57Pops, loci = loci89, dir = "Trees" , nboots = 1000)

Kodiak57Pops89lociFstTree <- dget("Trees/Kodiak57Pops89lociFstTree.txt")

Kodiak57Pops89lociFstTree$tree$edge.length <- pmax(0,Kodiak57Pops89lociFstTree$tree$edge.length) # Get rid of negative branches

Kodiak57PopsGroupVec8 <- c(Kodiak57PopsGroupVec9[1:16], Kodiak57PopsGroupVec9[17:57] - 1)
ColorTree <- treeColor.GCL(tree = Kodiak57Pops89lociFstTree$tree, currentnames = Kodiak57Pops, treenames = CommonNames473[167:(167+57)], groupvec = Kodiak57PopsGroupVec8, regioncol = match(Colors14[4:11],colors()), regionpch = NULL)
dput(x = ColorTree, file = "Trees/ColorTreeFstKodiak57Pops8RGRG89loci.txt")

emf(file = "Trees/Kodiak57Pops8RG89LociFstColorTree.emf", width = 4, height = 6)
par(oma = c(4, 0, 0, 0))
plot.phylo(x = ColorTree$tree, edge.color = ColorTree$color, edge.width = 3, use.edge.length = TRUE, 
           show.tip.label = FALSE, adj = 0.05, font = 0.2, cex = 0.6, no.margin = TRUE, type = "p", root.edge = TRUE)  # cex = 0.8
axisPhylo(1, las = 1, backward = FALSE)
mtext(text=expression(italic(F)[ST]), side = 1, cex = 1.5, outer = FALSE, line = 3, adj = 0.4)
dev.off()

ColorTree$tree$tip.label[36] <- "Moraine Creek"

emf(file = "Trees/Kodiak57Pops8RG89LociFstColorTreeLabels.emf", width = 4, height = 6)
par(oma = c(4, 0, 0, 0))
plot.phylo(x = ColorTree$tree, edge.color = ColorTree$color, edge.width = 3, use.edge.length = TRUE, 
           show.tip.label = TRUE, adj = 0.05, font = 0.2, cex = 0.5, no.margin = TRUE, type = "p", root.edge = TRUE)  # cex = 0.8
axisPhylo(1, las = 1, backward = FALSE)
mtext(text=expression(italic(F)[ST]), side = 1, cex = 1.5, outer = FALSE, line = 3, adj = 0.4)
legend("topright", legend = KMA14GroupsPC2[4:11], fill = Colors14[4:11], bty = "n", cex = 0.7)  
dev.off()

emf(file = "Trees/Legend14RG.emf", width = 4, height = 6)
plot.new()
legend("topright", legend = KMA14GroupsPC2, fill = Colors14, bty = "n", cex = 1)  
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Likelihood Profiles 89 Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("V:/WORK/Kyle/R Scripts/LeaveOneOutDistKS.GCL.R")
dir.create("Likelihood Profiles")

# KMA473Pops_89loci_Likelihood_Profile <- LeaveOneOutDistKS.GCL(sillyvec = KMA473Pops, loci = loci89, groupvec = KMA473PopsGroupVec15)
# dput(x = KMA473Pops_89loci_Likelihood_Profile, file = "Likelihood Profiles/KMA473Pops_89loci_Likelihood_Profile.txt")

KMA473Pops_89loci_Likelihood_Profile <- dget(file = "Likelihood Profiles/KMA473Pops_89loci_Likelihood_Profile.txt")
str(KMA473Pops_89loci_Likelihood_Profile[[1]], max.level = 1)

# Individual baseline genetic likelihood for RGs assigned (i.e. what is the genotype likelihood of all baseline individuals to RG X)
invisible(lapply(KMA473Pops_89loci_Likelihood_Profile[[1]], function(lst) {boxplot(lst, pch=16, cex = 0.5, notch=TRUE, col=Colors15[sort(KMA473PopsGroupVec15)], ylab="Probability", xlim = c(0, length(KMA473Pops) * 1.17), bty = "n", axes = FALSE, xlab = "Population", cex.lab = 1.5)
                                                                           axis(side = 2, cex.axis = 1.5)
                                                                           axis(side = 1, cex.axis = 1.5, at = c(seq(from = 0, to = length(KMA473Pops), by = 50), length(KMA473Pops)))
                                                                           
                                                                           segments(x0 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + 0.5, y0 = 0, x1 = c(0, cumsum(table(KMA473PopsGroupVec15))[-15]) + 0.5, y1 = 1, col = Colors15, lwd = 4)
                                                                           
                                                                           legend(x = length(KMA473Pops), y = 1, legend = Groups15, fill = Colors15, bty = "n", cex = 1.3)} ))  # Regional Flat Prior
# Confusion Matrices
KMA473Pops_89loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = KMA473Pops_89loci_Likelihood_Profile, groupnames = Groups15, groupvec = KMA473PopsGroupVec15, sillyvec = KMA473Pops)  # Regional Flat Prior
dput(x = KMA473Pops_89loci_Confusion, file = "Objects/KMA473Pops_89loci_Confusion.txt")
KMA473Pops_89loci_Confusion <- dget(file = "Objects/KMA473Pops_89loci_Confusion.txt")

KMA473Pops_89loci_Confusion[[1]]
dimnames(KMA473Pops_89loci_Confusion[[1]]) <- list(PCGroups15, PCGroups15)

require(lattice)
require(devEMF)

new.colors <- colorRampPalette(c("white", "black"))

emf(file = "Likelihood Profiles/KMA473Pops_89loci_Confusion.emf", width = 6.5, height = 6.5, family = "Times")
levelplot(KMA473Pops_89loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", 
          at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)))
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# New LeaveOneOutLikeProfile.GCL.R
source("V:/DATA/R_GEN/JJs GCL/LeaveOneOutlIkeProfile.GCL.R")
# KMA473Pops_89loci_Likelihood_Profile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = KMA473Pops, loci = loci89, groupvec = KMA473PopsGroupVec15, groupnames = Groups15, groupcomps = NULL, ncores = 6)
# dput(x = KMA473Pops_89loci_Likelihood_Profile_NEW, file = "Likelihood Profiles/KMA473Pops_89loci_Likelihood_Profile_NEW.txt")
KMA473Pops_89loci_Likelihood_Profile_NEW <- dget(file = "Likelihood Profiles/KMA473Pops_89loci_Likelihood_Profile_NEW.txt")
str(KMA473Pops_89loci_Likelihood_Profile_NEW)

source("V:/DATA/R_GEN/JJs GCL/PlotLikeProfile.GCL.R")
PlotLikeProfile.GCL(likeprof = KMA473Pops_89loci_Likelihood_Profile_NEW, popvec = KMA473Pops, loci = loci89, groupvec = KMA473PopsGroupVec15, groupnames = Groups15, dir = "Likelihood Profiles")
likeprof = KMA473Pops_89loci_Likelihood_Profile_NEW; popvec = KMA473Pops; loci = loci89; groupvec = KMA473PopsGroupVec15; groupnames = Groups15; dir = "Likelihood Profiles"
col = Colors15


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MDS 89 Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Want to confirm that the MDS look good before I spend a lot of time on repeated proof tests
## Dump a GENEPOP file for genind
gcl2Genepop.GCL(sillyvec = KMA473Pops, loci = loci89, path = "Genepop/KMA473Pops_89loci.gen", VialNums = TRUE)

require(adegenet)
genind <- read.genepop(file = "Genepop/KMA473Pops_89loci.gen")

genpop <- genind2genpop(genind)

AdegenetNei473Pop89loci <- dist.genpop(genpop, method = 1, diag = TRUE, upper = TRUE)
dput(x = AdegenetNei473Pop89loci, file = "Trees/AdegenetNei473Pop89loci.txt")
str(AdegenetNei473Pop89loci)

require(ape)
Nei473NJtree <- nj(AdegenetNei473Pop89loci)
str(Nei473NJtree)
plot.phylo(x = Nei473NJtree, cex = 0.5, no.margin = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MDS
library('rgl')

# MDS <- cmdscale(as.matrix(AdegenetNei473Pop89loci), k = 3)  # Did in base R as it kept crashing in RStudio...dunno why.
# dput(x = MDS, file = "Objects/MDSAdegenetNei473Pop89loci.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei473Pop89loci.txt")

str(MDS)
x <- as.vector(MDS[, 1])   
y <- as.vector(MDS[, 2])
z <- as.vector(MDS[, 3])

KMA473PopsGroupVec15 <- dget(file = "Objects/KMA473PopsGroupVec15.txt")
Colors15 <- dget(file = "Objects/Colors15.txt")


open3d()
par(family = "serif")
plot3d(x, y, z + abs(range(z)[1]), xlab = '', ylab = '', zlab = '', aspect = FALSE, col = Colors15[KMA473PopsGroupVec15], size = 0.5, type = 's', axes = TRUE, box = TRUE, top = TRUE, cex = 1)
box3d()
# plot3d(x, y, z + abs(range(z)[1]), aspect = FALSE, col = "black", size = 0.5, type = 'h', box = TRUE, axes = FALSE, top = FALSE, add = TRUE)  # adds pins to spheres
texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.2, 0), text = seq(KMA473Pops), font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)

# Plot individually
par(mar = c(4.1, 4.1, 0.6, 0.6))
plot(x, col = Colors15[KMA473PopsGroupVec15], pch = 16, xlab = "Pop", bty = "n")
CommonNames473[which.min(x)]; CommonNames473[which.max(x)]
plot(y, col = Colors15[KMA473PopsGroupVec15], pch = 16, xlab = "Pop", bty = "n")
CommonNames473[which.min(y)]; CommonNames473[which.max(y)]
plot(z, col = Colors15[KMA473PopsGroupVec15], pch = 16, xlab = "Pop", bty = "n")
CommonNames473[which.min(z)]; CommonNames473[which.max(z)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Baseline with all 89 Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA47389Baseline <- CreateBaseline.GCL(sillyvec = KMA473Pops, loci = loci89, dir = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/BAYES/Baseline",
                                       basename = "KMA473Pops89Markers", type = "BAYES")
dput(x = KMA47389Baseline, file = "Objects/KMA47389Baseline.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Repeated 100% Proof Tests ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Groups15.nospace <- gsub(pattern = "/", replacement = "", x = gsub(pattern = " ", replacement = "", x = Groups15))
# dput(x = Groups15.nospace, file = "Objects/Groups15.nospace.txt")

KMA473Pops15FlatPrior <- Prior.GCL(groupvec = KMA473PopsGroupVec15, groupweights = rep(1 / 15, 15), minval = 0.01)
dput(x = KMA473Pops15FlatPrior, file = "Objects/KMA473Pops15FlatPrior.txt")

KMA473PopsInits <- MultiChainInits.GCL(npops = 473, nchains = 5, prop = 0.9)
dput(x = KMA473PopsInits, file = "Objects/KMA473PopsInits.txt")

Groups15.nospace
CommonNames473

KMA473Pops.named <- setNames(object = KMA473Pops, nm = CommonNames473)
dput(x = KMA473Pops.named, file = "Objects/KMA473Pops.named.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get sample sizes for Proof Tests
# Sample size for each of the 473 Pops (pooled)
KMA473Pops.SampleSize <- sapply(paste(KMA473Pops.named, ".gcl", sep = ''), function(x) get(x)$n)
str(KMA473Pops.SampleSize)

# Sample size for each of the 15 RGs
KMA473Pops.15RG.SampleSize <- setNames(object = sapply(seq(Groups15.nospace), function(RG) {sum(KMA473Pops.SampleSize[which(KMA473PopsGroupVec15 == RG)])} ), nm = Groups15.nospace)

# Sample size for proof tests
ProofTest100.SampleSize <- pmin(floor(KMA473Pops.15RG.SampleSize / 2), 200) # 200 is standard, unless RG has < 400, then it is (n of RG/2)
dput(x = ProofTest100.SampleSize, file = "Objects/ProofTest100.SampleSize.txt")

# Create sample size matrix for repeated proof test function
ProofTest100.SampleSize.Matrix <- matrix(data = 0, nrow = length(Groups15.nospace), ncol = length(Groups15.nospace))
colnames(ProofTest100.SampleSize.Matrix) <- Groups15.nospace
rownames(ProofTest100.SampleSize.Matrix) <- Groups15.nospace
diag(ProofTest100.SampleSize.Matrix) <- ProofTest100.SampleSize
dput(x = ProofTest100.SampleSize.Matrix, file = "Objects/ProofTest100.SampleSize.Matrix.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create Proof Tests
KMA473PopsGroups15Repeated100ProofTests <- paste(rep(Groups15.nospace, each = 5), 1:5, sep = '')
dput(x = KMA473PopsGroups15Repeated100ProofTests, file = "Objects/KMA473PopsGroups15Repeated100ProofTests.txt")

KMA473Pops.15RG.SampleSize
ProofTest100.SampleSize.Matrix

dir.create("BAYES/100% Proof Tests")
dir.create("BAYES/100% Proof Tests/loci89")

## All RGs
for(RG in seq(Groups15.nospace)) {
  for(Proof in KMA473PopsGroups15Repeated100ProofTests[(RG * 5 - 4):(RG * 5)]) {
    assign(x = paste(Proof, "Proof", sep = ""),
           value = ProofTest.GCL(sillyvec = KMA473Pops.named, loci = loci89, groupnames = Groups15.nospace, groupvec = KMA473PopsGroupVec15,
                                 samplesize = ProofTest100.SampleSize.Matrix[, RG], prefix = Proof, dir = "BAYES/100% Proof Tests/loci89", 
                                 prprtnl = TRUE, type = "BAYES", suffix = "", nreps = 40000, nchains = 5, priorvec = KMA473Pops15FlatPrior, 
                                 initmat = KMA473PopsInits, thin = c(1, 1, 100), switches = "F T F T T T F"))
  }
}; rm(RG, Proof); beep(8)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dput proof test objects
dir.create("100%ProofTests objects")
objects(pattern = "Proof$")
invisible(sapply(KMA473PopsGroups15Repeated100ProofTests, function(proof) {dput(x = get(paste(proof, "Proof", sep = "")), file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = ""))} ))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create BAYES.output folders
sapply(KMA473PopsGroups15Repeated100ProofTests, function(proof) {dir <- paste(getwd(), "/BAYES/100% Proof Tests/loci89/BAYES.output/", proof, sep = "")
                                                                 dir.create(dir) })


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize 100% proof tests ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create("Estimates objects")



KMA473PopsGroups15Repeated100ProofTests

KMA473PopsGroups15Repeated100ProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, maindir = "BAYES/100% Proof Tests/loci89/BAYES.output", 
                                                                                 mixvec = KMA473PopsGroups15Repeated100ProofTests, prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1,
                                                                                 PosteriorOutput = TRUE); beep(8)
dput(x = KMA473PopsGroups15Repeated100ProofTestsEstimates, file = "Estimates objects/KMA473PopsGroups15Repeated100ProofTestsEstimates.txt"); beep(4)
dput(x = KMA473PopsGroups15Repeated100ProofTestsEstimates$Stats, file = "Estimates objects/KMA473PopsGroups15Repeated100ProofTestsEstimatesStats.txt"); beep(4)

str(KMA473PopsGroups15Repeated100ProofTestsEstimates, max.level = 1)

## Check Gelman-Rubin factor
lapply(KMA473PopsGroups15Repeated100ProofTestsEstimates$Stats, function(Proof) {Proof[, "GR"]} )
sapply(KMA473PopsGroups15Repeated100ProofTestsEstimates$Stats, function(Proof) {table(Proof[, "GR"] > 1.2)} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~
# Basic Bars

library(gplots)
library(plotrix)

KMA473PopsGroups15Repeated100ProofTestsEstimatesStats <- dget(file = "Estimates objects/loci89/KMA473PopsGroups15Repeated100ProofTestsEstimatesStats.txt"); beep(4)

str(KMA473PopsGroups15Repeated100ProofTestsEstimatesStats)

Groups15.nospace

ProofTest100.SampleSize

ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(5.1, 4.4, 2.1, 1.1))

sapply(Groups15.nospace, function(RG) {
  par(xpd = FALSE)
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {KMA473PopsGroups15Repeated100ProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {KMA473PopsGroups15Repeated100ProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {KMA473PopsGroups15Repeated100ProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", ProofTest100.SampleSize[RG], ")", sep = ''), xlab = "",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 1), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  abline(h = 0.9, lwd = 2)
  text(x = ProofPlot[3, ], y = -0.065, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.8)
  mtext(text = "Reporting Group", side = 1, line = 3.5)
} )

#~~~~~~~~~~~~~~~~~~
# Master plot of all proof tests
# Making a BIG posterior out of the 5 repeats
# This is NOT statistically appropriate
KMA473PopsGroups15Repeated100ProofTestsEstimates <- dget(file = "Estimates objects/KMA473PopsGroups15Repeated100ProofTestsEstimates.txt"); beep(4)
str(KMA473PopsGroups15Repeated100ProofTestsEstimates)
str(KMA473PopsGroups15Repeated100ProofTestsEstimates[[1]])

for(RG in Groups15.nospace) {
  BAYESEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, maindir = "BAYES/100% Proof Tests/loci89/BAYES.output", 
                                                 mixvec = grep(pattern = RG, x = KMA473PopsGroups15Repeated100ProofTests, value = TRUE), prior = "", ext = "RGN", nchains = 5, 
                                                 burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE)
  assign(x = paste(RG, "100ProofTestEstimates", sep = ""), value = BAYESEstimates, pos = 1)
  rm(BAYESEstimates)
  dput(x = get(paste(RG, "100ProofTestEstimates", sep = "")), file = paste("Estimates objects/", RG, "100ProofTestEstimates", ".txt.", sep = ""))
}; beep(4)


KMA473PopsGroups15Repeated100ProofTestsEstimatesAllChainsStats <- matrix(data = NA, nrow = length(Groups15.nospace), ncol = 3, dimnames = list(Groups15.nospace, c("median", "5%", "95%")))

for(i in seq(Groups15.nospace)) {
  RG <- Groups15.nospace[i]
  RGEstimatesAllChains <- Reduce(f = "rbind", x = get(paste(RG, "100ProofTestEstimates", sep = ""))$Output)
  KMA473PopsGroups15Repeated100ProofTestsEstimatesAllChainsStats[i, ] <- quantile(x = RGEstimatesAllChains[, i], probs = c(0.5, 0.05, 0.95))
}; rm(i, RG, RGEstimatesAllChains)

KMA473PopsGroups15Repeated100ProofTestsEstimatesAllChainsStats

require(gplots)
par(mar = c(3.1, 4.4, 2.1, 10.1))

PCGroups15 <- c("West of Chignik", "Black Lake", "Chignik Lake", "Upper Station / Akalura", "Frazer", "Ayakulik", "Karluk", "Uganik", "Northwest Kodiak", "Afognak", "Eastside Kodiak", "Saltery", "Cook Inlet", "Prince William Sound", "South of Cape Suckling")
ProofPlot <- barplot2(height = KMA473PopsGroups15Repeated100ProofTestsEstimatesAllChainsStats[, "median"],
                      ci.l = KMA473PopsGroups15Repeated100ProofTestsEstimatesAllChainsStats[, "5%"],
                      ci.u = KMA473PopsGroups15Repeated100ProofTestsEstimatesAllChainsStats[, "95%"],
                      names.arg = NULL, col = Colors15, main = "", xlab = "", cex.axis = 1.5,
                      ylab = "", plot.ci = TRUE, beside = TRUE, ylim = c(0, 1), ci.lwd = 2, ci.color = "black")
mtext(text = "Reporting Group", side = 1, line = 1, cex = 1.4)
mtext(text = "Proportion Correctly Allocated", side = 2, line = 3, cex = 1.5)
segments(x0 = -1, y0 = 0.9, x1 = 18, y1 = 0.9, col = "black", lwd = 5)
legend(x = 18, y = 0.9, legend = PCGroups15, fill = Colors15, bty = "n", xpd = TRUE)


#~~~~~~~~~~~~~~~~~~
# Master plot of all proof tests
# 1 bar per repeat, so lots of bars, but better representation of posteriors

library(gplots)
library(plotrix)

KMA473PopsGroups15Repeated100ProofTestsEstimatesStats <- dget(file = "Estimates objects/loci89/KMA473PopsGroups15loci89Repeated100ProofTestsEstimatesStats.txt"); beep(4)

str(KMA473PopsGroups15Repeated100ProofTestsEstimatesStats)

Groups15.nospace

ProofTest100.SampleSize

KMA473PopsGroups15Repeated100ProofTestsEstimatesAllRepeatsStats <- t(sapply(KMA473PopsGroups15Repeated100ProofTests, function(Mix) {
  c(KMA473PopsGroups15Repeated100ProofTestsEstimatesStats[[Mix]][names(unlist(sapply(Groups15.nospace, function(RG) {grep(pattern = RG, x = Mix)} ))), "median"],
  KMA473PopsGroups15Repeated100ProofTestsEstimatesStats[[Mix]][names(unlist(sapply(Groups15.nospace, function(RG) {grep(pattern = RG, x = Mix)} ))), "5%"],
  KMA473PopsGroups15Repeated100ProofTestsEstimatesStats[[Mix]][names(unlist(sapply(Groups15.nospace, function(RG) {grep(pattern = RG, x = Mix)} ))), "95%"])
} ))

colnames(KMA473PopsGroups15Repeated100ProofTestsEstimatesAllRepeatsStats) <- c("median", "5%", "95%")

KMA473PopsGroups15Repeated100ProofTestsEstimatesAllRepeatsStats

# ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(3.1, 4.4, 2.1, 10.1))

PCGroups15 <- c("West of Chignik", "Black Lake", "Chignik Lake", "Upper Station / Akalura", "Frazer", "Ayakulik", "Karluk", "Uganik", "Northwest Kodiak", "Afognak", "Eastside Kodiak", "Saltery", "Cook Inlet", "Prince William Sound", "South of Cape Suckling")
ProofPlot <- barplot2(height = KMA473PopsGroups15Repeated100ProofTestsEstimatesAllRepeatsStats[, "median"],
                      ci.l = KMA473PopsGroups15Repeated100ProofTestsEstimatesAllRepeatsStats[, "5%"],
                      ci.u = KMA473PopsGroups15Repeated100ProofTestsEstimatesAllRepeatsStats[, "95%"],
                      names.arg = NULL, col = rep(Colors15, each = 5), main = "", xlab = "", cex.axis = 1.5,
                      ylab = "", plot.ci = TRUE, beside = TRUE, ylim = c(0, 1), ci.lwd = 2, ci.color = "black", space = c(0.3, rep(c(rep(0.3, 4), 1), 14), rep(0.3, 4)))
mtext(text = "Reporting Group", side = 1, line = 1, cex = 1.4)
mtext(text = "Proportion Correctly Allocated", side = 2, line = 3, cex = 1.5)
segments(x0 = 0, y0 = 0.9, x1 = max(ProofPlot) + 0.5, y1 = 0.9, col = "black", lwd = 5)
legend(x = max(ProofPlot) + 2, y = 0.9, legend = PCGroups15, fill = Colors15, bty = "n", xpd = TRUE)

#~~~~~~~~~~~~~~~~~~
# Violinplot

require(vioplot)

substrRight <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

KMAGroups15Repeated100ProofTestsEstimates <- dget(file = "Estimates objects/KMAGroups15Repeated100ProofTestsEstimates.txt"); beep(4)


KMAGroups15Repeated100ProofTests
Groups15.nospace
RG <- "CookInlet"
proofs <- grep(pattern = RG, x = KMAGroups15Repeated100ProofTests)

lapply(KMAGroups15Repeated100ProofTestsEstimates$Stats[proofs], function(Proof) {Proof[, "GR"]} )
lapply(KMAGroups15Repeated100ProofTestsEstimates$Stats[proofs], function(Proof) {table(Proof[, "GR"] > 1.2)} )

par(mar = c(5.1, 4.6, 3.6, 6.6))


plot(y = t(sapply(KMAGroups15Repeated100ProofTestsEstimates$Stats[proofs], function(Proof) {Proof[, "mean"]} )), x = seq(5 * length(Groups15.nospace)), pch =16, col = rep(Colors15, each = 5), ylab = "", ylim = c(0, 1), xlab = "", axes = FALSE, main = "", cex = 2)
invisible(sapply(KMAGroups15Repeated100ProofTests[proofs], function(Proof) {sapply(1:15, function(i) {vioplot2(KMAGroups15Repeated100ProofTestsEstimates$Output[[Proof]][, i], at = as.numeric(substrRight(x = Proof, n = 1)) + ((i - 1) * 5), horizontal = FALSE, col = Colors15[i], border = TRUE, drawRect = FALSE, rectCol = Colors15[i], add = TRUE, wex = 1, lwd = 1)} )} ))
arrows(x0 = seq(5 * length(Groups15.nospace)), x1 = seq(5 * length(Groups15.nospace)), y0 = t(sapply(KMAGroups15Repeated100ProofTestsEstimates$Stats[proofs], function(Proof) {Proof[, "5%"]} )), y1 = t(sapply(KMAGroups15Repeated100ProofTestsEstimates$Stats[proofs], function(Proof) {Proof[, "95%"]} )), angle = 90, code = 3, length = 0.1, lwd = 2)
points(y = t(sapply(KMAGroups15Repeated100ProofTestsEstimates$Stats[proofs], function(Proof) {Proof[, "mean"]} )), x = seq(5 * length(Groups15.nospace)), cex = 1.5, pch = 21, col = "white", bg = rep(Colors15, each = 5), lwd = 2)
axis(side = 2, lwd = 3, cex.axis = 1.5, pos = 0)
mtext(text = "Proportion", side = 2, line = 1, cex = 2)
mtext(text = paste(RG, " (n = ", samplesz_proof_final[RG], ")", sep = ''), side = 3, cex = 2.5)
text(x = seq(from = 3, by = 5, length.out = 15) - 0.35, y = 0, labels = Groups15Short, srt = 45, pos = 1, offset = 2.5, xpd = TRUE, cex = 1.6)
segments(x0 = 0, y0 = 0, y1 = 0, x1 = 5 * length(Groups15.nospace), lwd = 3)
segments(x0 = 0, y0 = 0.9, y1 = 0.9, x1 = 5 * length(Groups15.nospace), lwd = 3)
#legend(x = 70, y = 0.9, xpd = TRUE, legend = Groups15Short, fill = Colors15, bty = "n")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Table 100% Proof Tests ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goal is to follow format of SEDM Table 42

# PCFisheryProofTestScenarioNames <- c("July Uyak", "July Alitak", "August Karluk", "June Cape Igvak", "June Ayakulik", "August Uganik", "Flat")
# dput(x = PCFisheryProofTestScenarioNames, file = "Objects/PCFisheryProofTestScenarioNames.txt")
# FisheryProofTestScenarioNamesShort <- c("Uyak", "Alitak", "Karluk", "Igvak", "Ayakulik", "Uganik", "Flat")
# dput(x = FisheryProofTestScenarioNamesShort, file = "Objects/FisheryProofTestScenarioNamesShort.txt")

# dir.create("Estimates tables")


# Inputs
one.hundred.percent <- TRUE
loci <- "loci89"
i <- 1  # 1:15 for 100% and 1:7 for fishery scenarios
table.file <- "Estimates tables/KMA473PopsGroups15RepeatedProofTests.xlsx"  # "V:/DOC/FINAL/SOCKEYE/Kodiak/KMA 2014-2016/KMA 2015 Baseline FMS/FMS 16-XX KMA Sockeye 2015 Baseline Tables.xlsx"


ProofTableMaker <- function(one.hundred.percent = TRUE, loci = "loci89", i, table.file = "Estimates appendix tables/Test.xlsx"){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function makes a 3 page table for repeated proof tests, for both 100% and fishery scenarios
  
  # one.hundred.percent - boolean - TRUE for 100% Proof Tests, FALSE for Fishery Scenario Tests
  # loci - character - "loci89", "loci46", "loci22", "loci24", "loci24Kodiak"
  # i - numeric - 1:15 for 100% reporting groups, 1:7 for fishery scenarios
  # table.file - character - "Estimates tables/KMA473PopsGroups15RepeatedProofTests.xlsx"
  
  # Output is to append a sheet to the "table.file" provided
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Helper functions
  ProofTableStats <- function(StatsObject) {
    ProofTableColnames <- c("median", "5%", "95%", "P=0", "mean", "sd")
    x <- formatC(x = t(t(StatsObject[, ProofTableColnames]) * c(rep(100, 3), 1, rep(100, 2))), digits = 2, format = "f")
    x[, -4] <- formatC(x = as.numeric(x[, -4]), digits = 1, format = "f")
    return(x)
  }
  
  AverageOfReplicates <- function(StatsObject) {
    x <- apply(sapply(StatsObject, function(rplct) {rplct[, "mean"]} ), 1, mean) * 100
    x <- formatC(x = x, digits = 1, format = "f")
    return(x)
  }
  
  
  # If else to sort between 100% Proof Test and Fishery Scenario Tests
  if(one.hundred.percent) {
    PCscenario <- PCGroups15[i]
    scenario <- Groups15.nospace[i]
    scenario.short <- Groups15Short[i]
    n <- ProofTest100.SampleSize[i]
    proofTableCaption <- paste("Table X.-Regional estimates of stock composition for 5 replicates of a 100% proof test of the ", PCscenario," reporting group included as part of the coast-wide sockeye salmon genetic baseline with ", as.numeric(gsub("[[:alpha:]]","",loci)), " loci. Each replicate was a sample of ", n," individuals removed from the genetic baseline. Estimates include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.", sep = '')  # used to estimate stock compositions of Kodiak Management Area harvests of sockeye salmon in 2014-2016
    proportions <- ProofTest100.SampleSize.Matrix[, i] / max(ProofTest100.SampleSize.Matrix[, i])
    stats.object <- dget(file = paste("Estimates objects/", loci, "/KMA473PopsGroups15", loci, "Repeated100ProofTestsEstimatesStats.txt", sep = ''))
  } else {
    PCscenario <- paste("Hypothetical", PCFisheryProofTestScenarioNames[i])
    scenario <- FisheryProofTestScenarioNames[i]
    scenario.short <- FisheryProofTestScenarioNamesShort[i]
    n <- rowSums(MixtureProofTest.SampleSize)[i]
    proofTableCaption <- paste("Table X.-Regional estimates of stock composition for 5 replicates of a ", PCscenario," sockeye salmon fishery in the Kodiak Management Area based on the coast-wide sockeye salmon genetic baseline with ", as.numeric(gsub("[[:alpha:]]","",loci)), " loci. Each replicate was a sample of ", n," individuals removed from the genetic baseline according to hypothetical proportions. Estimates include median, 90% credibility interval, the probability that the group estimate is equal to zero (P=0), mean and SD.", sep = '')
    MixtureProofTestProportions <- dget(file = "Objects/MixtureProofTestProportions.txt")
    proportions <- MixtureProofTestProportions[i, ]
    stats.object <- dget(file = paste("Estimates objects/", loci, "/KMA473PopsGroups15", loci, "RepeatedMixProofTestsEstimatesStats.txt", sep = ''))
  }
  
  disclaimer <- "Note: Stock composition estimates may not sum to 100% due to rounding error."
  proportions <- formatC(x = proportions * 100, digits = 1, format = "f")
  stats.object <- stats.object[grep(pattern = scenario, x = names(stats.object))]
  scenario.short <- gsub(pattern = "/", replacement = " ", x = scenario.short)
  
  # Building Table Page 1
  proofTablePage1 <- cbind(proofTableCaption, '', '', '', '', '' , '', '', '', '', '', '', '', '', '', '')
  proofTablePage1 <- rbind(proofTablePage1, cbind('', '', '', paste(PCscenario, "Replicate 1", sep = ' '), '', '', '', '', '', '', paste(PCscenario, "Replicate 2", sep = ' '), '', '', '', '', ''))
  proofTablePage1 <- rbind(proofTablePage1, cbind('', "True", '', '', "90% CI", '', '', '', '', '', '', "90% CI", '', '', '', ''))
  proofTablePage1 <- rbind(proofTablePage1, cbind("Reporting Group", "Percent", '', "Median", "5%", "95%", "P=0", "Mean", "SD", '', "Median", "5%", "95%", "P=0", "Mean", "SD"))
  proofTablePage1 <- rbind(proofTablePage1, cbind(PCGroups15, proportions, rep('', 15), ProofTableStats(StatsObject = stats.object[[1]]), rep('', 15), ProofTableStats(StatsObject = stats.object[[2]])))
  proofTablePage1 <- rbind(proofTablePage1, cbind(disclaimer, '', '', '', '', '' , '', '', '', '', '', '', '', '', '', ''))
  
  # Building Table Page 2
  proofTablePage2 <- cbind("Table X. Page 2 of 3.", '', '', '', '', '' , '', '', '', '', '', '', '', '', '', '')
  proofTablePage2 <- rbind(proofTablePage2, cbind('', '', '', paste(PCscenario, "Replicate 3", sep = ' '), '', '', '', '', '', '', paste(PCscenario, "Replicate 4", sep = ' '), '', '', '', '', ''))
  proofTablePage2 <- rbind(proofTablePage2, cbind('', "True", '', '', "90% CI", '', '', '', '', '', '', "90% CI", '', '', '', ''))
  proofTablePage2 <- rbind(proofTablePage2, cbind("Reporting Group", "Percent", '', "Median", "5%", "95%", "P=0", "Mean", "SD", '', "Median", "5%", "95%", "P=0", "Mean", "SD"))
  proofTablePage2 <- rbind(proofTablePage2, cbind(PCGroups15, proportions, rep('', 15), ProofTableStats(StatsObject = stats.object[[3]]), rep('', 15), ProofTableStats(StatsObject = stats.object[[4]])))
  proofTablePage2 <- rbind(proofTablePage2, cbind(disclaimer, '', '', '', '', '' , '', '', '', '', '', '', '', '', '', ''))
  
  # Building Table Page 3
  proofTablePage3 <- cbind("Table X. Page 3 of 3.", '', '', '', '', '' , '', '', '', '', '')
  proofTablePage3 <- rbind(proofTablePage3, cbind('', '', '', paste(PCscenario, "Replicate 5", sep = ' '), '', '', '', '', '', '', ''))
  proofTablePage3 <- rbind(proofTablePage3, cbind('', "True", '', '', "90% CI", '', '', '', '', '', "Average of 5 Replicates"))
  proofTablePage3 <- rbind(proofTablePage3, cbind("Reporting Group", "Percent", '', "Median", "5%", "95%", "P=0", "Mean", "SD", '', ''))
  proofTablePage3 <- rbind(proofTablePage3, cbind(PCGroups15, proportions, rep('', 15), ProofTableStats(StatsObject = stats.object[[5]]), rep('', 15), AverageOfReplicates(StatsObject = stats.object)))
  proofTablePage3 <- rbind(proofTablePage3, cbind(disclaimer, '', '', '', '', '' , '', '', '', '', ''))
  
  require(xlsx)
  suppressWarnings(write.xlsx(x = proofTablePage1, file = table.file, sheetName = paste(loci, scenario.short, "page1", sep = "_"), col.names = FALSE, row.names = FALSE, append = TRUE))
  suppressWarnings(write.xlsx(x = proofTablePage2, file = table.file, sheetName = paste(loci, scenario.short, "page2", sep = "_"), col.names = FALSE, row.names = FALSE, append = TRUE))
  suppressWarnings(write.xlsx(x = proofTablePage3, file = table.file, sheetName = paste(loci, scenario.short, "page3", sep = "_"), col.names = FALSE, row.names = FALSE, append = TRUE))
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Make Proof Tables!!!
# dir.create("Estimates appendix tables")

# 100%; loci89
if (.Platform$OS.type == "windows") flush.console()
pb <- txtProgressBar(min = 0, max = length(Groups15), style = 3)
for(i in seq(Groups15)){
  setTxtProgressBar(pb = pb, value = i)
  ProofTableMaker(one.hundred.percent = TRUE, loci = "loci89", i = i, table.file = "Estimates appendix tables/KMA473PopsGroups15Repeated100ProofTestsloci89Appendices.xlsx")
};# beep(4)

# Fishery Scenario; loci89
for(i in seq(FisheryProofTestScenarioNames)){
  ProofTableMaker(one.hundred.percent = FALSE, loci = "loci89", i = i, table.file = "Estimates appendix tables/KMA473PopsGroups15RepeatedFisheryTestsloci89Appendices.xlsx")
};# beep(4)


# 100%; loci46
for(i in seq(Groups15)){
  ProofTableMaker(one.hundred.percent = TRUE, loci = "loci46", i = i, table.file = "Estimates appendix tables/KMA473PopsGroups15Repeated100ProofTestsloci46Appendices.xlsx")
};# beep(4)

# Fishery Scenario; loci46
for(i in seq(FisheryProofTestScenarioNames)){
  ProofTableMaker(one.hundred.percent = FALSE, loci = "loci46", i = i, table.file = "Estimates appendix tables/KMA473PopsGroups15RepeatedFisheryTestsloci46Appendices.xlsx")
}; beep(4)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Fishery Scenario Tests ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA473Pops15FlatPrior
KMA473PopsInits
Groups15Short
CommonNames473
KMA473Pops

# Population sample sizes
KMA473Pops.SampleSize <- sapply(paste(KMA473Pops.named, ".gcl", sep = ''), function(x) get(x)$n)

# RG samples sizes
KMA473Pops.15RG.SampleSize <- setNames(object = sapply(seq(Groups15.nospace), function(RG) {sum(KMA473Pops.SampleSize[which(KMA473PopsGroupVec15 == RG)])} ), nm = Groups15.nospace)

# Proof test sample sizes, 100%
ProofTest100.SampleSize

## Mixture proof test
# Get fishery scenarios from managers
MixtureProofTestProportions <- read.table(file = "FisheryProofTestScenarios.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = Groups15Short)
str(MixtureProofTestProportions)
MixtureProofTestProportions <- t(data.matrix(MixtureProofTestProportions[, -1]))
dput(x = MixtureProofTestProportions, file = "Objects/MixtureProofTestProportions.txt")
dput(x = rownames(MixtureProofTestProportions), file = "Objects/FisheryProofTestScenarioNames.txt")

# Sample sizes are adjusted to get whole numbers
MixtureProofTest.SampleSize <- t(apply(MixtureProofTestProportions, 1, function(scenario) {floor(scenario * 400)} ))
apply(MixtureProofTest.SampleSize, 1, sum)
dput(x = MixtureProofTest.SampleSize, file = "Objects/MixtureProofTest.SampleSize.txt")

# Create directories
dir.create("BAYES/Mixture Proof Tests")
dir.create("BAYES/Mixture Proof Tests/loci89")
invisible(sapply(c("baseline", "control", "mixture", "output"), function(fldr) {dir.create(path = paste("BAYES/Mixture Proof Tests/loci89/BAYES.", fldr, sep = ""))} ))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create mixture proof tests
sapply(rownames(MixtureProofTest.SampleSize), function(Scenario) {sapply(1:5, function(Rpt) {
  assign(x = paste(Scenario, Rpt, "Proof", sep = ''), 
         value = ProofTest.GCL(sillyvec = KMA473Pops, loci = loci89, groupnames = Groups15Short, groupvec = KMA473PopsGroupVec15,
                               samplesize = MixtureProofTest.SampleSize[Scenario, ], prefix = paste(Scenario, Rpt, sep = ''), 
                               dir = "BAYES/Mixture Proof Tests/loci89", prprtnl = TRUE, type = "BAYES",
                               suffix = '', nreps = 40000, nchains = 5, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits, 
                               thin = c(1, 1, 100), switches = "F T F T T T F"), pos = 1)
} )} ); beep(8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dput proof test objects
objects(pattern = "Proof$")
KMA473PopsGroups15RepeatedMixProofTests <- paste(rep(rownames(MixtureProofTest.SampleSize), each = 5), 1:5, sep = '')
dput(x = KMA473PopsGroups15RepeatedMixProofTests, file = "Objects/KMA473PopsGroups15RepeatedMixProofTests.txt")

dir.create("MixtureProofTests objects")
invisible(sapply(KMA473PopsGroups15RepeatedMixProofTests, function(proof) {dput(x = get(paste(proof, "Proof", sep = "")), file = paste("MixtureProofTests objects/", proof, "Proof.txt", sep = ""))} ))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create BAYES.output folders
invisible(sapply(KMA473PopsGroups15RepeatedMixProofTests, function(proof) {dir <- paste(getwd(), "/BAYES/Mixture Proof Tests/loci89/BAYES.output/", proof, sep = "")
                                                                           dir.create(dir) }))

#~~~~~~~~~~~~~~~~~~
## RUN BAYES ##
#~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Summarize Fishery Scenario Tests ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Are they all done?
sapply(KMA473PopsGroups15RepeatedMixProofTests, function(Mix) {length(grep(pattern = Mix, x = list.files(path = "BAYES/Mixture Proof Tests/loci89/BAYES.output", pattern = ".RGN", recursive = TRUE), value = TRUE))} )



KMA473PopsGroups15RepeatedMixProofTests

# Just look at flat resutls (attractor/deflector RGs)
KMA473PopsGroups15RepeatedMixProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15Short), groupnames = Groups15Short, maindir = "BAYES/Mixture Proof Tests/loci89/BAYES.output", 
                                                                                 mixvec = KMA473PopsGroups15RepeatedMixProofTests, prior = "", ext = "RGN", 
                                                                                 nchains = 5, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE); beep(8)

dput(x = KMA473PopsGroups15RepeatedMixProofTestsEstimates, file = "Estimates objects/KMA473PopsGroups15RepeatedMixProofTestsEstimates.txt"); beep(4)
dput(x = KMA473PopsGroups15RepeatedMixProofTestsEstimates$Stats, file = "Estimates objects/KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats.txt"); beep(4)

str(KMA473PopsGroups15RepeatedMixProofTestsEstimates, max.level = 1)

## Check Gelman-Rubin factor
lapply(KMA473PopsGroups15RepeatedMixProofTestsEstimates$Stats, function(Proof) {Proof[, "GR"]} )
sapply(KMA473PopsGroups15RepeatedMixProofTestsEstimates$Stats, function(Proof) {table(Proof[, "GR"] > 1.2)} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~
# Basic Bars

library(gplots)
library(plotrix)

KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci89/KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats.txt"); beep(4)

str(KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats)

Groups15.nospace

MixtureProofTest.SampleSize.Final <- round(apply(MixtureProofTest.SampleSize, 1, sum), 0 )

MixtureProofTestProportions <- read.table(file = "FisheryProofTestScenarios.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = Groups15Short)
MixtureProofTestProportions <- MixtureProofTestProportions[, -1]
max(MixtureProofTestProportions)
ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(xpd = FALSE)
FisheryProofTestScenarioNames <- dget(file = "Objects/FisheryProofTestScenarioNames.txt")
# RG = "flat"
for(RG in FisheryProofTestScenarioNames) {
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", MixtureProofTest.SampleSize.Final[RG], ")", sep = ''), xlab = "Reporting Group",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 0.6), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  #abline(h = 1/15, lwd = 4)

  for(i in 1:15) {
    segments(x0 = ProofPlot[1, i] - 0.5, x1 = ProofPlot[5, i] + 0.5, y0 = MixtureProofTestProportions[i, RG], y1 = MixtureProofTestProportions[i, RG], lwd = 3)
  }
  
  legend(x = max(ProofPlot), y = 0.5, legend = PCGroups15, fill = Colors15, bty = "n", xpd = TRUE)
  
  text(x = ProofPlot[3, ], y = -0.03, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.7)
}

#~~~~~~~~~~~~~~~~~~
# Violinplot

require(vioplot)

substrRight <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

KMAGroups15RepeatedFlatMixProofTestsEstimates <- dget(file = "Estimates objects/KMAGroups15RepeatedFlatMixProofTestsEstimates.txt"); beep(4)


KMA473PopsGroups15RepeatedMixProofTests
Groups15.nospace
RG <- "flat"
proofs <- grep(pattern = RG, x = KMA473PopsGroups15RepeatedMixProofTests)

lapply(KMAGroups15RepeatedFlatMixProofTestsEstimates$Stats, function(Proof) {Proof[, "GR"]} )
lapply(KMAGroups15RepeatedFlatMixProofTestsEstimates$Stats, function(Proof) {table(Proof[, "GR"] > 1.2)} )

par(mar = c(5.1, 4.6, 3.6, 6.6))


plot(y = t(sapply(KMAGroups15RepeatedFlatMixProofTestsEstimates$Stats, function(Proof) {Proof[, "mean"]} )), x = seq(5 * length(Groups15.nospace)), pch =16, col = rep(Colors15, each = 5), ylab = "", ylim = c(0, 1), xlab = "", axes = FALSE, main = "", cex = 2)
invisible(sapply(KMA473PopsGroups15RepeatedMixProofTests[proofs], function(Proof) {sapply(1:15, function(i) {vioplot2(KMAGroups15RepeatedFlatMixProofTestsEstimates$Output[[Proof]][, i], at = as.numeric(substrRight(x = Proof, n = 1)) + ((i - 1) * 5), horizontal = FALSE, col = Colors15[i], border = TRUE, drawRect = FALSE, rectCol = Colors15[i], add = TRUE, wex = 1, lwd = 1)} )} ))
arrows(x0 = seq(5 * length(Groups15.nospace)), x1 = seq(5 * length(Groups15.nospace)), y0 = t(sapply(KMAGroups15RepeatedFlatMixProofTestsEstimates$Stats, function(Proof) {Proof[, "5%"]} )), y1 = t(sapply(KMAGroups15RepeatedFlatMixProofTestsEstimates$Stats, function(Proof) {Proof[, "95%"]} )), angle = 90, code = 3, length = 0.1, lwd = 2)
points(y = t(sapply(KMAGroups15RepeatedFlatMixProofTestsEstimates$Stats, function(Proof) {Proof[, "mean"]} )), x = seq(5 * length(Groups15.nospace)), cex = 1.5, pch = 21, col = "white", bg = rep(Colors15, each = 5), lwd = 2)
axis(side = 2, lwd = 3, cex.axis = 1.5, pos = 0)
mtext(text = "Proportion", side = 2, line = 1, cex = 2)
mtext(text = paste(RG, " (n = ", MixPTsamplesize_final[RG], ")", sep = ''), side = 3, cex = 2.5)
text(x = seq(from = 3, by = 5, length.out = 15) - 0.35, y = 0, labels = Groups15Short, srt = 45, pos = 1, offset = 2.5, xpd = TRUE, cex = 1.6)
segments(x0 = 0, y0 = 0, y1 = 0, x1 = 5 * length(Groups15.nospace), lwd = 3)
segments(x0 = 0, y0 = 1 / length(Groups15Short), y1 = 1 / length(Groups15Short), x1 = 5 * length(Groups15.nospace), lwd = 3)
#legend(x = 70, y = 0.9, xpd = TRUE, legend = Groups15Short, fill = Colors15, bty = "n")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Reduce to 24 SNPs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("V:/WORK/Kyle/R Scripts/VarCompKS.GCL.R")

# Varcomps for all 473 Pops, 15RGs, 96 loci
# KMA473Pops15Group96lociVarComps <- VarCompKS.GCL(sillyvec = KMA473Pops, loci = loci96, groupvec = KMA473PopsGroupVec15, fstatfile = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/Tables/fstatfile.dat"); beep(8)
# dput(x = KMA473Pops15Group96lociVarComps, file = "Objects/KMA473Pops15Group96lociVarComps.txt")

# Varcomps for all 473 Pops, 15RGs, 103 possible loci (all linked combinations)
KMA473Pops15Group103lociVarComps <- VarCompKS.GCL(sillyvec = KMA473Pops, loci = loci103, groupvec = KMA473PopsGroupVec15, dir = "Tables")
dput(x = KMA473Pops15Group103lociVarComps, file = "Variance Components/KMA473Pops15Group103lociVarComps.txt")
KMA473Pops15Group103lociVarComps <- dget(file = "Variance Components/KMA473Pops15Group103lociVarComps.txt")

# Varcomps for all 57 Pops (just Kodiak), 9RGs, 103 loci
Groups15
Kodiak57Pops
Kodiak57PopsGroupVec9

Kodiak57Pops9Group103lociVarComps <- VarCompKS.GCL(sillyvec = Kodiak57Pops, loci = loci103, groupvec = Kodiak57PopsGroupVec9, dir = "Tables"); beep(8)
dput(x = Kodiak57Pops9Group103lociVarComps, file = "Variance Components/Kodiak57Pops9Group103lociVarComps.txt")

Kodiak57Pops9Group103lociVarComps
round(Kodiak57Pops9Group103lociVarComps[order(Kodiak57Pops9Group103lociVarComps[, "P"]), ], 3)
round(Kodiak57Pops9Group103lociVarComps[order(Kodiak57Pops9Group103lociVarComps[, "S"]), ], 3)
round(Kodiak57Pops9Group103lociVarComps[order(Kodiak57Pops9Group103lociVarComps[, "S"] + Kodiak57Pops9Group103lociVarComps[, "P"]), ], 3)

# Varcomps for all 11 Pops (Frazer/Ayakulik), 2RGs, 103 loci
FrazerAyakulik12Pops <- KMA473Pops[which(KMA473PopsGroupVec15 >= 5 & KMA473PopsGroupVec15 <= 6)]
FrazerAyakulik12PopsGroupVec2 <- KMA473PopsGroupVec15[which(KMA473PopsGroupVec15 >= 5 & KMA473PopsGroupVec15 <= 6)] - 4

FrazerAyakulik12Pops2Group103lociVarComps <- VarCompKS.GCL(sillyvec = FrazerAyakulik12Pops, loci = loci103, groupvec = FrazerAyakulik12PopsGroupVec2, dir = "Tables"); beep(8)
dput(x = FrazerAyakulik12Pops2Group103lociVarComps, file = "Variance Components/FrazerAyakulik12Pops2Group103lociVarComps.txt")

FrazerAyakulik12Pops2Group103lociVarComps
round(FrazerAyakulik12Pops2Group103lociVarComps[order(FrazerAyakulik12Pops2Group103lociVarComps[, "P"]), ], 3)
round(FrazerAyakulik12Pops2Group103lociVarComps[order(FrazerAyakulik12Pops2Group103lociVarComps[, "S"]), ], 3)
round(FrazerAyakulik12Pops2Group103lociVarComps[order(FrazerAyakulik12Pops2Group103lociVarComps[, "S"] + FrazerAyakulik12Pops2Group103lociVarComps[, "P"]), ], 3)

# Varcomps for all 473 Pops, NO RGs, 103 loci (i.e. 2 level hierarchy)
sillyvec = KMA473Pops; loci = loci103; fstatfile = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/Tables/473pops103loci_fstatfile.dat"

require("hierfstat")

names(sillyvec) <- NULL

nloci=length(loci)
nalleles=LocusControl$nalleles[loci]
ploidy=LocusControl$ploidy[loci] 
alleles=LocusControl$alleles[loci]

my.gcl=sapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)},simplify=FALSE)

n=sapply(sillyvec,function(silly){my.gcl[[silly]]$n})

dat0=read.fstat.data(fstatfile)

dat=data.frame(dat0)

dimnames(dat)[[2]]=c("Pop",loci)

VarCompByLocus=array(0,c(nloci,3),dimnames=list(loci,c("S","I","G")))

cat("\nCalculating Variance Components for each locus\n", sep = '')
if (.Platform$OS.type == "windows") flush.console()
pb <- txtProgressBar(min = 0, max = nloci, style = 3)

for(locus in loci){
  setTxtProgressBar(pb = pb, value = which(loci == locus))
  
  data=data.frame(dat[,1],dat[,locus])
  
  vc=varcomp(data,diploid=ploidy[locus]>1)$overall
  
  if(!sum(is.na(vc))){
    
    VarCompByLocus[locus,1:length(vc)]=vc
    
  }
  
}

KMA473Pops0Group103lociVarComps <- VarCompByLocus
dput(x = KMA473Pops0Group103lociVarComps, file = "Variance Components/KMA473Pops0Group103lociVarComps.txt")

rm(sillyvec, loci, fstatfile, nloci, nalleles, ploidy, alleles, my.gcl, n, dat0, dat, VarCompByLocus, locus)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Table results
dir.create("Variance Components")
require(xlsx)
write.xlsx(x = KMA473Pops15Group103lociVarComps, file = "Variance Components/KMA473Pops15Group103lociVarComps.xlsx")
write.xlsx(x = Kodiak57Pops9Group103lociVarComps, file = "Variance Components/Kodiak57Pops9Group103lociVarComps.xlsx")
write.xlsx(x = FrazerAyakulik12Pops2Group103lociVarComps, file = "Variance Components/FrazerAyakulik12Pops2Group103lociVarComps.xlsx")
write.xlsx(x = KMA473Pops0Group103lociVarComps, file = "Variance Components/KMA473Pops0Group103lociVarComps.xlsx")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Trying a 4 level hierarchy: 1) Kodiak/nonKodiak, 2) 15RGs, 3) Pops, 4) Indiv in Pops, 5) Alleles in individuals
### 103 Loci
## Locus Control
LocusControl <- dget(file = "Objects/LocusControl103loci.txt")
loci103 <- dget(file = "Objects/loci103.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get Populations
KMA473Pops <- dget(file = "Objects/KMA473Pops.txt")

require(beepr)
sapply(KMA473Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPops103loci/", silly, ".txt", sep = "")), pos = 1)}); beep(2)
length(objects(pattern = "\\.gcl"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get objects
KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("OLD", "OriginalLocusControl.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUTS
sillyvec <- KMA473Pops
loci <- loci103
groupvec <- KMA473PopsGroupVec15
fstatfile = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/Tables/473pops103loci15RGfstatfile.dat"

# FUNCTION
require("hierfstat")

names(sillyvec) <- NULL

nloci=length(loci)
nalleles=LocusControl$nalleles[loci]
ploidy=LocusControl$ploidy[loci] 
alleles=LocusControl$alleles[loci]

my.gcl=sapply(sillyvec,function(silly){get(paste(silly,".gcl",sep=""),pos=1)},simplify=FALSE)

n=sapply(sillyvec,function(silly){my.gcl[[silly]]$n})

# Read in fstatfile
dat0=read.fstat.data(fstatfile)
dat=data.frame(,rep(groupvec,n),dat0)

dimnames(dat)[[2]]=c("Island", "RG", "Pop", loci)

G=max(groupvec)

VarCompByLocus=array(0,c(nloci,4),dimnames=list(loci,c("P","S","SS","I","G")))

for(locus in loci){
  
  data=data.frame(dat[,1:2],dat[,locus])
  
  vc=varcomp(data,diploid=ploidy[locus]>1)$overall
  
  if(!sum(is.na(vc))){
    
    VarCompByLocus[locus,1:length(vc)]=vc
    
  }
  
}

VarCompByLocus

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Table results
KMA473Pops103loci4levelVarCompByLocus_list <- lapply(seq(list.files(path = "Variance Components", pattern = "4level")), function(i) {dget(file = paste(getwd(), "/Variance Components/KMA473Pops103loci4levelVarCompByLocus_", i, ".txt", sep = ""))} )

KMA473Pops103loci4levelVarCompByLocus <- do.call(what = "rbind", KMA473Pops103loci4levelVarCompByLocus_list)

str(KMA473Pops103loci4levelVarCompByLocus)

dput(x = KMA473Pops103loci4levelVarCompByLocus, file = "Variance Components/KMA473Pops103loci4levelVarCompByLocus.txt")
write.xlsx(x = KMA473Pops103loci4levelVarCompByLocus, file = "Variance Components/KMA473Pops103loci4levelVarCompByLocus.xlsx")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# I picked 24 markers based on Variance componenets
# P -> Ayakulik/Frazer 2RG 11Pops
# P -> Kodiak 9RG 50Pops
# P -> KMA 15RG 473Pops
# S -> KMA 0RG 473Pops
# P -> KMA 2 Regions (Kodiak/non-Kodiak) 15RG 473Pops
# S -> KMA 2 Regions (Kodiak/non-Kodiak) 15RG 473Pops

# Also checked out Ho and Fst
# Allele frequency plots were examined
# As were Fluidigm plots to make sure these are "easy" to score
# Kept two sets of linked markers "One_CO1.One_Cytb_26" & "One_Tf_ex10-750.One_Tf_ex3-182"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Attempt 1
loci22 <- dget(file = "OLD 746 Collections/Objects/loci22.txt")  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Attempt 2
# No linked markers
loci24 <- dget(file = "OLD 746 Collections/Objects/loci24.txt")  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Attempt 3
# No linked markers, Kodiak heavy
loci24Kodiak <- dget(file = "OLD 746 Collections/Objects/loci24Kodiak.txt")  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Attempt 4
# 48 markers
loci46 <- dget(file = "OLD 746 Collections/Objects/loci46.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Likelihood Profiles 24 Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Want to confirm that the likelihood profiles look good before I spend a lot of time on repeated proof tests

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Testing Kyle's custom version with progress bar
# Appears to work! Nice to see how long things are taking as it goes.
# Kodiak 50 Populations 9 RGs
rm(list = ls(all = TRUE))
setwd("V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline")
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")
source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")
source("V:/DATA/R_GEN/JJs GCL/LeaveOneOutlIkeProfile.GCL.R")

Kodiak57Pops <- dget(file = "Objects/Kodiak57Pops.txt")
Kodiak57PopsGroupVec9 <- dget(file = "Objects/Kodiak57PopsGroupVec9.txt")
Groups9.nospace.Kodiak <- dget(file = "Objects/Groups9.nospace.Kodiak.txt")
Groups9Short.Kodiak <- dget(file = "Objects/Groups9Short.Kodiak.txt")
Colors9 <- dget(file = "Objects/Colors9.txt")

LocusControl <- dget(file = "Objects/LocusControl103.txt")
loci103 <- dget(file = "Objects/loci103.txt")
invisible(sapply(Kodiak57Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPopsloci103/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)


# loci22 ~~~~~~~~~~
loci22 <- dget(file = "OLD 746 Collections/Objects/loci22.txt")  
Kodiak57Pops_22loci_Likelihood_Profile <- LeaveOneOutDistKS.GCL(sillyvec = Kodiak57Pops, loci = loci22, groupvec = Kodiak57PopsGroupVec9)
dput(x = Kodiak57Pops_22loci_Likelihood_Profile, file = "Likelihood Profiles/Kodiak57Pops_22loci_Likelihood_Profile.txt")

Kodiak57Pops_22loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = Kodiak57Pops_22loci_Likelihood_Profile, groupnames = Groups9.nospace.Kodiak, groupvec = Kodiak57PopsGroupVec9, sillyvec = Kodiak57Pops); beep(8)
dput(x = Kodiak57Pops_22loci_Confusion, file = "Objects/Kodiak57Pops_22loci_Confusion.txt")

Kodiak57Pops_22loci_Likelihood_Profile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = Kodiak57Pops, loci = loci22, groupvec = Kodiak57PopsGroupVec9, groupnames = Groups9.nospace.Kodiak, groupcomps = NULL, ncores = 7); beep(5)
dput(x = Kodiak57Pops_22loci_Likelihood_Profile_NEW, file = "Likelihood Profiles/Kodiak57Pops_22loci_Likelihood_Profile_NEW.txt")
PlotLikeProfileKS.GCL(likeprof = Kodiak57Pops_22loci_Likelihood_Profile_NEW, popvec = Kodiak57Pops, loci = loci22, groupvec = Kodiak57PopsGroupVec9, groupnames = Groups9.nospace.Kodiak, dir = "Likelihood Profiles", filename = "Kodiak57Pops_22loci_", col = Colors9)


# loci24 ~~~~~~~~~~
loci24 <- dget(file = "OLD 746 Collections/Objects/loci24.txt")  
Kodiak57Pops_24loci_Likelihood_Profile <- LeaveOneOutDistKS.GCL(sillyvec = Kodiak57Pops, loci = loci24, groupvec = Kodiak57PopsGroupVec9)
dput(x = Kodiak57Pops_24loci_Likelihood_Profile, file = "Likelihood Profiles/Kodiak57Pops_24loci_Likelihood_Profile.txt")

Kodiak57Pops_24loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = Kodiak57Pops_24loci_Likelihood_Profile, groupnames = Groups9.nospace.Kodiak, groupvec = Kodiak57PopsGroupVec9, sillyvec = Kodiak57Pops); beep(8)
dput(x = Kodiak57Pops_24loci_Confusion, file = "Objects/Kodiak57Pops_24loci_Confusion.txt")

Kodiak57Pops_24loci_Likelihood_Profile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = Kodiak57Pops, loci = loci24, groupvec = Kodiak57PopsGroupVec9, groupnames = Groups9.nospace.Kodiak, groupcomps = NULL, ncores = 7); beep(5)
dput(x = Kodiak57Pops_24loci_Likelihood_Profile_NEW, file = "Likelihood Profiles/Kodiak57Pops_24loci_Likelihood_Profile_NEW.txt")
PlotLikeProfileKS.GCL(likeprof = Kodiak57Pops_24loci_Likelihood_Profile_NEW, popvec = Kodiak57Pops, loci = loci24, groupvec = Kodiak57PopsGroupVec9, groupnames = Groups9.nospace.Kodiak, dir = "Likelihood Profiles", filename = "Kodiak57Pops_24loci_", col = Colors9)


# loci24Kodiak ~~~~~~~~~~
loci24Kodiak <- dget(file = "OLD 746 Collections/Objects/loci24Kodiak.txt")  
Kodiak57Pops_24lociKodiak_Likelihood_Profile <- LeaveOneOutDistKS.GCL(sillyvec = Kodiak57Pops, loci = loci24Kodiak, groupvec = Kodiak57PopsGroupVec9)
dput(x = Kodiak57Pops_24lociKodiak_Likelihood_Profile, file = "Likelihood Profiles/Kodiak57Pops_24lociKodiak_Likelihood_Profile.txt")

Kodiak57Pops_24lociKodiak_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = Kodiak57Pops_24lociKodiak_Likelihood_Profile, groupnames = Groups9.nospace.Kodiak, groupvec = Kodiak57PopsGroupVec9, sillyvec = Kodiak57Pops); beep(8)
dput(x = Kodiak57Pops_24lociKodiak_Confusion, file = "Objects/Kodiak57Pops_24lociKodiak_Confusion.txt")


# loci46 ~~~~~~~~~~
loci46 <- dget(file = "OLD 746 Collections/Objects/loci46.txt")
Kodiak57Pops_46loci_Likelihood_Profile <- LeaveOneOutDistKS.GCL(sillyvec = Kodiak57Pops, loci = loci46, groupvec = Kodiak57PopsGroupVec9)
dput(x = Kodiak57Pops_46loci_Likelihood_Profile, file = "Likelihood Profiles/Kodiak57Pops_46loci_Likelihood_Profile.txt")

Kodiak57Pops_46loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = Kodiak57Pops_46loci_Likelihood_Profile, groupnames = Groups9.nospace.Kodiak, groupvec = Kodiak57PopsGroupVec9, sillyvec = Kodiak57Pops); beep(8)
dput(x = Kodiak57Pops_46loci_Confusion, file = "Objects/Kodiak57Pops_46loci_Confusion.txt")


# loci89 ~~~~~~~~~~
loci89 <- dget(file = "Objects/loci89.txt")
Kodiak57Pops_89loci_Likelihood_Profile <- LeaveOneOutDistKS.GCL(sillyvec = Kodiak57Pops, loci = loci89, groupvec = Kodiak57PopsGroupVec9)
dput(x = Kodiak57Pops_89loci_Likelihood_Profile, file = "Likelihood Profiles/Kodiak57Pops_89loci_Likelihood_Profile.txt")

Kodiak57Pops_89loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = Kodiak57Pops_89loci_Likelihood_Profile, groupnames = Groups9.nospace.Kodiak, groupvec = Kodiak57PopsGroupVec9, sillyvec = Kodiak57Pops); beep(8)
dput(x = Kodiak57Pops_89loci_Confusion, file = "Objects/Kodiak57Pops_89loci_Confusion.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Kodiak57Pops_22loci_Likelihood_Profile <- dget(file = "Likelihood Profiles/Kodiak57Pops_22loci_Likelihood_Profile.txt")
Kodiak57Pops_24loci_Likelihood_Profile <- dget(file = "Likelihood Profiles/Kodiak57Pops_24loci_Likelihood_Profile.txt")
Kodiak57Pops_24lociKodiak_Likelihood_Profile <- dget(file = "Likelihood Profiles/Kodiak57Pops_24lociKodiak_Likelihood_Profile.txt")
Kodiak57Pops_46loci_Likelihood_Profile <- dget(file = "Likelihood Profiles/Kodiak57Pops_46loci_Likelihood_Profile.txt")
Kodiak57Pops_89loci_Likelihood_Profile <- dget(file = "Likelihood Profiles/Kodiak57Pops_89loci_Likelihood_Profile.txt")

Kodiak57Pops_22loci_Confusion <- dget(file = "Objects/Kodiak57Pops_22loci_Confusion.txt")
Kodiak57Pops_24loci_Confusion <- dget(file = "Objects/Kodiak57Pops_24loci_Confusion.txt")
Kodiak57Pops_24lociKodiak_Confusion <- dget(file = "Objects/Kodiak57Pops_24lociKodiak_Confusion.txt")
Kodiak57Pops_46loci_Confusion <- dget(file = "Objects/Kodiak57Pops_46loci_Confusion.txt")
Kodiak57Pops_89loci_Confusion <- dget(file = "Objects/Kodiak57Pops_89loci_Confusion.txt")

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(Kodiak57Pops_22loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci22", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))
levelplot(Kodiak57Pops_24loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci24", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))
levelplot(Kodiak57Pops_24lociKodiak_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci24Kodiak", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))
levelplot(Kodiak57Pops_46loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci46", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))
levelplot(Kodiak57Pops_89loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci89", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))

# Individual baseline genetic likelihood for RGs assigned (i.e. what is the genotype likelihood of all baseline individuals to RG X)
suppressWarnings(invisible(lapply(Kodiak57Pops_22loci_Likelihood_Profile[[1]], function(lst) {boxplot(lst, pch=16, notch=TRUE, col=Colors9[Kodiak57PopsGroupVec9], ylab="Probability")} )))  # Regional Flat Prior
suppressWarnings(invisible(lapply(Kodiak57Pops_24loci_Likelihood_Profile[[1]], function(lst) {boxplot(lst, pch=16, notch=TRUE, col=Colors9[Kodiak57PopsGroupVec9], ylab="Probability")} )))  # Regional Flat Prior
suppressWarnings(invisible(lapply(Kodiak57Pops_24lociKodiak_Likelihood_Profile[[1]], function(lst) {boxplot(lst, pch=16, notch=TRUE, col=Colors9[Kodiak57PopsGroupVec9], ylab="Probability")} )))  # Regional Flat Prior
suppressWarnings(invisible(lapply(Kodiak57Pops_46loci_Likelihood_Profile[[1]], function(lst) {boxplot(lst, pch=16, notch=TRUE, col=Colors9[Kodiak57PopsGroupVec9], ylab="Probability")} )))  # Regional Flat Prior
suppressWarnings(invisible(lapply(Kodiak57Pops_89loci_Likelihood_Profile[[1]], function(lst) {boxplot(lst, pch=16, notch=TRUE, col=Colors9[Kodiak57PopsGroupVec9], ylab="Probability")} )))  # Regional Flat Prior

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline")
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")
source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")
source("V:/DATA/R_GEN/JJs GCL/LeaveOneOutlIkeProfile.GCL.R")

KMA473Pops <- dget(file = "Objects/KMA473Pops.txt")
KMA473PopsGroupVec15 <- dget(file = "Objects/KMA473PopsGroupVec15.txt")
Groups15.nospace <- dget(file = "Objects/Groups15.nospace.txt")
Groups15 <- dget(file = "Objects/Groups15.txt")
Colors15 <- dget(file = "Objects/Colors15.txt")

LocusControl <- dget(file = "Objects/LocusControl103.txt")
loci103 <- dget(file = "Objects/loci103.txt")
invisible(sapply(KMA473Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPopsloci103/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)

# loci46 ~~~~~~~~~~
loci46 <- dget(file = "OLD 746 Collections/Objects/loci46.txt")

KMA473Pops_46loci_Likelihood_Profile <- LeaveOneOutDistKS.GCL(sillyvec = KMA473Pops, loci = loci46, groupvec = KMA473PopsGroupVec15)
dput(x = KMA473Pops_46loci_Likelihood_Profile, file = "Likelihood Profiles/KMA473Pops_46loci_Likelihood_Profile.txt")

KMA473Pops_46loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = KMA473Pops_46loci_Likelihood_Profile, groupnames = Groups15, groupvec = KMA473PopsGroupVec15, sillyvec = KMA473Pops); beep(8)
dput(x = KMA473Pops_46loci_Confusion, file = "Objects/KMA473Pops_46loci_Confusion.txt")


KMA473Pops_46loci_Confusion <- dget(file = "Objects/KMA473Pops_46loci_Confusion.txt")
require(lattice)
require(devEMF)

new.colors <- colorRampPalette(c("white", "black"))

emf(file = "Likelihood Profiles/KMA473Pops_46loci_Confusion.emf", width = 6.5, height = 6.5, family = "Times")
levelplot(KMA473Pops_46loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", 
          at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 90)))
dev.off()

KMA473Pops_46loci_Likelihood_Profile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = KMA473Pops, loci = loci46, groupvec = KMA473PopsGroupVec15, groupnames = Groups15.nospace, groupcomps = NULL, ncores = 7); beep(5)
dput(x = KMA473Pops_46loci_Likelihood_Profile_NEW, file = "Likelihood Profiles/KMA473Pops_46loci_Likelihood_Profile_NEW.txt")
PlotLikeProfileKS.GCL(likeprof = KMA473Pops_46loci_Likelihood_Profile_NEW, popvec = KMA473Pops, loci = loci46, groupvec = KMA473PopsGroupVec15, groupnames = Groups15.nospace, dir = "Likelihood Profiles", filename = "KMA473Pops_46loci_", col = Colors15)
likeprof = KMA473Pops_46loci_Likelihood_Profile_NEW; popvec = KMA473Pops; loci = loci46; groupvec = KMA473PopsGroupVec15; groupnames = Groups15.nospace; dir = "Likelihood Profiles"; filename = "KMA473Pops_46loci_"; col = Colors15

# loci24Kodiak ~~~~~~~~~~
loci24Kodiak <- dget(file = "OLD 746 Collections/Objects/loci24Kodiak.txt")

KMA473Pops_24lociKodiak_Likelihood_Profile <- LeaveOneOutDistKS.GCL(sillyvec = KMA473Pops, loci = loci24Kodiak, groupvec = KMA473PopsGroupVec15)
dput(x = KMA473Pops_24lociKodiak_Likelihood_Profile, file = "Likelihood Profiles/KMA473Pops_24lociKodiak_Likelihood_Profile.txt")

KMA473Pops_24lociKodiak_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = KMA473Pops_24lociKodiak_Likelihood_Profile, groupnames = Groups15, groupvec = KMA473PopsGroupVec15, sillyvec = KMA473Pops); beep(8)
dput(x = KMA473Pops_24lociKodiak_Confusion, file = "Objects/KMA473Pops_24lociKodiak_Confusion.txt")

KMA473Pops_24lociKodiak_Likelihood_Profile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = KMA473Pops, loci = loci24Kodiak, groupvec = KMA473PopsGroupVec15, groupnames = Groups15.nospace, groupcomps = NULL, ncores = 7); beep(5)
dput(x = KMA473Pops_24lociKodiak_Likelihood_Profile_NEW, file = "Likelihood Profiles/KMA473Pops_24lociKodiak_Likelihood_Profile_NEW.txt")
PlotLikeProfileKS.GCL(likeprof = KMA473Pops_24lociKodiak_Likelihood_Profile_NEW, popvec = KMA473Pops, loci = loci24Kodiak, groupvec = KMA473PopsGroupVec15, groupnames = Groups15.nospace, dir = "Likelihood Profiles", filename = "KMA473Pops_24lociKodiak_", col = Colors15)

# loci24 ~~~~~~~~~~
loci24 <- dget(file = "OLD 746 Collections/Objects/loci24.txt")

KMA473Pops_24loci_Likelihood_Profile <- LeaveOneOutDistKS.GCL(sillyvec = KMA473Pops, loci = loci24, groupvec = KMA473PopsGroupVec15)
dput(x = KMA473Pops_24loci_Likelihood_Profile, file = "Likelihood Profiles/KMA473Pops_24loci_Likelihood_Profile.txt")

KMA473Pops_24loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = KMA473Pops_24loci_Likelihood_Profile, groupnames = Groups15, groupvec = KMA473PopsGroupVec15, sillyvec = KMA473Pops); beep(8)
dput(x = KMA473Pops_24loci_Confusion, file = "Objects/KMA473Pops_24loci_Confusion.txt")

KMA473Pops_24loci_Likelihood_Profile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = KMA473Pops, loci = loci24, groupvec = KMA473PopsGroupVec15, groupnames = Groups15.nospace, groupcomps = NULL, ncores = 7); beep(5)
dput(x = KMA473Pops_24loci_Likelihood_Profile_NEW, file = "Likelihood Profiles/KMA473Pops_24loci_Likelihood_Profile_NEW.txt")
PlotLikeProfileKS.GCL(likeprof = KMA473Pops_24loci_Likelihood_Profile_NEW, popvec = KMA473Pops, loci = loci24, groupvec = KMA473PopsGroupVec15, groupnames = Groups15.nospace, dir = "Likelihood Profiles", filename = "KMA473Pops_24loci_", col = Colors15)

# loci22 ~~~~~~~~~~
loci22 <- dget(file = "OLD 746 Collections/Objects/loci22.txt")

KMA473Pops_22loci_Likelihood_Profile <- LeaveOneOutDistKS.GCL(sillyvec = KMA473Pops, loci = loci22, groupvec = KMA473PopsGroupVec15)
dput(x = KMA473Pops_22loci_Likelihood_Profile, file = "Likelihood Profiles/KMA473Pops_22loci_Likelihood_Profile.txt")

KMA473Pops_22loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = KMA473Pops_22loci_Likelihood_Profile, groupnames = Groups15, groupvec = KMA473PopsGroupVec15, sillyvec = KMA473Pops); beep(8)
dput(x = KMA473Pops_22loci_Confusion, file = "Objects/KMA473Pops_22loci_Confusion.txt")

KMA473Pops_22loci_Likelihood_Profile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = KMA473Pops, loci = loci22, groupvec = KMA473PopsGroupVec15, groupnames = Groups15.nospace, groupcomps = NULL, ncores = 7); beep(5)
dput(x = KMA473Pops_22loci_Likelihood_Profile_NEW, file = "Likelihood Profiles/KMA473Pops_22loci_Likelihood_Profile_NEW.txt")
PlotLikeProfileKS.GCL(likeprof = KMA473Pops_22loci_Likelihood_Profile_NEW, popvec = KMA473Pops, loci = loci22, groupvec = KMA473PopsGroupVec15, groupnames = Groups15.nospace, dir = "Likelihood Profiles", filename = "KMA473Pops_22loci_", col = Colors15)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(lattice)
new.colors <- colorRampPalette(c("white", "black"))

KMA473Pops_24loci_Confusion <- dget(file = "Objects/KMA473Pops_24loci_Confusion.txt")
KMA473Pops_24lociKodiak_Confusion <- dget(file = "Objects/KMA473Pops_24lociKodiak_Confusion.txt")
KMA473Pops_46loci_Confusion <- dget(file = "Objects/KMA473Pops_46loci_Confusion.txt")
KMA473Pops_89loci_Confusion <- dget(file = "Objects/KMA473Pops_89loci_Confusion.txt")

levelplot(KMA473Pops_24loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci24", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))
levelplot(KMA473Pops_24lociKodiak_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci24Kodiak", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))
levelplot(KMA473Pops_46loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci46", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))
levelplot(KMA473Pops_89loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci89", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))

range(KMA473Pops_46loci_Confusion[[1]] - KMA473Pops_89loci_Confusion[[1]])
levelplot(KMA473Pops_46loci_Confusion[[1]] - KMA473Pops_89loci_Confusion[[1]], col.regions = colorRampPalette(c("red", "white", "blue")), xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "loci89", at = seq(-0.2, 0.2, length.out = 100), scales = list(x = list(rot = 45)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### MDS 24 Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Want to confirm that the MDS look good before I spend a lot of time on repeated proof tests
## Dump a GENEPOP file for genind
gcl2Genepop.GCL(sillyvec = KMA473Pops, loci = loci22, path = "Genepop/KMA473Pops_22loci.gen", VialNums = TRUE)

require(adegenet)
genind <- read.genepop(file = "Genepop/KMA473Pops_22loci.gen")

genpop <- genind2genpop(genind)

AdegenetNei473Pop22loci <- dist.genpop(genpop, method = 1, diag = TRUE, upper = TRUE)
dput(x = AdegenetNei473Pop22loci, file = "Trees/AdegenetNei473Pop22loci.txt")
str(AdegenetNei473Pop22loci)

require(ape)
Nei473NJtree <- nj(AdegenetNei473Pop22loci)
str(Nei473NJtree)
plot.phylo(x = Nei473NJtree, cex = 0.5, no.margin = TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MDS
library('rgl')

# MDS <- cmdscale(as.matrix(AdegenetNei473Pop22loci), k = 3)  # Did in base R as it kept crashing in RStudio...dunno why.
# MDS <- cmdscale(as.matrix(AdegenetNei473Pop22loci), k = 463, eig = TRUE)  # Do with all possible dimensions
# dput(x = MDS, file = "Objects/MDSAdegenetNei473Pop22loci.txt")
# dput(x = MDS, file = "Objects/MDSAdegenetNei473Pop22loci_alldim.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei473Pop22loci.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei473Pop93nuclearloci.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei473Pop22loci_alldim.txt")

x <- as.vector(MDS$points[, 1])   
y <- as.vector(MDS$points[, 2])
z <- as.vector(MDS$points[, 3])

KMA473PopsGroupVec15 <- dget(file = "Objects/KMA473PopsGroupVec15.txt")
Colors15 <- dget(file = "Objects/Colors15.txt")


open3d()
par(family = "serif")
plot3d(x, y, z + abs(range(z)[1]), xlab = '', ylab = '', zlab = '', aspect = FALSE, col = Colors15[KMA473PopsGroupVec15], size = 0.5, type = 's', axes = TRUE, box = TRUE, top = TRUE, cex = 1)
box3d()
# plot3d(x, y, z + abs(range(z)[1]), aspect = FALSE, col = "black", size = 0.5, type = 'h', box = TRUE, axes = FALSE, top = FALSE, add = TRUE)  # adds pins to spheres
texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.2, 0), text = seq(KMA473Pops), font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)

# Plot individually
par(mar = c(4.1, 4.1, 0.6, 0.6))
plot(x, col = Colors15[KMA473PopsGroupVec15], pch = 16, xlab = "Pop", bty = "n")
CommonNames473[which.min(x)]; CommonNames473[which.max(x)]
plot(y, col = Colors15[KMA473PopsGroupVec15], pch = 16, xlab = "Pop", bty = "n")
CommonNames473[which.min(y)]; CommonNames473[which.max(y)]
plot(z, col = Colors15[KMA473PopsGroupVec15], pch = 16, xlab = "Pop", bty = "n")
CommonNames473[which.min(z)]; CommonNames473[which.max(z)]

# Look at MDS
plot(MDS$eig, type = "l", xlab = "Dimension", ylab = "Eigen Value", bty = "l", lwd = 3)
plot(MDS$eig[1:20], type = "l", xlab = "Dimension", ylab = "Eigen Value", bty = "l", lwd = 3, ylim = c(0, 2))
hist(MDS$eig[MDS$eig > 0], xlab = "Eigenvalue", breaks = 20, col = 8)

sapply(1:20, function(i) {plot(MDS$points[, i], col = Colors15[KMA473PopsGroupVec15], pch = 16, xlab = "Pop", bty = "n", ylab = paste("Dimensions", i), main = round(MDS$eig[i], 3), ylim = c(-0.2, 0.4))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Karluk

Karluk100ProofTests <- grep(pattern = "Karluk", x = KMAGroups15Repeated100ProofTests, value = TRUE)
Karluk100ProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, maindir = "BAYES/BAYES.output", 
                                                             mixvec = Karluk100ProofTests, prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1,
                                                             PosteriorOutput = TRUE); beep(8)

dput(x = Karluk100ProofTestsEstimates, file = "Estimates objects/Karluk100ProofTestsEstimates.txt")
Karluk100ProofTestsEstimates <- dget(file = "Estimates objects/Karluk100ProofTestsEstimates.txt")

str(Karluk100ProofTestsEstimates)
Karluk100ProofTestsEstimates$Stats

## Check Gelman-Rubin factor
lapply(Karluk100ProofTestsEstimates$Stats, function(Proof) {Proof[, "GR"]} )
lapply(Karluk100ProofTestsEstimates$Stats, function(Proof) {table(Proof[, "GR"] > 1.2)} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot results
# Basic Bars
#~~~~~~~~~~~~~~~~~~
library(gplots)
library(plotrix)

str(Karluk100ProofTestsEstimates)

Groups15.nospace

samplesz_proof_final

## Karluk
par(family = "serif")
ProofPlot <- barplot2(height = t(matrix(data = c(Karluk100ProofTestsEstimates$Stats$Karluk1[, 3], Karluk100ProofTestsEstimates$Stats$Karluk2[, 3],
                                                 Karluk100ProofTestsEstimates$Stats$Karluk3[, 3], Karluk100ProofTestsEstimates$Stats$Karluk4[, 3],
                                                 Karluk100ProofTestsEstimates$Stats$Karluk5[, 3]), nrow = 15, ncol = 5, byrow = FALSE)),
                      ci.l = t(matrix(data = c(Karluk100ProofTestsEstimates$Stats$Karluk1[, 4], Karluk100ProofTestsEstimates$Stats$Karluk2[, 4],
                                               Karluk100ProofTestsEstimates$Stats$Karluk3[, 4], Karluk100ProofTestsEstimates$Stats$Karluk4[, 4],
                                               Karluk100ProofTestsEstimates$Stats$Karluk5[, 4]), nrow = 15, ncol = 5, byrow = FALSE)),
                      ci.u = t(matrix(data = c(Karluk100ProofTestsEstimates$Stats$Karluk1[, 5], Karluk100ProofTestsEstimates$Stats$Karluk2[, 5],
                                               Karluk100ProofTestsEstimates$Stats$Karluk3[, 5], Karluk100ProofTestsEstimates$Stats$Karluk4[, 5],
                                               Karluk100ProofTestsEstimates$Stats$Karluk5[, 5]), nrow = 15, ncol = 5, byrow = FALSE)), names.arg = NULL, 
                      col = 'blue', main = paste("Karluk (n = ", samplesz_proof_final["Karluk"], ")", sep = ''), xlab = "Reporting Group", ylab = "Proportion",
                      plot.ci = TRUE, beside = TRUE, ylim = c(0, 1))
abline(h = 0)
abline(h = 0.9, lwd = 2)
text(x = ProofPlot[3, ], y = -0.065, labels = Groups15.nospace, srt = 45, xpd = NA, cex = 0.75)

#~~~~~~~~~~~~~~~~~~
# Violinplot

require(vioplot)

substrRight <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

Karluk100ProofTests

par(mar = c(5.1, 4.6, 3.6, 6.6))

plot(y = t(sapply(Karluk100ProofTestsEstimates$Stats, function(Proof) {Proof[, "mean"]} )), x = seq(5 * length(Groups15.nospace)), pch =16, col = rep(Colors15, each = 5), ylab = "", ylim = c(0, 1), xlab = "", axes = FALSE, main = "", cex = 2)
invisible(sapply(Karluk100ProofTests, function(Proof) {sapply(1:15, function(i) {vioplot2(Karluk100ProofTestsEstimates$Output[[Proof]][, i], at = as.numeric(substrRight(x = Proof, n = 1)) + ((i - 1) * 5), horizontal = FALSE, col = Colors15[i], border = TRUE, drawRect = FALSE, rectCol = Colors15[i], add = TRUE, wex = 1, lwd = 1)} )} ))
arrows(x0 = seq(5 * length(Groups15.nospace)), x1 = seq(5 * length(Groups15.nospace)), y0 = t(sapply(Karluk100ProofTestsEstimates$Stats, function(Proof) {Proof[, "5%"]} )), y1 = t(sapply(Karluk100ProofTestsEstimates$Stats, function(Proof) {Proof[, "95%"]} )), angle = 90, code = 3, length = 0.1, lwd = 2)
points(y = t(sapply(Karluk100ProofTestsEstimates$Stats, function(Proof) {Proof[, "mean"]} )), x = seq(5 * length(Groups15.nospace)), cex = 1.5, pch = 21, col = "white", bg = rep(Colors15, each = 5), lwd = 2)
axis(side = 2, lwd = 3, cex.axis = 1.5, pos = 0)
mtext(text = "Proportion", side = 2, line = 1, cex = 2)
mtext(text = paste("Karluk (n = ", samplesz_proof_final["Karluk"], ")", sep = ''), side = 3, cex = 2.5)
text(x = seq(from = 3, by = 5, length.out = 15) - 0.35, y = 0, labels = Groups15Short, srt = 45, pos = 1, offset = 2.5, xpd = TRUE, cex = 1.6)
segments(x0 = 0, y0 = 0, y1 = 0, x1 = 5 * length(Groups15.nospace), lwd = 3)
segments(x0 = 0, y0 = 0.9, y1 = 0.9, x1 = 5 * length(Groups15.nospace), lwd = 3)
#legend(x = 70, y = 0.9, xpd = TRUE, legend = Groups15Short, fill = Colors15, bty = "n")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Re-do Proof Tests with Different Markerset ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Feed the proof test objects with fish IDs to re-proof test to run the same fish for different markerset
# Switching from 5 Repeats with 5 Chains to only 1 Chain in order to save time.
# Only difference between chains is the inits.
# How to pick inits for only 1 Chain? Random? Or just use the first chain inits from the 5 chains?

# code from Repeated proof test
sillyvec = KMA473Pops.named
nchains = 1
groupnames = Groups15.nospace
G = length(groupnames)
groupvec = KMA473PopsGroupVec15


initmat=array(NA,c(length(sillyvec),nchains),dimnames=list(sillyvec,paste("Chain",1:nchains,sep="")))
for(chain in 1:nchains){
  for(g in 1:G){
    IND=groupvec==g
    gg=rgamma(sum(IND),1/G,1) 
    initmat[sillyvec[IND],chain]=gg 
  }
  initmat[,chain]=initmat[,chain]/sum(initmat[,chain])   
}
plot(initmat[, 1], type = "h")

plot(round(sapply(seq(G), function(g) {sum(initmat[groupvec == g, 1])} ), 3), type = "h")

plot(table(groupvec) / length(sillyvec), type = "h")

sum(initmat[groupvec == 1, 1])
sum(initmat[, 1])

# Rather than do a random init, just use the first part of the multichain inits so that we can compare apples to apples
# This means that the first 1/5th populations carry 90%

# Confirm that chains are similar in their 90% CI
# Note that this is due to both initial starting values AND the random seeds
KMA473PopsGroups15Repeated100ProofTestsEstimatesFrazer <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, maindir = "BAYES/100% Proof Tests/loci89/BAYES.output", 
                                                                                       mixvec = KMA473PopsGroups15Repeated100ProofTests[21:25], prior = "", ext = "RGN", nchains = 5, burn = 0.5, alpha = 0.1,
                                                                                       PosteriorOutput = TRUE); beep(8)

myquant <- function(x) {
  quantile(x = x, probs = c(0.05, 0.5, 0.95))
}

lapply(KMA473PopsGroups15Repeated100ProofTestsEstimatesFrazer$Output, function(peat) {
  rbind(myquant(peat[1:20000, 5]), myquant(peat[20001:40000, 5]), myquant(peat[40001:60000, 5]), myquant(peat[60001:80000, 5]), myquant(peat[80001:100000, 5]), myquant(peat[, 5]))
})

# Close enough, just use "KMA473PopsInits[, 1, drop = FALSE]"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### loci46 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

loci46 <- dget(file = "OLD 746 Collections/objects/loci46.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 100% Proof Tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Just doing Frazer and Ayakulik, as these were the "worst" performing 100% proof tests
dir.create("BAYES/100% Proof Tests/loci46")

KMA473PopsGroups15Repeated100ProofTests

#~~~~~~~~~~~~~~~~~~
## Get Proof Objects
# Frazer
for(proof in grep(pattern = "Frazer", x = KMA473PopsGroups15Repeated100ProofTests, value = TRUE)) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

# Ayakulik
for(proof in grep(pattern = "Ayakulik", x = KMA473PopsGroups15Repeated100ProofTests, value = TRUE)) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

#~~~~~~~~~~~~~~~~~~
## Create BAYES files
for (proof in paste(rep(x = c("Ayakulik", "Frazer"), each = 5), 1:5, sep = "")) {
  ReProofTestKS.GCL(sillyvec = KMA473Pops.named, loci = loci46, groupnames = Groups15.nospace, groupvec = KMA473PopsGroupVec15, 
                    ProofTestIDs = get(paste(proof, "Proof", sep = "")), prefix = proof, dir = "BAYES/100% Proof Tests/loci46", 
                    suffix = "", nreps = 40000, nchains = 1, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits[, 1, drop = FALSE], 
                    type = "BAYES", thin = c(1, 1, 100), switches = "F T F T T T F")
}

# Create BAYES Output folders
dir.create("BAYES/100% Proof Tests/loci46/BAYES.output")
sapply(paste(rep(x = c("Frazer", "Ayakulik"), each = 5), 1:5, sep = ""), function(proof) {dir <- paste(getwd(), "/BAYES/100% Proof Tests/loci46/BAYES.output/", proof, sep = "")
                                                                                          dir.create(dir) })

#~~~~~~~~~~~~~~~~~~
## Summarize Results
FrazerAyakulikloci46Repeated100ProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, maindir = "BAYES/100% Proof Tests/loci46/BAYES.output", 
                                                                                   mixvec = paste(rep(x = c("Frazer", "Ayakulik"), each = 5), 1:5, sep = ""), prior = "", ext = "RGN", nchains = 1, burn = 0.5, 
                                                                                   alpha = 0.1, PosteriorOutput = TRUE); beep(5)

str(FrazerAyakulikloci46Repeated100ProofTestsEstimates)
str(list(Stats = FrazerAyakulikloci46Repeated100ProofTestsEstimates$Stats[1:5], Output = FrazerAyakulikloci46Repeated100ProofTestsEstimates$Output[1:5]))

dir.create("Estimates objects/loci46")
dput(x = list(Stats = FrazerAyakulikloci46Repeated100ProofTestsEstimates$Stats[1:5], Output = FrazerAyakulikloci46Repeated100ProofTestsEstimates$Output[1:5]), 
     file = "Estimates objects/loci46/Frazer100ProofTestEstimates.txt")
dput(x = list(Stats = FrazerAyakulikloci46Repeated100ProofTestsEstimates$Stats[6:10], Output = FrazerAyakulikloci46Repeated100ProofTestsEstimates$Output[6:10]), 
     file = "Estimates objects/loci46/Ayakulik100ProofTestEstimates.txt")


#~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~

library(gplots)
library(plotrix)

Groups15.nospace

ProofTest100.SampleSize

ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(5.1, 4.4, 2.1, 1.1))

sapply(c("Frazer", "Ayakulik"), function(RG) {
  par(xpd = FALSE)
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {FrazerAyakulikloci46Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {FrazerAyakulikloci46Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {FrazerAyakulikloci46Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", ProofTest100.SampleSize[RG], ")", sep = ''), xlab = "",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 1), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  abline(h = 0.9, lwd = 2)
  text(x = ProofPlot[3, ], y = -0.065, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.8)
  mtext(text = "Reporting Group", side = 1, line = 3.5)
} )

dir.create("BAYES/100% Proof Tests/loci46/Figures")


#~~~~~~~~~~~~~~~~~~
# Expanding to all RGs
#~~~~~~~~~~~~~~~~~~
## Get Proof Objects
# All non-Frazer/Ayakulik
for(proof in paste(rep(x = Groups15.nospace[-c(5:6)], each = 5), 1:5, sep = "")) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

#~~~~~~~~~~~~~~~~~~
## Create BAYES files
for (proof in paste(rep(x = Groups15.nospace[-c(5:6)], each = 5), 1:5, sep = "")) {
  ReProofTestKS.GCL(sillyvec = KMA473Pops.named, loci = loci46, groupnames = Groups15.nospace, groupvec = KMA473PopsGroupVec15, 
                    ProofTestIDs = get(paste(proof, "Proof", sep = "")), prefix = proof, dir = "BAYES/100% Proof Tests/loci46", 
                    suffix = "", nreps = 40000, nchains = 1, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits[, 1, drop = FALSE], 
                    type = "BAYES", thin = c(1, 1, 100), switches = "F T F T T T F")
}


# Create BAYES Output folders
sapply(paste(rep(x = Groups15.nospace[-c(5:6)], each = 5), 1:5, sep = ""), function(proof) {dir <- paste(getwd(), "/BAYES/100% Proof Tests/loci46/BAYES.output/", proof, sep = "")
                                                                                            dir.create(dir) })

#~~~~~~~~~~~~~~~~~~
## Summarize Results
KMA473PopsGroups15loci46Repeated100ProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, maindir = "BAYES/100% Proof Tests/loci46/BAYES.output", 
                                                                                       mixvec = KMA473PopsGroups15Repeated100ProofTests, prior = "", ext = "RGN", nchains = 1, burn = 0.5, 
                                                                                       alpha = 0.1, PosteriorOutput = TRUE); beep(5)
str(KMA473PopsGroups15loci46Repeated100ProofTestsEstimates)

dput(x = KMA473PopsGroups15loci46Repeated100ProofTestsEstimates, file = "Estimates objects/loci46/KMA473PopsGroups15loci46Repeated100ProofTestsEstimates.txt")
dput(x = KMA473PopsGroups15loci46Repeated100ProofTestsEstimates$Stats, file = "Estimates objects/loci46/KMA473PopsGroups15loci46Repeated100ProofTestsEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~

library(gplots)
library(plotrix)

Groups15.nospace

ProofTest100.SampleSize

ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(5.1, 4.4, 2.1, 1.1))

sapply(Groups15.nospace, function(RG) {
  par(xpd = FALSE)
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {KMA473PopsGroups15loci46Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {KMA473PopsGroups15loci46Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {KMA473PopsGroups15loci46Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", ProofTest100.SampleSize[RG], ")", sep = ''), xlab = "",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 1), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  abline(h = 0.9, lwd = 2)
  text(x = ProofPlot[3, ], y = -0.065, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.8)
  mtext(text = "Reporting Group", side = 1, line = 3.5)
} )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Fishery Scenario Tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create("BAYES/Mixture Proof Tests/loci46")

KMA473PopsGroups15RepeatedMixProofTests

#~~~~~~~~~~~~~~~~~~
## Get Proof Objects
for(proof in KMA473PopsGroups15RepeatedMixProofTests) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("MixtureProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

#~~~~~~~~~~~~~~~~~~
## Create BAYES files
for (proof in KMA473PopsGroups15RepeatedMixProofTests){
  ReProofTestKS.GCL(sillyvec = KMA473Pops, loci = loci46, groupnames = Groups15Short, groupvec = KMA473PopsGroupVec15, 
                    ProofTestIDs = get(paste(proof, "Proof", sep = "")), prefix = proof, dir = "BAYES/Mixture Proof Tests/loci46", 
                    suffix = "", nreps = 40000, nchains = 1, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits[, 1, drop = FALSE], 
                    type = "BAYES", thin = c(1, 1, 100), switches = "F T F T T T F")
}

# Create BAYES Output folders
dir.create("BAYES/Mixture Proof Tests/loci46/BAYES.output")
invisible(sapply(KMA473PopsGroups15RepeatedMixProofTests, function(proof) {dir <- paste(getwd(), "/BAYES/Mixture Proof Tests/loci46/BAYES.output/", proof, sep = "")
                                                                           dir.create(dir) } ))

#~~~~~~~~~~~~~~~~~~
## Summarize Results
KMA473PopsGroups15loci46RepeatedMixProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, 
                                                                                       maindir = "BAYES/Mixture Proof Tests/loci46/BAYES.output", 
                                                                                       mixvec = KMA473PopsGroups15RepeatedMixProofTests, prior = "", 
                                                                                       ext = "RGN", nchains = 1, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE); beep(5)
str(KMA473PopsGroups15loci46RepeatedMixProofTestsEstimates)

dput(x = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimates, file = "Estimates objects/loci46/KMA473PopsGroups15loci46RepeatedMixProofTestsEstimates.txt")
dput(x = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimates$Stats, file = "Estimates objects/loci46/KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~

library(gplots)
library(plotrix)

KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci46/KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats.txt")

str(KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats)

Groups15.nospace

MixtureProofTest.SampleSize.Final <- round(apply(MixtureProofTest.SampleSize, 1, sum), 0 )

MixtureProofTestProportions <- read.table(file = "FisheryProofTestScenarios.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = Groups15Short)
MixtureProofTestProportions <- MixtureProofTestProportions[, -1]
max(MixtureProofTestProportions)

# PCGroups15 <- c("West of Chignik", "Black Lake", "Chignik Lake", "Upper Station / Akalura", "Frazer", "Ayakulik", "Karluk", "Uganik", "Northwest Kodiak", "Afognak", "Eastside Kodiak", "Saltery", "Cook Inlet", "Prince William Sound", "South of Cape Suckling")
# dput(x = PCGroups15, file = "Objects/PCGroups15.txt")
PCGroups15 <- dget(file = "Objects/PCGroups15.txt")
ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(3.1, 4.4, 2.1, 10.1))
par(xpd = FALSE)
FisheryProofTestScenarioNames <- dget(file = "Objects/FisheryProofTestScenarioNames.txt")
# RG = "flat"
for(RG in FisheryProofTestScenarioNames) {
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", MixtureProofTest.SampleSize.Final[RG], ")", sep = ''), xlab = "Reporting Group",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 0.6), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  #abline(h = 1/15, lwd = 4)
  
  for(i in 1:15) {
    segments(x0 = ProofPlot[1, i] - 0.5, x1 = ProofPlot[5, i] + 0.5, y0 = MixtureProofTestProportions[i, RG], y1 = MixtureProofTestProportions[i, RG], lwd = 3)
  }
  
  legend(x = max(ProofPlot), y = 0.5, legend = PCGroups15, fill = Colors15, bty = "n", xpd = TRUE)
  
  text(x = ProofPlot[3, ], y = -0.03, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.7)
}

dir.create("BAYES/Mixture Proof Tests/loci46/Figures")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### loci24Kodiak ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

loci24Kodiak <- dget(file = "OLD 746 Collections/objects/loci24Kodiak.txt")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 100% Proof Tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Just doing Frazer and Ayakulik, as these were the "worst" performing 100% proof tests
dir.create("BAYES/100% Proof Tests/loci24Kodiak")

KMA473PopsGroups15Repeated100ProofTests

#~~~~~~~~~~~~~~~~~~
## Get Proof Objects
# Frazer
for(proof in grep(pattern = "Frazer", x = KMA473PopsGroups15Repeated100ProofTests, value = TRUE)) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

# Ayakulik
for(proof in grep(pattern = "Ayakulik", x = KMA473PopsGroups15Repeated100ProofTests, value = TRUE)) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

#~~~~~~~~~~~~~~~~~~
## Create BAYES files
for (proof in paste(rep(x = c("Ayakulik", "Frazer"), each = 5), 1:5, sep = "")) {
  ReProofTestKS.GCL(sillyvec = KMA473Pops.named, loci = loci24Kodiak, groupnames = Groups15.nospace, groupvec = KMA473PopsGroupVec15, 
                    ProofTestIDs = get(paste(proof, "Proof", sep = "")), prefix = proof, dir = "BAYES/100% Proof Tests/loci24Kodiak", 
                    suffix = "", nreps = 40000, nchains = 1, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits[, 1, drop = FALSE], 
                    type = "BAYES", thin = c(1, 1, 100), switches = "F T F T T T F")
}


# Create BAYES Output folders
dir.create("BAYES/100% Proof Tests/loci24Kodiak/BAYES.output")
sapply(paste(rep(x = c("Frazer", "Ayakulik"), each = 5), 1:5, sep = ""), function(proof) {dir <- paste(getwd(), "/BAYES/100% Proof Tests/loci24Kodiak/BAYES.output/", proof, sep = "")
                                                                                          dir.create(dir) })

#~~~~~~~~~~~~~~~~~~
## Summarize Results
FrazerAyakulikloci24KodiakRepeated100ProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, maindir = "BAYES/100% Proof Tests/loci24Kodiak/BAYES.output", 
                                                                                   mixvec = paste(rep(x = c("Frazer", "Ayakulik"), each = 5), 1:5, sep = ""), prior = "", ext = "RGN", nchains = 1, burn = 0.5, 
                                                                                   alpha = 0.1, PosteriorOutput = TRUE); beep(5)
str(FrazerAyakulikloci24KodiakRepeated100ProofTestsEstimates)

dir.create("Estimates objects/loci24Kodiak")
dput(x = list(Stats = FrazerAyakulikloci24KodiakRepeated100ProofTestsEstimates$Stats[1:5], Output = FrazerAyakulikloci24KodiakRepeated100ProofTestsEstimates$Output[1:5]), 
     file = "Estimates objects/loci24Kodiak/Frazer100ProofTestEstimates.txt")
dput(x = list(Stats = FrazerAyakulikloci24KodiakRepeated100ProofTestsEstimates$Stats[6:10], Output = FrazerAyakulikloci24KodiakRepeated100ProofTestsEstimates$Output[6:10]), 
     file = "Estimates objects/loci24Kodiak/Ayakulik100ProofTestEstimates.txt")

#~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~

library(gplots)
library(plotrix)

Groups15.nospace

ProofTest100.SampleSize

ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(5.1, 4.4, 2.1, 1.1))

sapply(c("Frazer", "Ayakulik"), function(RG) {
  par(xpd = FALSE)
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {FrazerAyakulikloci24KodiakRepeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {FrazerAyakulikloci24KodiakRepeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {FrazerAyakulikloci24KodiakRepeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", ProofTest100.SampleSize[RG], ")", sep = ''), xlab = "",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 1), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  abline(h = 0.9, lwd = 2)
  text(x = ProofPlot[3, ], y = -0.065, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.8)
  mtext(text = "Reporting Group", side = 1, line = 3.5)
} )

dir.create("BAYES/100% Proof Tests/loci24Kodiak/Figures")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Fishery Scenario Tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create("BAYES/Mixture Proof Tests/loci24Kodiak")

KMA473PopsGroups15RepeatedMixProofTests

#~~~~~~~~~~~~~~~~~~
## Get Proof Objects
for(proof in KMA473PopsGroups15RepeatedMixProofTests) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("MixtureProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

#~~~~~~~~~~~~~~~~~~
## Create BAYES files
for (proof in KMA473PopsGroups15RepeatedMixProofTests){
  ReProofTestKS.GCL(sillyvec = KMA473Pops, loci = loci24Kodiak, groupnames = Groups15Short, groupvec = KMA473PopsGroupVec15, 
                    ProofTestIDs = get(paste(proof, "Proof", sep = "")), prefix = proof, dir = "BAYES/Mixture Proof Tests/loci24Kodiak", 
                    suffix = "", nreps = 40000, nchains = 1, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits[, 1, drop = FALSE], 
                    type = "BAYES", thin = c(1, 1, 100), switches = "F T F T T T F")
}

# Create BAYES Output folders
dir.create("BAYES/Mixture Proof Tests/loci24Kodiak/BAYES.output")
invisible(sapply(KMA473PopsGroups15RepeatedMixProofTests, function(proof) {dir <- paste(getwd(), "/BAYES/Mixture Proof Tests/loci24Kodiak/BAYES.output/", proof, sep = "")
                                                                           dir.create(dir) } ))


#~~~~~~~~~~~~~~~~~~
## Summarize Results
KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, 
                                                                                       maindir = "BAYES/Mixture Proof Tests/loci24Kodiak/BAYES.output", 
                                                                                       mixvec = KMA473PopsGroups15RepeatedMixProofTests, prior = "", 
                                                                                       ext = "RGN", nchains = 1, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE); beep(5)
str(KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimates)

dput(x = KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimates, file = "Estimates objects/loci24Kodiak/KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimates.txt")
dput(x = KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimates$Stats, file = "Estimates objects/loci24Kodiak/KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~

library(gplots)
library(plotrix)

KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci24Kodiak/KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats.txt")

str(KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats)

Groups15.nospace

MixtureProofTest.SampleSize.Final <- round(apply(MixtureProofTest.SampleSize, 1, sum), 0 )

MixtureProofTestProportions <- read.table(file = "FisheryProofTestScenarios.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = Groups15Short)
MixtureProofTestProportions <- MixtureProofTestProportions[, -1]
max(MixtureProofTestProportions)

# PCGroups15 <- c("West of Chignik", "Black Lake", "Chignik Lake", "Upper Station / Akalura", "Frazer", "Ayakulik", "Karluk", "Uganik", "Northwest Kodiak", "Afognak", "Eastside Kodiak", "Saltery", "Cook Inlet", "Prince William Sound", "South of Cape Suckling")
# dput(x = PCGroups15, file = "Objects/PCGroups15.txt")
PCGroups15 <- dget(file = "Objects/PCGroups15.txt")
ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(3.1, 4.4, 2.1, 10.1))
par(xpd = FALSE)
FisheryProofTestScenarioNames <- dget(file = "Objects/FisheryProofTestScenarioNames.txt")
# RG = "flat"
for(RG in FisheryProofTestScenarioNames) {
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", MixtureProofTest.SampleSize.Final[RG], ")", sep = ''), xlab = "Reporting Group",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 0.6), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  #abline(h = 1/15, lwd = 4)
  
  for(i in 1:15) {
    segments(x0 = ProofPlot[1, i] - 0.5, x1 = ProofPlot[5, i] + 0.5, y0 = MixtureProofTestProportions[i, RG], y1 = MixtureProofTestProportions[i, RG], lwd = 3)
  }
  
  legend(x = max(ProofPlot), y = 0.5, legend = PCGroups15, fill = Colors15, bty = "n", xpd = TRUE)
  
  text(x = ProofPlot[3, ], y = -0.03, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.7)
}

dir.create("BAYES/Mixture Proof Tests/loci24Kodiak/Figures")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### loci24 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

loci24 <- dget(file = "OLD 746 Collections/objects/loci24.txt")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 100% Proof Tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Just doing Frazer and Ayakulik, as these were the "worst" performing 100% proof tests
dir.create("BAYES/100% Proof Tests/loci24")

KMA473PopsGroups15Repeated100ProofTests

#~~~~~~~~~~~~~~~~~~
## Get Proof Objects
# Frazer
for(proof in grep(pattern = "Frazer", x = KMA473PopsGroups15Repeated100ProofTests, value = TRUE)) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

# Ayakulik
for(proof in grep(pattern = "Ayakulik", x = KMA473PopsGroups15Repeated100ProofTests, value = TRUE)) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

#~~~~~~~~~~~~~~~~~~
## Create BAYES files
for (proof in paste(rep(x = c("Ayakulik", "Frazer"), each = 5), 1:5, sep = "")) {
  ReProofTestKS.GCL(sillyvec = KMA473Pops.named, loci = loci24, groupnames = Groups15.nospace, groupvec = KMA473PopsGroupVec15, 
                    ProofTestIDs = get(paste(proof, "Proof", sep = "")), prefix = proof, dir = "BAYES/100% Proof Tests/loci24", 
                    suffix = "", nreps = 40000, nchains = 1, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits[, 1, drop = FALSE], 
                    type = "BAYES", thin = c(1, 1, 100), switches = "F T F T T T F")
}


# Create BAYES Output folders
dir.create("BAYES/100% Proof Tests/loci24/BAYES.output")
sapply(paste(rep(x = c("Frazer", "Ayakulik"), each = 5), 1:5, sep = ""), function(proof) {dir <- paste(getwd(), "/BAYES/100% Proof Tests/loci24/BAYES.output/", proof, sep = "")
                                                                                          dir.create(dir) })

#~~~~~~~~~~~~~~~~~~
## Summarize Results
FrazerAyakulikloci24Repeated100ProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, maindir = "BAYES/100% Proof Tests/loci24/BAYES.output", 
                                                                                         mixvec = paste(rep(x = c("Frazer", "Ayakulik"), each = 5), 1:5, sep = ""), prior = "", ext = "RGN", nchains = 1, burn = 0.5, 
                                                                                         alpha = 0.1, PosteriorOutput = TRUE); beep(5)
str(FrazerAyakulikloci24Repeated100ProofTestsEstimates)

dir.create("Estimates objects/loci24")
dput(x = list(Stats = FrazerAyakulikloci24Repeated100ProofTestsEstimates$Stats[1:5], Output = FrazerAyakulikloci24Repeated100ProofTestsEstimates$Output[1:5]), 
     file = "Estimates objects/loci24/Frazer100ProofTestEstimates.txt")
dput(x = list(Stats = FrazerAyakulikloci24Repeated100ProofTestsEstimates$Stats[6:10], Output = FrazerAyakulikloci24Repeated100ProofTestsEstimates$Output[6:10]), 
     file = "Estimates objects/loci24/Ayakulik100ProofTestEstimates.txt")

#~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~

library(gplots)
library(plotrix)

Groups15.nospace

ProofTest100.SampleSize

ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(5.1, 4.4, 2.1, 1.1))

sapply(c("Frazer", "Ayakulik"), function(RG) {
  par(xpd = FALSE)
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {FrazerAyakulikloci24Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {FrazerAyakulikloci24Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {FrazerAyakulikloci24Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", ProofTest100.SampleSize[RG], ")", sep = ''), xlab = "",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 1), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  abline(h = 0.9, lwd = 2)
  text(x = ProofPlot[3, ], y = -0.065, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.8)
  mtext(text = "Reporting Group", side = 1, line = 3.5)
} )

dir.create("BAYES/100% Proof Tests/loci24/Figures")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Fishery Scenario Tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create("BAYES/Mixture Proof Tests/loci24")

KMA473PopsGroups15RepeatedMixProofTests

#~~~~~~~~~~~~~~~~~~
## Get Proof Objects
for(proof in KMA473PopsGroups15RepeatedMixProofTests) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("MixtureProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

#~~~~~~~~~~~~~~~~~~
## Create BAYES files
for (proof in KMA473PopsGroups15RepeatedMixProofTests){
  ReProofTestKS.GCL(sillyvec = KMA473Pops, loci = loci24, groupnames = Groups15Short, groupvec = KMA473PopsGroupVec15, 
                    ProofTestIDs = get(paste(proof, "Proof", sep = "")), prefix = proof, dir = "BAYES/Mixture Proof Tests/loci24", 
                    suffix = "", nreps = 40000, nchains = 1, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits[, 1, drop = FALSE], 
                    type = "BAYES", thin = c(1, 1, 100), switches = "F T F T T T F")
}

# Create BAYES Output folders
dir.create("BAYES/Mixture Proof Tests/loci24/BAYES.output")
invisible(sapply(KMA473PopsGroups15RepeatedMixProofTests, function(proof) {dir <- paste(getwd(), "/BAYES/Mixture Proof Tests/loci24/BAYES.output/", proof, sep = "")
                                                                           dir.create(dir) } ))


#~~~~~~~~~~~~~~~~~~
## Summarize Results
KMA473PopsGroups15loci24RepeatedMixProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, 
                                                                                             maindir = "BAYES/Mixture Proof Tests/loci24/BAYES.output", 
                                                                                             mixvec = KMA473PopsGroups15RepeatedMixProofTests, prior = "", 
                                                                                             ext = "RGN", nchains = 1, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE); beep(5)
str(KMA473PopsGroups15loci24RepeatedMixProofTestsEstimates)

dput(x = KMA473PopsGroups15loci24RepeatedMixProofTestsEstimates, file = "Estimates objects/loci24/KMA473PopsGroups15loci24RepeatedMixProofTestsEstimates.txt")
dput(x = KMA473PopsGroups15loci24RepeatedMixProofTestsEstimates$Stats, file = "Estimates objects/loci24/KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~

library(gplots)
library(plotrix)

KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci24/KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats.txt")

str(KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats)

Groups15.nospace

MixtureProofTest.SampleSize.Final <- round(apply(MixtureProofTest.SampleSize, 1, sum), 0 )

MixtureProofTestProportions <- read.table(file = "FisheryProofTestScenarios.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = Groups15Short)
MixtureProofTestProportions <- MixtureProofTestProportions[, -1]
max(MixtureProofTestProportions)

# PCGroups15 <- c("West of Chignik", "Black Lake", "Chignik Lake", "Upper Station / Akalura", "Frazer", "Ayakulik", "Karluk", "Uganik", "Northwest Kodiak", "Afognak", "Eastside Kodiak", "Saltery", "Cook Inlet", "Prince William Sound", "South of Cape Suckling")
# dput(x = PCGroups15, file = "Objects/PCGroups15.txt")
PCGroups15 <- dget(file = "Objects/PCGroups15.txt")
ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(3.1, 4.4, 2.1, 10.1))
par(xpd = FALSE)
FisheryProofTestScenarioNames <- dget(file = "Objects/FisheryProofTestScenarioNames.txt")
# RG = "flat"
for(RG in FisheryProofTestScenarioNames) {
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", MixtureProofTest.SampleSize.Final[RG], ")", sep = ''), xlab = "Reporting Group",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 0.6), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  #abline(h = 1/15, lwd = 4)
  
  for(i in 1:15) {
    segments(x0 = ProofPlot[1, i] - 0.5, x1 = ProofPlot[5, i] + 0.5, y0 = MixtureProofTestProportions[i, RG], y1 = MixtureProofTestProportions[i, RG], lwd = 3)
  }
  
  legend(x = max(ProofPlot), y = 0.5, legend = PCGroups15, fill = Colors15, bty = "n", xpd = TRUE)
  
  text(x = ProofPlot[3, ], y = -0.03, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.7)
}

dir.create("BAYES/Mixture Proof Tests/loci24/Figures")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### loci22 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

loci22 <- dget(file = "OLD 746 Collections/objects/loci22.txt")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 100% Proof Tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Just doing Frazer and Ayakulik, as these were the "worst" performing 100% proof tests
dir.create("BAYES/100% Proof Tests/loci22")

KMA473PopsGroups15Repeated100ProofTests

#~~~~~~~~~~~~~~~~~~
## Get Proof Objects
# Frazer
for(proof in grep(pattern = "Frazer", x = KMA473PopsGroups15Repeated100ProofTests, value = TRUE)) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

# Ayakulik
for(proof in grep(pattern = "Ayakulik", x = KMA473PopsGroups15Repeated100ProofTests, value = TRUE)) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

#~~~~~~~~~~~~~~~~~~
## Create BAYES files
for (proof in paste(rep(x = c("Ayakulik", "Frazer"), each = 5), 1:5, sep = "")) {
  ReProofTestKS.GCL(sillyvec = KMA473Pops.named, loci = loci22, groupnames = Groups15.nospace, groupvec = KMA473PopsGroupVec15, 
                    ProofTestIDs = get(paste(proof, "Proof", sep = "")), prefix = proof, dir = "BAYES/100% Proof Tests/loci22", 
                    suffix = "", nreps = 40000, nchains = 1, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits[, 1, drop = FALSE], 
                    type = "BAYES", thin = c(1, 1, 100), switches = "F T F T T T F")
}


# Create BAYES Output folders
dir.create("BAYES/100% Proof Tests/loci22/BAYES.output")
sapply(paste(rep(x = c("Frazer", "Ayakulik"), each = 5), 1:5, sep = ""), function(proof) {dir <- paste(getwd(), "/BAYES/100% Proof Tests/loci22/BAYES.output/", proof, sep = "")
                                                                                          dir.create(dir) })

#~~~~~~~~~~~~~~~~~~
## Summarize Results
FrazerAyakulikloci22Repeated100ProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, maindir = "BAYES/100% Proof Tests/loci22/BAYES.output", 
                                                                                   mixvec = paste(rep(x = c("Frazer", "Ayakulik"), each = 5), 1:5, sep = ""), prior = "", ext = "RGN", nchains = 1, burn = 0.5, 
                                                                                   alpha = 0.1, PosteriorOutput = TRUE); beep(5)
str(FrazerAyakulikloci22Repeated100ProofTestsEstimates)

dir.create("Estimates objects/loci22")
dput(x = list(Stats = FrazerAyakulikloci22Repeated100ProofTestsEstimates$Stats[1:5], Output = FrazerAyakulikloci22Repeated100ProofTestsEstimates$Output[1:5]), 
     file = "Estimates objects/loci22/Frazer100ProofTestEstimates.txt")
dput(x = list(Stats = FrazerAyakulikloci22Repeated100ProofTestsEstimates$Stats[6:10], Output = FrazerAyakulikloci22Repeated100ProofTestsEstimates$Output[6:10]), 
     file = "Estimates objects/loci22/Ayakulik100ProofTestEstimates.txt")

#~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~

library(gplots)
library(plotrix)

Groups15.nospace

ProofTest100.SampleSize

ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(5.1, 4.4, 2.1, 1.1))

sapply(c("Frazer", "Ayakulik"), function(RG) {
  par(xpd = FALSE)
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {FrazerAyakulikloci22Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {FrazerAyakulikloci22Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {FrazerAyakulikloci22Repeated100ProofTestsEstimates$Stats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", ProofTest100.SampleSize[RG], ")", sep = ''), xlab = "",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 1), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  abline(h = 0.9, lwd = 2)
  text(x = ProofPlot[3, ], y = -0.065, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.8)
  mtext(text = "Reporting Group", side = 1, line = 3.5)
} )

dir.create("BAYES/100% Proof Tests/loci22/Figures")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Fishery Scenario Tests
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create("BAYES/Mixture Proof Tests/loci22")

KMA473PopsGroups15RepeatedMixProofTests

#~~~~~~~~~~~~~~~~~~
## Get Proof Objects
for(proof in KMA473PopsGroups15RepeatedMixProofTests) {
  assign(x = paste(proof, "Proof", sep = ""), value = dget(file = paste("MixtureProofTests objects/", proof, "Proof.txt", sep = "")), pos = 1)
}; rm(proof)

#~~~~~~~~~~~~~~~~~~
## Create BAYES files
for (proof in KMA473PopsGroups15RepeatedMixProofTests){
  ReProofTestKS.GCL(sillyvec = KMA473Pops, loci = loci22, groupnames = Groups15Short, groupvec = KMA473PopsGroupVec15, 
                    ProofTestIDs = get(paste(proof, "Proof", sep = "")), prefix = proof, dir = "BAYES/Mixture Proof Tests/loci22", 
                    suffix = "", nreps = 40000, nchains = 1, priorvec = KMA473Pops15FlatPrior, initmat = KMA473PopsInits[, 1, drop = FALSE], 
                    type = "BAYES", thin = c(1, 1, 100), switches = "F T F T T T F")
}

# Create BAYES Output folders
dir.create("BAYES/Mixture Proof Tests/loci22/BAYES.output")
invisible(sapply(KMA473PopsGroups15RepeatedMixProofTests, function(proof) {dir <- paste(getwd(), "/BAYES/Mixture Proof Tests/loci22/BAYES.output/", proof, sep = "")
                                                                           dir.create(dir) } ))


#~~~~~~~~~~~~~~~~~~
## Summarize Results
KMA473PopsGroups15loci22RepeatedMixProofTestsEstimates <- CustomCombineBAYESOutput.GCL(groupvec = seq(Groups15.nospace), groupnames = Groups15.nospace, 
                                                                                       maindir = "BAYES/Mixture Proof Tests/loci22/BAYES.output", 
                                                                                       mixvec = KMA473PopsGroups15RepeatedMixProofTests, prior = "", 
                                                                                       ext = "RGN", nchains = 1, burn = 0.5, alpha = 0.1, PosteriorOutput = TRUE); beep(5)
str(KMA473PopsGroups15loci22RepeatedMixProofTestsEstimates)

dput(x = KMA473PopsGroups15loci22RepeatedMixProofTestsEstimates, file = "Estimates objects/loci22/KMA473PopsGroups15loci22RepeatedMixProofTestsEstimates.txt")
dput(x = KMA473PopsGroups15loci22RepeatedMixProofTestsEstimates$Stats, file = "Estimates objects/loci22/KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats.txt")

#~~~~~~~~~~~~~~~~~~
## Plot results
#~~~~~~~~~~~~~~~~~~

library(gplots)
library(plotrix)

KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci22/KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats.txt")

str(KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats)

Groups15.nospace

MixtureProofTest.SampleSize.Final <- round(apply(MixtureProofTest.SampleSize, 1, sum), 0 )

MixtureProofTestProportions <- read.table(file = "FisheryProofTestScenarios.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = Groups15Short)
MixtureProofTestProportions <- MixtureProofTestProportions[, -1]
max(MixtureProofTestProportions)

# PCGroups15 <- c("West of Chignik", "Black Lake", "Chignik Lake", "Upper Station / Akalura", "Frazer", "Ayakulik", "Karluk", "Uganik", "Northwest Kodiak", "Afognak", "Eastside Kodiak", "Saltery", "Cook Inlet", "Prince William Sound", "South of Cape Suckling")
# dput(x = PCGroups15, file = "Objects/PCGroups15.txt")
PCGroups15 <- dget(file = "Objects/PCGroups15.txt")
ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)

par(mar = c(3.1, 4.4, 2.1, 10.1))
par(xpd = FALSE)
FisheryProofTestScenarioNames <- dget(file = "Objects/FisheryProofTestScenarioNames.txt")
# RG = "flat"
for(RG in FisheryProofTestScenarioNames) {
  ProofPlot <- barplot2(height = t(sapply(1:5, function(i) {KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "median"]} )),
                        ci.l = t(sapply(1:5, function(i) {KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "5%"]} )),
                        ci.u = t(sapply(1:5, function(i) {KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats[[paste(RG, i, sep = "")]][, "95%"]} )),
                        names.arg = NULL, col = rep(Colors15, each = 5), main = paste(RG, " (n = ", MixtureProofTest.SampleSize.Final[RG], ")", sep = ''), xlab = "Reporting Group",
                        ylab = "Proportion", plot.ci = TRUE, beside = TRUE, ylim = c(0, 0.6), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5))
  abline(h = 0)
  #abline(h = 1/15, lwd = 4)
  
  for(i in 1:15) {
    segments(x0 = ProofPlot[1, i] - 0.5, x1 = ProofPlot[5, i] + 0.5, y0 = MixtureProofTestProportions[i, RG], y1 = MixtureProofTestProportions[i, RG], lwd = 3)
  }
  
  legend(x = max(ProofPlot), y = 0.5, legend = PCGroups15, fill = Colors15, bty = "n", xpd = TRUE)
  
  text(x = ProofPlot[3, ], y = -0.03, labels = Groups15Short, srt = 45, xpd = NA, cex = 0.7)
}

dir.create("BAYES/Mixture Proof Tests/loci22/Figures")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### SEDM Style Plot of Mixture Proof Tests ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~
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
plot.ci.extract <- function(x, scenario, group, stat.col) {
  list.index <- grep(pattern = scenario, x = names(x))  
  sapply(list.index, function(i) {x[[i]][group, stat.col]} )
}

# Inputs for extractor function
str(KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats)
FisheryProofTestScenarioNames
Groups15.nospace
c("5%", "median", "95%")

#~~~~~~~~~~~~~~~~~~
# Create a figure wrapup funtion
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")

library('gplots')
library('devEMF')

PlotProof.GCL <- function(x, scenario, groups, filedir, percent = TRUE) {
  
  MixtureProofTestProportions <- t(dget(file = "Objects/MixtureProofTestProportions.txt"))
  
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
  layout(mat = FiveByThreeMatrix, widths = c(0.3, rep(1, 3)), heights = c(1, 1, 1, 1, 1, 0.4))
  
  
  # Plot 1 - y-axis title
  par(mar = c(0, 0, 0, 0))
  par(srt = 90)
  plot.new()
  text(labels = y.label, x = 0.2, y = 0.5, adj = c(0.5, 0.5), font = 1, cex = 2)
  
  # Plot 2-16 - Groups
  sapply(seq(groups), function(i) {
    par(mar = c(0, 0, 0, 0))
    par(srt = 0)
    if(i %in% c(1, 4, 7, 10)) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 5.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", cex = 1.5)
    }
    
    if(i == 13) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 5.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', cex = 1.5)
    }
    
    if(i %in% c(14, 15)) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 5.5), ylim = c(-0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', yaxt = "n", cex = 1.5)
    }
    
    if(i %in% c(2, 3, 5, 6, 8, 9, 11, 12)) {
      plotCI(x = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "median"),
             ui = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "95%"),
             li = plot.ci.extract(x = x, scenario = scenario, group = groups[i], stat.col = "5%"),
             xlim = c(0.5, 5.5), ylim = c(0, 0.65 * percent.multiplier), gap = 0, pch = 16, ylab = '', xlab = '', xaxt = "n", yaxt = "n", cex = 1.5)
    }
    
    abline(h = MixtureProofTestProportions[i, scenario], col = "red", lwd = 2)
    abline(h = 0, lwd = 1, lty = 2)
    text(PCGroups15[i], x = 5, y = 0.6 * percent.multiplier, adj = 1, cex = 1)
  } )
  
  # Plot 17 - Blank Corner
  plot.new()
  
  ## Plot 18 - x-axis title
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(labels = "Replicate", x = 0.5, y = 0.25, adj = c(0.5, 0.5), font = 1, cex = 2)
  
  dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# loci89
sapply(FisheryProofTestScenarioNames, function(scenario){
  PlotProof.GCL(x = dget(file = "Estimates objects/loci89/KMA473PopsGroups15loci89RepeatedMixProofTestsEstimatesStats.txt"), 
                scenario = scenario, groups = Groups15Short,
                filedir = "BAYES/Mixture Proof Tests/loci89/Figures",
                percent = TRUE)
})

# loci46
sapply(FisheryProofTestScenarioNames, function(scenario){
  PlotProof.GCL(x = dget(file = "Estimates objects/loci46/KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats.txt"), 
                scenario = scenario, groups = Groups15.nospace,
                filedir = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/BAYES/Mixture Proof Tests/loci46/Figures")
})

# loci22
sapply(FisheryProofTestScenarioNames, function(scenario){
  PlotProof.GCL(x = dget(file = "Estimates objects/loci22/KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats.txt"), 
                scenario = scenario, groups = Groups15.nospace,
                filedir = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/BAYES/Mixture Proof Tests/loci22/Figures")
})

# loci24
sapply(FisheryProofTestScenarioNames, function(scenario){
  PlotProof.GCL(x = dget(file = "Estimates objects/loci24/KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats.txt"), 
                scenario = scenario, groups = Groups15.nospace,
                filedir = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/BAYES/Mixture Proof Tests/loci24/Figures")
})

# loci24Kodiak
sapply(FisheryProofTestScenarioNames, function(scenario){
  PlotProof.GCL(x = dget(file = "Estimates objects/loci24Kodiak/KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats.txt"), 
                scenario = scenario, groups = Groups15.nospace,
                filedir = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/BAYES/Mixture Proof Tests/loci24Kodiak/Figures")
})

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create one figure with all 5 markersets!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Inputs
# x, scenario, groups, filedir


x.loci89 <- dget(file = "Estimates objects/loci89/KMA473PopsGroups15loci89RepeatedMixProofTestsEstimatesStats.txt")
for(k in seq(x.loci89)) {rownames(x.loci89[[k]]) <- Groups15.nospace}  # Need to standardize dimnames across estimate lists
x.loci46 <- dget(file = "Estimates objects/loci46/KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats.txt")
x.loci22 <- dget(file = "Estimates objects/loci22/KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats.txt")
x.loci24 <- dget(file = "Estimates objects/loci24/KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats.txt")
x.loci24Kodiak <- dget(file = "Estimates objects/loci24Kodiak/KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats.txt")
scenario <- FisheryProofTestScenarioNames[1]
groups <- Groups15.nospace


#~~~~~~~~~~~~~~~~~~
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
plot.ci.extract <- function(x, scenario, group, stat.col) {
  list.index <- grep(pattern = scenario, x = names(x))  
  sapply(list.index, function(i) {x[[i]][group, stat.col]} )
}

# Inputs for extractor function
# str(KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats)
# FisheryProofTestScenarioNames
# Groups15.nospace
# c("5%", "median", "95%")

#~~~~~~~~~~~~~~~~~~
# Create a figure wrapup funtion
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")

library('gplots')
library('devEMF')

percent = TRUE

col = c("black", "grey70")
cex.pch = 1.5
ci.lwd = 2
sfrac = 0.02


for(scenario in FisheryProofTestScenarioNames) {
  
  x.loci89 <- dget(file = "Estimates objects/loci89/KMA473PopsGroups15loci89RepeatedMixProofTestsEstimatesStats.txt")
  for(k in seq(x.loci89)) {rownames(x.loci89[[k]]) <- Groups15.nospace}  # Need to standardize dimnames across estimate lists
  x.loci46 <- dget(file = "Estimates objects/loci46/KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats.txt")
  i <- which(scenario == FisheryProofTestScenarioNames)
  
  if(percent){
    percent.multiplier <- 100
    x.loci89 <- sapply(x.loci89, function(rpt) {rpt[, 1:5] * percent.multiplier}, simplify = FALSE)
    x.loci46 <- sapply(x.loci46, function(rpt) {rpt[, 1:5] * percent.multiplier}, simplify = FALSE)
    file = paste("BAYES/Mixture Proof Tests/Master.", i, scenario, ".Percent.emf", sep = '')
    y.label = "Percent"
  } else {
    percent.multiplier <- 1
    file = paste("BAYES/Mixture Proof Tests/Master.", i, scenario, ".emf", sep = '')
    y.label = "Proportion"
  }
  
  MixtureProofTestProportions <- t(dget(file = "Objects/MixtureProofTestProportions.txt")) * percent.multiplier
    
  emf(file = file, width = 6.5, height = 8, family = "Times")
  #png(file = paste(filedir, "/", scenario, ".png", sep = ""), width = 6.5, height = 7.5, units = "in", res = 1200, pointsize = 4, family = "Times")
  layout(mat = FiveByThreeMatrix, widths = c(0.3, rep(1, 3)), heights = c(1, 1, 1, 1, 1, 0.4))
  
  # Plot 1 - y-axis title
  par(mar = c(0, 0, 0, 0))
  par(srt = 90)
  plot.new()
  text(labels = y.label, x = 0.2, y = 0.5, adj = c(0.5, 0.5), font = 1, cex = 2)
    
  # Plot 2-16 - Groups
  sapply(seq(groups), function(j) {
    
    y <- c(t(matrix(data = c(plot.ci.extract(x = x.loci89, scenario = scenario, group = groups[j], stat.col = "median"),
                             plot.ci.extract(x = x.loci46, scenario = scenario, group = groups[j], stat.col = "median")),
                    nrow = 5)))
    ui <- c(t(matrix(data = c(plot.ci.extract(x = x.loci89, scenario = scenario, group = groups[j], stat.col = "95%"),
                              plot.ci.extract(x = x.loci46, scenario = scenario, group = groups[j], stat.col = "95%")),
                     nrow = 5)))
    li <- c(t(matrix(data = c(plot.ci.extract(x = x.loci89, scenario = scenario, group = groups[j], stat.col = "5%"),
                              plot.ci.extract(x = x.loci46, scenario = scenario, group = groups[j], stat.col = "5%")),
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
    
    abline(h = MixtureProofTestProportions[j, scenario], col = "red", lwd = 2)
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
#### Final markerset for KMA Sockeye Mixtures loci46 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(x = sort(unlist(strsplit(loci46, split = "\\."))), quote = FALSE, row.names = FALSE, col.names = FALSE,
            file = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Mixtures/48SNP loci for KMA Sockeye Mixtures S155.txt")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Create Baseline with 46 Loci ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA47346Baseline <- CreateBaseline.GCL(sillyvec = KMA473Pops, loci = loci46, dir = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/BAYES/Baseline",
                                       basename = "KMA473Pops46Markers", type = "BAYES")
dput(x = KMA47346Baseline, file = "Objects/KMA47346Baseline.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Stats for Abstract ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Number of individuals in final baseline
pop.samp.size <- pop.samp.size
sum(pop.samp.size)  # 65,332

# Average number of individuals per population
mean(pop.samp.size)  # 138
median(pop.samp.size)  # 95
hist(pop.samp.size)
range(pop.samp.size)  # 40-504
sort(pop.samp.size, decreasing = TRUE)  # 40-504
sum(pop.samp.size < 75)  # 17 pops below minimum sample size
table(Groups15[KMA473PopsGroupVec15[pop.samp.size < 75]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 100% Proof Tests
# Average correct allocations
loci <- "loci46"

StatsObject <- dget(file = paste("Estimates objects/", loci, "/KMA473PopsGroups15", loci, "Repeated100ProofTestsEstimatesStats.txt", sep = ''))
# Mean across RGs across repeats of the median from each posterior
round(mean(apply(sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "median"]} )} ), 2, mean)) * 100, 1)
# Range of medians
round(range(sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "median"]} )} )) * 100, 1)
# How many repeats didn't meet a median of 0.9
sum(sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "median"]} )} ) < 0.9)
# How many repeats total
length(sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "median"]} )} ))
# Mean for RGs across repeats of the median from each posterior
round(apply(sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "median"]} )} ), 2, mean) * 100, 1)
# Median allocations for 2 Frazer repeats that failed
lapply(StatsObject[which(sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "median"]} )} ) < 0.9)], function(rpt) {round(rpt[5:6, "median"] * 100, 1)} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fishery Scenario Tests
# Bias and Root Mean Square Error

# Get Fishery Scenario Stats for loci89
MixtureProofTestProportions <- read.table(file = "FisheryProofTestScenarios.txt", header = TRUE, sep = "\t", as.is = TRUE, row.names = NULL)

KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci89/KMA473PopsGroups15loci89RepeatedMixProofTestsEstimatesStats.txt"); beep(4)
str(KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats, max.level = 1)

# View median
KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats$flat1[, "median"]

# Extract median values from all flat mixtures
medians.loci89.flat <- sapply(KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats[31:35], function(proof) {proof[, "median"]})

# Calculate RMSE for flat mixtures
sqrt(rowSums((medians.loci89.flat - MixtureProofTestProportions[, "flat"])^2) / 5)

## Need to create function
# Inputs
stats = KMA473PopsGroups15RepeatedMixProofTestsEstimatesStats
scenario = "flat"
proportions = t(dget(file = "Objects/MixtureProofTestProportions"))
estimator = "median"

#~~~~~~~~~~~~~~~~~~
# Creating a function for Bias and RMSE (percent not proportions)
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


#~~~~~~~~~~~~~~~~~~
### Tabling Output
## Creating a function to table

MasterProofTableMaker <- function(one.hundred.percent = TRUE, loci = "loci89", BiasRMSElist, table.file, percent = TRUE) {
  
  require(xlsx)
  
  if(percent) {
    digits = 1
  } else {
    digits = 2
  }
  
  if(one.hundred.percent) {
    x <- NULL
    x.mat <- matrix(data = seq(PCGroups15), ncol = 3, byrow = TRUE)
    scenario <- Groups15.nospace
    
    for(i in seq(nrow(x.mat))){
      x <- rbind(x, cbind('', PCGroups15[x.mat[i, 1]], '', '', '', '', PCGroups15[x.mat[i, 2]], '', '', '', '', PCGroups15[x.mat[i, 3]], '', '', ''),
                 cbind("Reporting Group", "Average", "Bias", "RMSE", "CI Width", "", "Average", "Bias", "RMSE", "CI Width", "", "Average", "Bias", "RMSE", "CI Width"),
                 cbind(PCGroups15,
                       formatC(x = BiasRMSElist[[scenario[x.mat[i, 1]]]], digits = digits, format = "f"),
                       rep('', length(Groups15)),
                       formatC(x = BiasRMSElist[[scenario[x.mat[i, 2]]]], digits = digits, format = "f"),
                       rep('', length(Groups15)),
                       formatC(x = BiasRMSElist[[scenario[x.mat[i, 3]]]], digits = digits, format = "f"))
      )
    }
    suppressWarnings(write.xlsx(x = x, file = table.file, sheetName = paste("100 Proof Summary", loci), col.names = FALSE, row.names = FALSE, append = TRUE))
    print(ProofTest100.SampleSize)
  } else {
    x <- NULL
    x.mat <- suppressWarnings(matrix(data = seq(FisheryProofTestScenarioNames), ncol = 2, byrow = TRUE))
    scenario <- FisheryProofTestScenarioNames
    proportions <- dget(file = "Objects/MixtureProofTestProportions.txt")
    if(percent) {proportions <- proportions * 100}
    
    for(i in seq(nrow(x.mat))){
      x <- rbind(x, cbind('', paste("Hypothetical", PCFisheryProofTestScenarioNames[x.mat[i, 1]], "Scenario"), '', '', '', '', '', paste("Hypothetical", PCFisheryProofTestScenarioNames[x.mat[i, 2]], "Scenario"), '', '', '', ''),
                 cbind("Reporting Group", "True", "Average", "Bias", "RMSE", "CI Width", "", "True", "Average", "Bias", "RMSE", "CI Width"),
                 cbind(PCGroups15,
                       formatC(x = proportions[x.mat[i, 1], ], digits = digits, format = "f"),
                       formatC(x = BiasRMSElist[[scenario[x.mat[i, 1]]]], digits = digits, format = "f"),
                       rep('', length(Groups15)),
                       formatC(x = proportions[x.mat[i, 2], ], digits = digits, format = "f"),
                       formatC(x = BiasRMSElist[[scenario[x.mat[i, 2]]]], digits = digits, format = "f"))
      )
    }
    suppressWarnings(write.xlsx(x = x, file = table.file, sheetName = paste("Fishery Proof Summary", loci), col.names = FALSE, row.names = FALSE, append = TRUE))
    print(apply(MixtureProofTest.SampleSize, 1, sum))
  }
}


## 100% Proof
# Get Inputs
Groups15.nospace <- dget(file = "Objects/Groups15.nospace.txt")
ProofTestProportions <- array(data = 0, dim = rep(length(Groups15.nospace), 2), dimnames = list(Groups15Short, Groups15.nospace))
diag(ProofTestProportions) <- 1

# loci89, all scenarios
KMA473PopsGroups15loci89Repeated100ProofTestsEstimatesStats <- dget(file = "Estimates objects/loci89/KMA473PopsGroups15loci89Repeated100ProofTestsEstimatesStats.txt")
Proof.BiasRMSE.89loci.list <- sapply(Groups15.nospace, function(scenario) {BiasRMSE.GCL(stats = KMA473PopsGroups15loci89Repeated100ProofTestsEstimatesStats, scenario = scenario, proportions = ProofTestProportions, estimator = "median")}, simplify = FALSE )
Proof.BiasRMSE.89loci.mat <- t(sapply(Groups15.nospace, function(i) {Proof.BiasRMSE.89loci.list[[i]][i, ]} ))
round(Proof.BiasRMSE.89loci.mat[, 1], 1)
MasterProofTableMaker(one.hundred.percent = TRUE, loci = "loci89", BiasRMSElist = Proof.BiasRMSE.89loci.list, table.file = "Estimates tables/KMA473PopsGroups15RepeatedProofTestsTables.xlsx")


# loci46, all scenarios
KMA473PopsGroups15loci46Repeated100ProofTestsEstimatesStats <- dget(file = "Estimates objects/loci46/KMA473PopsGroups15loci46Repeated100ProofTestsEstimatesStats.txt")
Proof.BiasRMSE.46loci.list <- sapply(Groups15.nospace, function(scenario) {BiasRMSE.GCL(stats = KMA473PopsGroups15loci46Repeated100ProofTestsEstimatesStats, scenario = scenario, proportions = ProofTestProportions, estimator = "median")}, simplify = FALSE )
Proof.BiasRMSE.46loci.mat <- t(sapply(Groups15.nospace, function(i) {Proof.BiasRMSE.46loci.list[[i]][i, ]} ))
MasterProofTableMaker(one.hundred.percent = TRUE, loci = "loci46", BiasRMSElist = Proof.BiasRMSE.46loci.list, table.file = "Estimates tables/KMA473PopsGroups15RepeatedProofTestsTables.xlsx")


# Compare RMSE between loci 89 and 46
plot(Proof.BiasRMSE.46loci.mat[, "rmse"], pch = 16, ylim = c(0, 0.12), cex = 2, xlab = "Reporting Group", ylab = "RMSE", bty = 'n', axes = FALSE)
points(Proof.BiasRMSE.89loci.mat[, "rmse"], pch = 16, col = 2, cex = 2)
axis(side = 2)
axis(side = 1, at = 1:15, labels = rep("", 15), pos = 0)
text(x = 1:15, y = -0.009, labels = Groups15Short, srt = 45, xpd = TRUE, pos = 1, adj = c(1, NA))
legend("topright", legend = c("loci46", "loci89"), pch = 16, col = c(1, 2), bty = 'n', cex = 2)

# Compare RMSE between loci 89 and 46
plot(Proof.BiasRMSE.46loci.mat[, "bias"], pch = 16, ylim = c(-0.1, 0), cex = 2, xlab = "Reporting Group", ylab = "Bias", bty = 'n', axes = FALSE)
points(Proof.BiasRMSE.89loci.mat[, "bias"], pch = 16, col = 2, cex = 2)
axis(side = 2)
axis(side = 1, at = 1:15, labels = rep("", 15), pos = 0)
text(x = 1:15, y = 0.01, labels = Groups15Short, srt = 45, xpd = TRUE, pos = 1, adj = c(0, NA))
legend("bottomright", legend = c("loci46", "loci89"), pch = 16, col = c(1, 2), bty = 'n', cex = 2)

#~~~~~~~~
## Fishery Scenarios
# Get Inputs
MixtureProofTestProportions <- t(dget(file = "Objects/MixtureProofTestProportions.txt"))
# loci89, all scenarios
KMA473PopsGroups15loci89RepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci89/KMA473PopsGroups15loci89RepeatedMixProofTestsEstimatesStats.txt"); beep(4)
Fishery.BiasRMSE.89loci.list <- sapply(FisheryProofTestScenarioNames, function(scenario) {BiasRMSE.GCL(stats = KMA473PopsGroups15loci89RepeatedMixProofTestsEstimatesStats, scenario = scenario, proportions = MixtureProofTestProportions, estimator = "median")}, simplify = FALSE )
# across all scenarios
x.array.89 <- array(unlist(Fishery.BiasRMSE.89loci.list), dim = c(15,4,7), dimnames = list(Groups15.nospace, c("average", "bias", "rmse", "ci.width"), FisheryProofTestScenarioNames))
round(apply(X = x.array.89[, , 1:6], MARGIN = 1:2, FUN = mean)[, 2:4], 1)
round(apply(X = x.array.89[, , 1:6], MARGIN = 1:2, FUN = min)[, 2:4], 1)
round(apply(X = x.array.89[, , 1:6], MARGIN = 1:2, FUN = max)[, 2:4], 1)
round(x.array.89[, 2:4, "flat"], 1)
# Table
MasterProofTableMaker(one.hundred.percent = FALSE, loci = "loci89", BiasRMSElist = Fishery.BiasRMSE.89loci.list, table.file = "Estimates tables/KMA473PopsGroups15RepeatedProofTestsTables.xlsx")

# loci46, all scenarios
KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci46/KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats.txt"); beep(4)
Fishery.BiasRMSE.46loci.list <- sapply(FisheryProofTestScenarioNames, function(scenario) {BiasRMSE.GCL(stats = KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats, scenario = scenario, proportions = MixtureProofTestProportions, estimator = "median")}, simplify = FALSE )
# across all scenarios
x.array.46 <- array(unlist(Fishery.BiasRMSE.46loci.list), dim = c(15,4,7), dimnames = list(Groups15.nospace, c("average", "bias", "rmse", "ci.width"), FisheryProofTestScenarioNames))
round(apply(X = x.array.46[c("Frazer", "Ayakulik"), , 1:6], MARGIN = 1:2, FUN = mean)[, 2:4], 1)
round(apply(X = x.array.46[c("Frazer", "Ayakulik"), , 1:6], MARGIN = 1:2, FUN = min)[, 2:4], 1)
round(apply(X = x.array.46[c("Frazer", "Ayakulik"), , 1:6], MARGIN = 1:2, FUN = max)[, 2:4], 1)
round(x.array.46[, 2:4, "flat"], 1)

round(apply(X = x.array.46[, , 1:6], MARGIN = 1:2, FUN = mean)[, 2:4], 1)

lapply(KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats, function(rpt) {rpt[5:6, 1:5]})
MixtureProofTestProportions

# Table
MasterProofTableMaker(one.hundred.percent = FALSE, loci = "loci46", BiasRMSElist = Fishery.BiasRMSE.46loci.list, table.file = "Estimates tables/KMA473PopsGroups15RepeatedProofTestsTables.xlsx")

# loci22, all scenarios
KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci22/KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats.txt"); beep(4)
Fishery.BiasRMSE.22loci.list <- sapply(FisheryProofTestScenarioNames, function(scenario) {BiasRMSE.GCL(stats = KMA473PopsGroups15loci22RepeatedMixProofTestsEstimatesStats, scenario = scenario, proportions = MixtureProofTestProportions, estimator = "median")}, simplify = FALSE )

# loci24, all scenarios
KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci24/KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats.txt"); beep(4)
Fishery.BiasRMSE.24loci.list <- sapply(FisheryProofTestScenarioNames, function(scenario) {BiasRMSE.GCL(stats = KMA473PopsGroups15loci24RepeatedMixProofTestsEstimatesStats, scenario = scenario, proportions = MixtureProofTestProportions, estimator = "median")}, simplify = FALSE )

# loci24Kodiak, all scenarios
KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats <- dget(file = "Estimates objects/loci24Kodiak/KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats.txt"); beep(4)
Fishery.BiasRMSE.24Kodiakloci.list <- sapply(FisheryProofTestScenarioNames, function(scenario) {BiasRMSE.GCL(stats = KMA473PopsGroups15loci24KodiakRepeatedMixProofTestsEstimatesStats, scenario = scenario, proportions = MixtureProofTestProportions, estimator = "median")}, simplify = FALSE )



# Compare RMSE between markersets
type = "l"
lwd = 5

for(scenario in FisheryProofTestScenarioNames) {
  
  ymax <- max(Fishery.BiasRMSE.89loci.list[[scenario]][, "rmse"], Fishery.BiasRMSE.46loci.list[[scenario]][, "rmse"], Fishery.BiasRMSE.22loci.list[[scenario]][, "rmse"], Fishery.BiasRMSE.24loci.list[[scenario]][, "rmse"], Fishery.BiasRMSE.24Kodiakloci.list[[scenario]][, "rmse"])
  
  plot(NA, pch = 16, xlim = c(1, 19), ylim = c(0, ymax), cex = 2, xlab = "Reporting Group", ylab = "RMSE", bty = 'n', axes = FALSE, main = scenario)
  axis(side = 2)
  axis(side = 1, at = 1:15, labels = rep("", 15), pos = 0)
  text(x = 1:15, y = 0 - ymax/25, labels = Groups15Short, srt = 45, xpd = TRUE, adj = c(1, NA))
  legend(x = 15.5, y = ymax, legend = c("89", "46", "22", "24", "24Kodiak"), pch = 16, col = 1:5, bty = 'n', cex = 1.4)
  points(Fishery.BiasRMSE.89loci.list[[scenario]][, "rmse"], pch = 16, col = 1, cex = 2, type = type, lwd = lwd)
  points(Fishery.BiasRMSE.46loci.list[[scenario]][, "rmse"], pch = 16, col = 2, cex = 2, type = type, lwd = lwd)
  points(Fishery.BiasRMSE.22loci.list[[scenario]][, "rmse"], pch = 16, col = 3, cex = 2, type = type, lwd = lwd)
  points(Fishery.BiasRMSE.24loci.list[[scenario]][, "rmse"], pch = 16, col = 4, cex = 2, type = type, lwd = lwd)
  points(Fishery.BiasRMSE.24Kodiakloci.list[[scenario]][, "rmse"], pch = 16, col = 5, cex = 2, type = type, lwd = lwd)
}

# Compare biases between markersets
type = "o"
lwd = 5

for(scenario in FisheryProofTestScenarioNames) {
  
  ymin <- min(Fishery.BiasRMSE.89loci.list[[scenario]][, "bias"], Fishery.BiasRMSE.46loci.list[[scenario]][, "bias"], Fishery.BiasRMSE.22loci.list[[scenario]][, "bias"], Fishery.BiasRMSE.24loci.list[[scenario]][, "bias"], Fishery.BiasRMSE.24Kodiakloci.list[[scenario]][, "bias"])
  ymax <- max(Fishery.BiasRMSE.89loci.list[[scenario]][, "bias"], Fishery.BiasRMSE.46loci.list[[scenario]][, "bias"], Fishery.BiasRMSE.22loci.list[[scenario]][, "bias"], Fishery.BiasRMSE.24loci.list[[scenario]][, "bias"], Fishery.BiasRMSE.24Kodiakloci.list[[scenario]][, "bias"])
  
  plot(NA, pch = 16, xlim = c(1, 19), ylim = c(ymin, ymax), cex = 2, xlab = "Reporting Group", ylab = "Mean Bias", bty = 'n', axes = FALSE, main = scenario)
  axis(side = 2)
  axis(side = 1, at = 1:15, labels = rep("", 15), pos = 0, tck = -0.01)
  axis(side = 1, at = 1:15, labels = rep("", 15), pos = 0, tck = 0.01)
  segments(x0 = 0, y0 = 0, x1 = 15, y1 = 0)
  text(x = 1:15, y = ymin - ((ymax - ymin) / 25), labels = Groups15Short, srt = 45, xpd = TRUE, adj = c(1, NA))
  legend(x = 15.5, y = ymax, legend = c("89", "46", "22", "24", "24Kodiak"), pch = 16, col = 1:5, bty = 'n', cex = 1.4)
  points(Fishery.BiasRMSE.89loci.list[[scenario]][, "bias"], pch = 16, col = 1, cex = 2, type = type, lwd = lwd)
  points(Fishery.BiasRMSE.46loci.list[[scenario]][, "bias"], pch = 16, col = 2, cex = 2, type = type, lwd = lwd)
  points(Fishery.BiasRMSE.22loci.list[[scenario]][, "bias"], pch = 16, col = 3, cex = 2, type = type, lwd = lwd)
  points(Fishery.BiasRMSE.24loci.list[[scenario]][, "bias"], pch = 16, col = 4, cex = 2, type = type, lwd = lwd)
  points(Fishery.BiasRMSE.24Kodiakloci.list[[scenario]][, "bias"], pch = 16, col = 5, cex = 2, type = type, lwd = lwd)
}

#~~~~~~~~~~~~~~~~~~
# Comparison with WASSIP collections
WASSIP294Pops <- dget(file = "V:\\Analysis\\_WORK_JunkDrawer\\WASSIP\\zzzSockeye\\Baseline\\Objects\\WASSIP294Pops.txt")
WASSIPCollections <- sapply(WASSIP294Pops, function(pop) {unlist(strsplit(x = pop, split = "\\."))} )
KodiakCollections <- sapply(Kodiak57Pops, function(pop) {unlist(strsplit(x = pop, split = "\\."))} )

table(unlist(KodiakCollections) %in% unlist(WASSIPCollections))
table(Kodiak57Pops %in% WASSIP294Pops)

unlist(KodiakCollections)[!unlist(KodiakCollections) %in% unlist(WASSIPCollections)]
unlist(KodiakCollections)[unlist(KodiakCollections) %in% unlist(WASSIPCollections)]

KodiakCollections
table(sapply(KodiakCollections, function(pop) {sum(pop %in% unlist(WASSIPCollections))} ))  # new pops are 0
table(sapply(KodiakCollections, function(pop) {length(table(pop %in% unlist(WASSIPCollections)))} ))  # new collections to pops = 2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figures for FMS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Master plot of all 100% proof tests (a la WASSIP Figure 18)

# Libraries
library(devEMF)
library(gplots)
library(plotrix)

# Inputs
loci <- "loci46"

percent = 100

Groups15.nospace <- dget(file = "Objects/Groups15.nospace.txt")
ProofTest100.SampleSize <- dget(file = "Objects/ProofTest100.SampleSize.txt")
ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)
StatsObject <- dget(file = paste("Estimates objects/", loci, "/KMA473PopsGroups15", loci, "Repeated100ProofTestsEstimatesStats.txt", sep = ''))
StatsObject <- sapply(StatsObject, function(rpt) {rpt[, 1:5] * percent}, simplify = FALSE)
PCGroups15 <- dget(file = "Objects/PCGroups15.txt")

# Create .emf file
emf(file = paste("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/BAYES/100% Proof Tests/", loci, "/Figures/MasterPlotPercent.emf", sep = ''), width = 9.5, height = 5.5, family = "Times")

par(mar = c(2.6, 4.1, 1.1, 9.1))

ProofPlot <- barplot2(height = sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "median"]} )} ),
                      ci.l = sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "5%"]} )} ),
                      ci.u = sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "95%"]} )} ),
                      names.arg = NULL, col = rep(Colors15, each = 5), main = '', xlab = '', ylab = "Percent Correctly Allocated", plot.ci = TRUE, 
                      beside = TRUE, ylim = c(0, 1 * percent), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5)
)
abline(h = 0)
abline(h = 0.9 * percent, lwd = 2)
mtext(text = "Reporting Group", side = 1, line = 1)
legend(x = max(ProofPlot), y = 0.8 * percent, legend = PCGroups15, fill = Colors15, bty = "n", xpd = TRUE)

dev.off()

## Making shorter for CommFish Poster
# Create .emf file
emf(file = paste("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/BAYES/100% Proof Tests/", loci, "/Figures/MasterPlotPercentShort.emf", sep = ''), width = 13.5, height = 5.25, family = "Times")

par(mar = c(2.1, 4.1, 1.1, 12.1))

ProofPlot <- barplot2(height = sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "median"]} )} ),
                      ci.l = sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "5%"]} )} ),
                      ci.u = sapply(Groups15.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "95%"]} )} ),
                      names.arg = NULL, col = rep(Colors15, each = 5), main = '', xlab = '', ylab = '', plot.ci = TRUE, cex.axis = 1.5,
                      beside = TRUE, ylim = c(0, 1 * percent), ci.lwd = 2, ci.color = rep(ci.bar.colors.Groups15, each = 5)
)
abline(h = 0)
abline(h = 0.9 * percent, lwd = 2)
mtext(text = "Reporting Group", side = 1, line = 0.5, cex = 1.5)
mtext(text = "Percentage Correctly Allocated", side = 2, line = 2.5, cex = 1.5)
legend(x = max(ProofPlot), y = 0.9 * percent, legend = PCGroups15, fill = Colors15, bty = "n", xpd = TRUE, cex = 1.3)

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Figures for COMFISH Poster ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Example plot of "results" a la WASSIP and SEDM (barplots with maps)
# Decided to go with late.Karluk scenario as it performed well, and doesn't have a lot of CI fish

# Libraries
library(devEMF)
library(gplots)
library(plotrix)

KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats <- dget("Estimates objects/loci46/KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats.txt")


plot.ci.extract <- function(x, scenario, group, stat.col) {
  list.index <- grep(pattern = scenario, x = names(x))  
  sapply(list.index, function(i) {x[[i]][group, stat.col]} )
}

# # Inputs for extractor function
# str(KMA473PopsGroups15loci46RepeatedMixProofTestsEstimatesStats)
# FisheryProofTestScenarioNames
# Groups15.nospace
# c("5%", "median", "95%")




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
#### Look at attractor/deflector RGs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~
# Comparing Likelihood Profiles and Confusion Matrices
round(diag(KMA473Pops_89loci_Confusion$GroupByGroup), 3)
round(diag(KMA473Pops_46loci_Confusion$GroupByGroup), 3)

round(summary(diag(KMA473Pops_89loci_Confusion$GroupByGroup)), 3)
round(summary(diag(KMA473Pops_46loci_Confusion$GroupByGroup)), 3)











KMAGroups15RepeatedFlatMixProofTestsEstimatesStats
deviation.from.mean.flat <- sapply(Groups15Short, function(RG) {mean(sapply(KMAGroups15RepeatedFlatMixProofTestsEstimatesStats, function(rpt) {rpt[RG, "mean"]} ))} ) - (1 / length(Groups15Short))
percent.deviation.from.mean.flat <- deviation.from.mean.flat / (1 / 15) * 100

par(mar = c(2.1, 6.1, 2.1, 4.1))
plot(percent.deviation.from.mean.flat, pch = 16, axes = FALSE, cex = 3, cex.lab = 2, col = Colors15, xlab = '', ylab = '')
abline(h = 0, lwd = 3)
axis(side = 2, cex.axis = 2, lwd = 3)
# text(x = seq(Groups15Short), y = min(deviation.from.mean.flat) - (0.05 * sum(abs(range(deviation.from.mean.flat)))), pos = 1, srt = 45, labels = Groups15Short, xpd = TRUE, cex = 2)
mtext(text = "Percent mean deviation over repreats (%)", side = 2, line = 4, cex = 2)
text(x = seq(Groups15Short), y = percent.deviation.from.mean.flat, srt = 0, labels = Groups15Short, xpd = TRUE, cex = 1.8, pos = c(rep(4, 8), 3, rep(4, 5)), offset = 0.8)
text(x = 2, y = range(percent.deviation.from.mean.flat) + c((0.15 * sum(abs(range(percent.deviation.from.mean.flat)))), -(0.15 * sum(abs(range(percent.deviation.from.mean.flat))))), labels = c("Deflector", "Attractor"), cex = 3, col = c("red", "black"))


sapply(Groups15Short, function(RG) {mean(sapply(KMAGroups15RepeatedFlatMixProofTestsEstimatesStats, function(rpt) {rpt[RG, "sd"] / (1 / length(Groups15Short))} ))} ) * 100
CV.flat <- sapply(Groups15Short, function(RG) {mean(sapply(KMAGroups15RepeatedFlatMixProofTestsEstimatesStats, function(rpt) {rpt[RG, "sd"] / rpt[RG, "mean"]} ))} ) * 100

par(mar = c(7.1, 6.1, 2.1, 4.1))
plot(CV.flat, col = Colors15, lwd = 20, type = "h", bty = "n", axes = FALSE, xlab = '', ylab = '', ylim = c(0, max(CV.flat)))
axis(side = 1, lwd = 3, labels = NA, at = c(0, seq(Groups15Short)), pos = 0)
axis(side = 2, lwd = 3, cex.axis = 2.5)
mtext(text = "CV (%)", side = 2, line = 4, cex = 2)
text(x = seq(Groups15Short), y = -3, pos = 2, srt = 45, labels = Groups15Short, xpd = TRUE, cex = 2, offset = -2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Quick look at Karluk / Frazer / Ayakulik PostQC Collections MDS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline")
source("V:\\DATA\\R_GEN\\GCL Source Scripts\\Functions.GCL.r")
source("V:/WORK/Kyle/R Source Scripts/Functions.GCL_KS.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 96 Loci
## Locus Control
LocusControl <- dget(file = "Objects/OriginalLocusControl.txt")
loci96 <- dget(file = "Objects/loci96.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get PostQC Collections
# FrazerKarlukAyakulik41Collections <- readClipboard()
# dput(x = FrazerKarlukAyakulik41Collections, file = "Objects/FrazerKarlukAyakulik41Collections.txt")
FrazerKarlukAyakulik41Collections <- dget(file = "Objects/FrazerKarlukAyakulik41Collections.txt")

require(beepr)
invisible(sapply(FrazerKarlukAyakulik41Collections, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCCollections/", silly, ".txt", sep = "")), pos = 1)})); beep(2)
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

# FrazerKarlukAyakulik41CollectionsCommonNames <- readClipboard()
# dput(x = FrazerKarlukAyakulik41CollectionsCommonNames, file = "Objects/FrazerKarlukAyakulik41CollectionsCommonNames.txt")
FrazerKarlukAyakulik41CollectionsCommonNames <- dget(file = "Objects/FrazerKarlukAyakulik41CollectionsCommonNames.txt")

Groups3 <- c("Frazer", "Ayakulik", "Karluk")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Combine Loci
CombineLoci.GCL(sillyvec = FrazerKarlukAyakulik41Collections, markerset = c("One_CO1", "One_Cytb_17", "One_Cytb_26"), update = TRUE); beep(2)
CombineLoci.GCL(sillyvec = FrazerKarlukAyakulik41Collections, markerset = c("One_Tf_ex10-750", "One_Tf_ex3-182"), update = TRUE); beep(2)

loci89

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MDS 89 loci
# gcl2Genepop.GCL(sillyvec = FrazerKarlukAyakulik41Collections, loci = loci89, path = "Genepop/FrazerKarlukAyakulik41Collections_89loci.gen", VialNums = TRUE)

require(adegenet)
genind <- read.genepop(file = "Genepop/FrazerKarlukAyakulik41Collections_89loci.gen")

genpop <- genind2genpop(genind)

# AdegenetNei41Col89loci <- dist.genpop(genpop, method = 1, diag = TRUE, upper = TRUE)
# dput(x = AdegenetNei41Col89loci, file = "Trees/AdegenetNei41Col89loci.txt")
AdegenetNei41Col89loci <- dget(file = "Trees/AdegenetNei41Col89loci.txt")
str(AdegenetNei41Col89loci)

require(ape)
Nei41NJtree <- nj(AdegenetNei41Col89loci)
str(Nei41NJtree)
Nei41NJtree$tip.label <- readClipboard()
plot.phylo(x = Nei41NJtree, cex = 0.5, no.margin = TRUE)

library('rgl')

MDS <- cmdscale(as.matrix(AdegenetNei41Col89loci), k = 3)
# MDS <- cmdscale(as.matrix(AdegenetNei41Col89loci), k = 40, eig = TRUE)  # Do with all possible dimensions
# dput(x = MDS, file = "Objects/MDSAdegenetNei41Col89loci.txt")
# dput(x = MDS, file = "Objects/MDSAdegenetNei41Col89loci_alldim.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei41Col89loci.txt")

x <- as.vector(MDS[, 1])   
y <- as.vector(MDS[, 2])
z <- as.vector(MDS[, 3])

FrazerKarlukAyakulik41CollectionsGroupVec3
Colors15 <- dget(file = "Objects/Colors15.txt")


open3d()
par(family = "serif")
plot3d(x, y, z + abs(range(z)[1]), xlab = '', ylab = '', zlab = '', aspect = FALSE, col = Colors15[FrazerKarlukAyakulik41CollectionsGroupVec3], size = 0.5, type = 's', axes = TRUE, box = TRUE, top = TRUE, cex = 1)
box3d()
# plot3d(x, y, z + abs(range(z)[1]), aspect = FALSE, col = "black", size = 0.5, type = 'h', box = TRUE, axes = FALSE, top = FALSE, add = TRUE)  # adds pins to spheres
# texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.8, 0.8), text = seq(FrazerKarlukAyakulik41Collections), font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)
texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.2, 0.2), text = FrazerKarlukAyakulik41Collections, font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)

rgl.snapshot("MDS/MDSAdegenetNei41ColFrazerKarlukAyakulik89loci2.png", fmt="png", top=TRUE )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MDS 88 loci (no One_U1004-183)
grep(pattern = "One_U1004-183", x = loci89)
# gcl2Genepop.GCL(sillyvec = FrazerKarlukAyakulik41Collections, loci = loci89[-54], path = "Genepop/FrazerKarlukAyakulik41Collections_88loci.gen", VialNums = TRUE)

require(adegenet)
genind <- read.genepop(file = "Genepop/FrazerKarlukAyakulik41Collections_88loci.gen")

genpop <- genind2genpop(genind)

# AdegenetNei41Col88loci <- dist.genpop(genpop, method = 1, diag = TRUE, upper = TRUE)
# dput(x = AdegenetNei41Col88loci, file = "Trees/AdegenetNei41Col88loci.txt")
AdegenetNei41Col88loci <- dget(file = "Trees/AdegenetNei41Col88loci.txt")
str(AdegenetNei41Col88loci)

require(ape)
Nei41NJtree <- nj(AdegenetNei41Col88loci)
str(Nei41NJtree)
Nei41NJtree$tip.label <- readClipboard()
plot.phylo(x = Nei41NJtree, cex = 0.5, no.margin = TRUE)

library('rgl')

MDS <- cmdscale(as.matrix(AdegenetNei41Col88loci), k = 3)
# MDS <- cmdscale(as.matrix(AdegenetNei41Col88loci), k = 40, eig = TRUE)  # Do with all possible dimensions
# dput(x = MDS, file = "Objects/MDSAdegenetNei41Col88loci.txt")
# dput(x = MDS, file = "Objects/MDSAdegenetNei41Col88loci_alldim.txt")
MDS <- dget(file = "Objects/MDSAdegenetNei41Col88loci.txt")

x <- as.vector(MDS[, 1])   
y <- as.vector(MDS[, 2])
z <- as.vector(MDS[, 3])

FrazerKarlukAyakulik41CollectionsGroupVec3
Colors15 <- dget(file = "Objects/Colors15.txt")


open3d()
par(family = "serif")
plot3d(x, y, z + abs(range(z)[1]), xlab = '', ylab = '', zlab = '', aspect = FALSE, col = Colors15[FrazerKarlukAyakulik41CollectionsGroupVec3], size = 0.5, type = 's', axes = TRUE, box = TRUE, top = TRUE, cex = 1)
box3d()
# plot3d(x, y, z + abs(range(z)[1]), aspect = FALSE, col = "black", size = 0.5, type = 'h', box = TRUE, axes = FALSE, top = FALSE, add = TRUE)  # adds pins to spheres
# texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.8, 0.8), text = seq(FrazerKarlukAyakulik41Collections), font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)
texts3d(x, y, z + abs(range(z)[1]), adj = c(-0.2, 0.2), text = FrazerKarlukAyakulik41Collections, font = 2, cex = 0.8, add = TRUE, top = TRUE, axes = FALSE)  # adds numbers to points(adj moves the numbers around the points)

rgl.snapshot("MDS/MDSAdegenetNei41ColFrazerKarlukAyakulik88loci2.png", fmt="png", top=TRUE )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Exploring DAPC ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline")
# source("V:\\Analysis\\R files\\Scripts\\PROD\\Functions.GCL.r")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get objects
KMAobjects <- list.files(path = "Objects", recursive = FALSE)
KMAobjects <- KMAobjects[!KMAobjects %in% c("OLD", "Likelihood Profiles", "OriginalLocusControl.txt")]
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "Objects", objct, sep = "/")), pos = 1) })); beep(2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 103 Loci
## Locus Control
LocusControl <- dget(file = "Objects/LocusControl103.txt")
loci103 <- dget(file = "Objects/loci103.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get Populations
KMA473Pops <- dget(file = "Objects/KMA473Pops.txt")
Karluk16Pops <- KMA473Pops[which(KMA473PopsGroupVec15 == 7)]

require(beepr)
invisible(sapply(Karluk16Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPopsloci103/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
length(objects(pattern = "\\.gcl"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Dump a GENEPOP file for genind

loci88diploid <- sort(c(loci89[1:87], "One_Tf_ex3-182"))

gcl2Genepop.GCL(sillyvec = Karluk16Pops, loci = loci88diploid, path = "Genepop/Karluk16Pops_88diploidloci.gen", VialNums = FALSE)

require(adegenet)
genind <- read.genepop(file = "Genepop/Karluk16Pops_88diploidloci.gen")
genind

grp <- find.clusters(x = genind, max.n.clust = 40)
str(grp)

grp$size
sapply(Karluk16Pops, function(silly) {get(paste(silly, ".gcl", sep = ''))$n} )

table(pop(genind), grp$grp)
table.value(table(pop(genind), grp$grp), col.labels = paste("infer", seq(unique(grp$grp))), row.labels = paste("original", seq(unique(pop(genind)))))


grp2 <- find.clusters(x = genind, max.n.clust = 40)
str(grp2)

grp2$size
sapply(Karluk16Pops, function(silly) {get(paste(silly, ".gcl", sep = ''))$n} )

table(pop(genind), grp2$grp)
table.value(table(pop(genind), grp2$grp), col.labels = paste("infer", seq(unique(grp2$grp))), row.labels = paste("original", seq(unique(pop(genind)))))

dapc8 <- dapc(x = genind, grp = grp$grp)
dapc8
scatter(dapc8)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get Populations
KMA473Pops <- dget(file = "Objects/KMA473Pops.txt")
SWKodiak37Pops <- KMA473Pops[which(KMA473PopsGroupVec15 %in% 4:7)]

require(beepr)
invisible(sapply(SWKodiak37Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste(getwd(), "/Raw genotypes/PostQCPopsloci103/", silly, ".txt", sep = "")), pos = 1)} )); beep(2)
length(objects(pattern = "\\.gcl"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Dump a GENEPOP file for genind

loci88diploid <- sort(c(loci89[1:87], "One_Tf_ex3-182"))

gcl2Genepop.GCL(sillyvec = SWKodiak37Pops, loci = loci88diploid, path = "Genepop/SWKodiak37Pops_88diploidloci.gen", VialNums = FALSE)

require(adegenet)
genind <- read.genepop(file = "Genepop/SWKodiak37Pops_88diploidloci.gen")
genind

grp <- find.clusters(x = genind, max.n.clust = 40)
str(grp)

grp$size
sapply(Karluk16Pops, function(silly) {get(paste(silly, ".gcl", sep = ''))$n} )

table(pop(genind), grp$grp)
table.value(table(pop(genind), grp$grp), col.labels = paste("infer", seq(unique(grp$grp))), row.labels = paste("original", seq(unique(pop(genind)))))
sapply(PCGroups15[4:7], function(group) {KMA473Pops[which(KMA473PopsGroupVec15 == which(PCGroups15 == group))]} )


dapc8 <- dapc(x = genind, grp = grp$grp)
dapc8
scatter(dapc8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Ayakulik/Frazer 100% Proof Test re-do with loci46 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create Proof Tests
Groups14.nospace <- c(Groups15.nospace[1:4], "AyakulikFrazer", Groups15.nospace[7:15])


## Ayakulik/Frazer
for(Proof in paste0("AyakulikFrazer", 1:5)) {
  assign(x = paste(Proof, "Proof", sep = ""),
         ProofTest.GCL(sillyvec = KMA473Pops.named, loci = loci46, groupnames = Groups14.nospace, groupvec = KMA473PopsGroupVec14,
                       samplesize = setNames(object = c(rep(0,4), 200, rep(0, 9)), nm = Groups14.nospace), prefix = Proof,
                       dir = "BAYES/100% Proof Tests/loci46", 
                       prprtnl = TRUE, type = "BAYES", suffix = "", nreps = 40000, nchains = 5, priorvec = KMA473Pops14FlatPrior, 
                       initmat = KMA473PopsInits, thin = c(1, 1, 100), switches = "F T F T T T F"))
}; rm(Proof)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Dput proof test objects
objects(pattern = "Proof$")
invisible(sapply(paste0("AyakulikFrazer", 1:5), function(proof) {dput(x = get(paste(proof, "Proof", sep = "")), file = paste("100%ProofTests objects/", proof, "Proof.txt", sep = ""))} ))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create BAYES.output folders
sapply(paste0("AyakulikFrazer", 1:5), function(proof) {dir <- paste(getwd(), "/BAYES/100% Proof Tests/loci46/BAYES.output/", proof, sep = "")
dir.create(dir) })



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Summarize 100% Proof Tests
AyakulikFrazerloci46Repeated100ProofTestsEstimates <- 
  CustomCombineBAYESOutput.GCL(
    groupvec = seq(Groups14.nospace), groupnames = Groups14.nospace, maindir = "BAYES/100% Proof Tests/loci46/BAYES.output", 
    mixvec = paste0("AyakulikFrazer", 1:5), prior = "", ext = "RGN", nchains = 4, burn = 0.5, alpha = 0.1,
    PosteriorOutput = TRUE); beep(2)

dput(x = AyakulikFrazerloci46Repeated100ProofTestsEstimates, file = "Estimates objects/loci46/AyakulikFrazerloci46Repeated100ProofTestsEstimates.txt"); beep(4)
dput(x = AyakulikFrazerloci46Repeated100ProofTestsEstimates$Stats, file = "Estimates objects/loci46/AyakulikFrazerloci46Repeated100ProofTestsEstimatesStats.txt"); beep(4)

str(AyakulikFrazerloci46Repeated100ProofTestsEstimates, max.level = 1)

## Check Gelman-Rubin factor
lapply(AyakulikFrazerloci46Repeated100ProofTestsEstimates$Stats, function(Proof) {Proof[, "GR"]} )
sapply(AyakulikFrazerloci46Repeated100ProofTestsEstimates$Stats, function(Proof) {table(Proof[, "GR"] > 1.2)} )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Master plot of all 100% proof tests (a la WASSIP Figure 18)

# Libraries
library(devEMF)
library(gplots)
library(plotrix)

# Inputs
loci <- "loci46"

percent = 100

Groups15.nospace <- dget(file = "Objects/Groups15.nospace.txt")
ProofTest100.SampleSize <- dget(file = "Objects/ProofTest100.SampleSize.txt")
ci.bar.colors.Groups15 <- setNames(object = c("black", "grey50", rep("black", 3), "grey50", rep("black", 9)), nm = Groups15.nospace)
StatsObject <- dget(file = paste("Estimates objects/", loci, "/KMA473PopsGroups15", loci, "Repeated100ProofTestsEstimatesStats.txt", sep = ''))
StatsObject <- sapply(StatsObject, function(rpt) {rpt[, 1:5] * percent}, simplify = FALSE)

AyakulikFrazerObject <- dget(file = paste0("Estimates objects/", loci, "/AyakulikFrazer", loci, "Repeated100ProofTestsEstimatesStats.txt"))
AyakulikFrazerObject <- sapply(AyakulikFrazerObject, function(rpt) {rpt[, 1:5] * percent}, simplify = FALSE)

StatsObject <- c(StatsObject, AyakulikFrazerObject)

PCGroups15 <- dget(file = "Objects/PCGroups15.txt")
KMA14GroupsPC2 <- dget(file = "V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures/Objects/KMA14GroupsPC2.txt")

# Create .emf file
emf(file = paste("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/BAYES/100% Proof Tests/", loci, "/Figures/MasterPlotPercent14RG.emf", sep = ''), width = 9.5, height = 5.5, family = "Times")

par(mar = c(2.6, 4.1, 1.1, 9.1))

ProofPlot <- barplot2(height = sapply(Groups14.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "median"]} )} ),
                      ci.l = sapply(Groups14.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "5%"]} )} ),
                      ci.u = sapply(Groups14.nospace, function(RG) {sapply(1:5, function(i) {StatsObject[[paste(RG, i, sep = '')]][RG, "95%"]} )} ),
                      names.arg = NULL, col = rep(Colors14, each = 5), main = '', xlab = '', ylab = "Percent Correctly Allocated", plot.ci = TRUE, 
                      beside = TRUE, ylim = c(0, 1 * percent), ci.lwd = 2  # , ci.color = rep(ci.bar.colors.Groups15, each = 5)
)
abline(h = 0)
abline(h = 0.9 * percent, lwd = 2)
mtext(text = "Reporting Group", side = 1, line = 1)
legend(x = max(ProofPlot), y = 0.8 * percent, legend = KMA14GroupsPC2, fill = Colors14, bty = "n", xpd = TRUE)

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Likelihood Profile 31 Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA473Pops31GroupsLikelihoodProfile <- LeaveOneOutDist.GCL(sillyvec = KMA473Pops, loci = loci46, groupvec = KMA473PopsGroupVec31)
KMA473Pops_31Groups_46loci_LikelihoodProfile <- KMA473Pops31GroupsLikelihoodProfile
dput(x = KMA473Pops_31Groups_46loci_LikelihoodProfile, file = "Likelihood Profiles/KMA473Pops_31Groups_46loci_LikelihoodProfile.txt")

KMA473Pops31Groups_46loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = KMA473Pops_31Groups_46loci_LikelihoodProfile, groupnames = KMA31GroupsPC, groupvec = KMA473PopsGroupVec31, sillyvec = KMA473Pops)
dput(x = KMA473Pops31Groups_46loci_Confusion, file = "Objects/KMA473Pops31Groups_46loci_Confusion.txt")

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(KMA473Pops31Groups_46loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "31 Fine Scale Groups, 46 loci", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))


KMA473Pops_31Groups_46loci_LikelihoodProfile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = KMA473Pops, loci = loci46, groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, groupcomps = NULL, ncores = 7)
dput(x = KMA473Pops_31Groups_46loci_LikelihoodProfile_NEW, file = "Likelihood Profiles/KMA473Pops_31Groups_46loci_LikelihoodProfile_NEW.txt")
PlotLikeProfile.GCL(likeprof = KMA473Pops_31Groups_46loci_LikelihoodProfile_NEW, popvec = KMA473Pops, loci = loci46, groupvec = KMA473PopsGroupVec31, groupnames = KMA31GroupsPC, dir = "Likelihood Profiles", filename = "KMA473Pops_31Groups_46loci_LikelihoodProfile_NEW", col = )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Likelihood Profile 17 Groups ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KMA473Pops_17UCIGroups_46loci_LikelihoodProfile <- LeaveOneOutDist.GCL(sillyvec = KMA473Pops, loci = loci46, groupvec = KMA473PopsGroupVec17UCI)
dput(x = KMA473Pops_17UCIGroups_46loci_LikelihoodProfile, file = "Likelihood Profiles/KMA473Pops_17UCIGroups_46loci_LikelihoodProfile.txt")

KMA473Pops17UCIGroups_46loci_Confusion <- ConfusionMatrices.GCL(LeaveOneOutDist = KMA473Pops_17UCIGroups_46loci_LikelihoodProfile, groupnames = KMA17UCIGroups, groupvec = KMA473PopsGroupVec17UCI, sillyvec = KMA473Pops)
dput(x = KMA473Pops17UCIGroups_46loci_Confusion, file = "Objects/KMA473Pops17UCIGroups_46loci_Confusion.txt")

require(lattice)
new.colors <- colorRampPalette(c("white", "black"))
levelplot(KMA473Pops17UCIGroups_46loci_Confusion[[1]], col.regions = new.colors, xlab = "Known Origin", ylab = "Mean Genotype Likelihood", main = "17 Fine Scale Groups with UCI, 46 loci", at = seq(0, 1, length.out = 100), scales = list(x = list(rot = 45)))



KMA473Pops_17UCIGroups_46loci_LikelihoodProfile_NEW <- LeaveOneOutLikeProfile.GCL(popvec = KMA473Pops, loci = loci46, groupvec = KMA473PopsGroupVec17UCI, groupnames = KMA17UCIGroups, groupcomps = NULL, ncores = 8)
dput(x = KMA473Pops_17UCIGroups_46loci_LikelihoodProfile_NEW, file = "Likelihood Profiles/KMA473Pops_17UCIGroups_46loci_LikelihoodProfile_NEW.txt")
PlotLikeProfile.GCL(likeprof = KMA473Pops_17UCIGroups_46loci_LikelihoodProfile_NEW, popvec = KMA473Pops, loci = loci46, groupvec = KMA473PopsGroupVec17UCI, groupnames = KMA17UCIGroups, dir = "Likelihood Profiles", filename = "KMA473Pops_17UCIGroups_46loci_LikelihoodProfile_NEW_", col = c(Colors14[1:11], "blue", "cyan", "pink2", "grey80", Colors14[13:14]))

str(KMA473Pops_17UCIGroups_46loci_LikelihoodProfile_NEW)
OtherUCI <- sapply(KMA473Pops[KMA473PopsGroupVec17UCI == 12], function(silly) {which(KMA473Pops_17UCIGroups_46loci_LikelihoodProfile_NEW$Attributes$FromPop == silly)})

sort(sapply(OtherUCI, function(pop) {mean(KMA473Pops_17UCIGroups_46loci_LikelihoodProfile_NEW$IndividualByGroup[pop, "Other Cook Inlet"])}))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Repeated mixture tests with genetic_msa ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get proof objects from original loci89, 15 RG, 7 scenarios (these have 5 chains each for 5 repeats)
KMA473PopsGroups15RepeatedMixProofTests

KMAobjects <- list.files("MixtureProofTests objects", recursive = FALSE, pattern = ".txt", full.names = FALSE)
KMAobjects

invisible(sapply(KMAobjects, function(objct) {assign(x = unlist(strsplit(x = objct, split = ".txt")), value = dget(file = paste(getwd(), "MixtureProofTests objects", objct, sep = "/")), pos = 1) }))

str(early.Ayakulik1Proof, max.level = 1)


sum(sapply(KMA473Pops, function(silly) {get(paste0(silly, ".gcl"))$n}))

#~~~~~~~~~~~~~~~~~~
# Create baseline and mixture sillyvec
ProofTestIDs <- early.Ayakulik1Proof
sillyvec <- KMA473Pops
loci <- loci89
groupvec <- KMA473PopsGroupVec15
groupnames <- Groups15Short

sapply(ProofTestIDs, function(RG) {sum(sapply(RG, function(pop) {length(pop)}))})  # Confirm mixture proportions


# Code from ReProofTest.GCL
IncludedGroups <- names(ProofTestIDs)
baselinesillyvec <- sillyvec
names(baselinesillyvec) <- sillyvec
mixsillyvec <- NULL
baselineSillysIncluded <- NULL

for(group in IncludedGroups){
  
  g <- match(group,groupnames)
  groupsillys <- sillyvec[groupvec == g]
  newnames <- paste(group, groupsillys, sep = ".")
  names(newnames) <- groupsillys
  names(ProofTestIDs[[group]]) <- groupsillys
  IND <- sapply(groupsillys, function(silly) {!is.null(ProofTestIDs[[group]][[silly]])} )
  newnames <- newnames[IND]
  
  for(silly in groupsillys[IND]){
    newname <- newnames[silly]
    gclname <- paste(silly, "gcl", sep = ".")
    silly.gcl <- get(gclname, pos = 1)
    baselineSillysIncluded <- c(newname, baselineSillysIncluded)
    baselinesillyvec[silly] <- newname
    newgclname <- paste(newname, "gcl", sep = ".")
    newmixname <- paste("mix", newname, sep = ".")
    mixsillyvec <- c(newmixname, mixsillyvec)
    newmixgclname <- paste(newmixname, "gcl", sep = ".")
    
    assign(newgclname, silly.gcl, pos = 1)
    assign(newmixgclname, silly.gcl, pos = 1)
    
    AllIDs <- dimnames(silly.gcl$scores)[[1]]
    baselineIDs2remove <- as.character(sort(as.numeric(ProofTestIDs[[group]][[silly]])))
    mixIDs2remove <- AllIDs[is.na(match(AllIDs,baselineIDs2remove))]
    
    invisible(RemoveIDs.GCL(silly = newname, IDs = list(baselineIDs2remove)))
    invisible(RemoveIDs.GCL(silly = newmixname, IDs = list(mixIDs2remove)))
    
  }  # silly
  
}  # group
#~~~~~~~~~~~~~~~~~~

mixsillyvec
baselinesillyvec

sum(sapply(mixsillyvec, function(silly) {get(paste0(silly, ".gcl"))$n}))
sum(sapply(baselinesillyvec, function(silly) {get(paste0(silly, ".gcl"))$n}))


source("C:/Users/krshedd/Documents/R/GCL-R-Scripts/genetic_msa.R")
# dir.create("genetic_msa")

# Create initial starting values matrix
group_inits <- matrix(data = rep(c(rep(0.3, 3), rep(0.1/12, 15)), 5), nrow = 15, ncol = 5, byrow = FALSE)  # had to define since nchains does not == ngroups

# Pool all mixture sillys into one silly
PoolCollections.GCL(collections = mixsillyvec, loci = loci89, newname = "collection_mix")  # create one mixture silly object (.gcl)
str(collection_mix.gcl)

#~~~~~~~~
collection_mix = "collection_mix"
collections_base = baselinesillyvec
loci = loci89
groupnames = Groups15Short
groups = KMA473PopsGroupVec15
nchains = 5
nits = 4e4
nburn = 2e4
thin = 1
group_prior = rep(1 / max(groups), max(groups))
group_inits = group_inits
out_dir = "genetic_msa"
level = 0.1
#~~~~~~~~

early.Ayakulik1Proof_genetic_msa_Estimates <- genetic_msa(collection_mix = "collection_mix", collections_base = baselinesillyvec, loci = loci89, groupnames = Groups15Short, groups = KMA473PopsGroupVec15, nchains = 5, nits = 4e4, nburn = 2e4, thin = 1, group_prior = rep(1 / max(groups), max(groups)), group_inits = group_inits, out_dir = "genetic_msa", level = 0.1)
# 6.89 hours

str(early.Ayakulik1Proof_genetic_msa_Estimates)
colnames(early.Ayakulik1Proof_genetic_msa_Estimates)[6] <- "GR"

dput(x = early.Ayakulik1Proof_genetic_msa_Estimates, file = "early.Ayakulik1Proof_genetic_msa_Estimates.txt")
early.Ayakulik1Proof_genetic_msa_Estimates <- dget(file = "genetic_msa/early.Ayakulik1Proof_genetic_msa_Estimates.txt")

early.Ayakulik1Proof_BAYES_Estimates <- dget(file = "Estimates objects/loci89/KMA473PopsGroups15loci89RepeatedMixProofTestsEstimatesStats.txt")[["early.Ayakulik1"]]
str(early.Ayakulik1Proof_BAYES_Estimates)
