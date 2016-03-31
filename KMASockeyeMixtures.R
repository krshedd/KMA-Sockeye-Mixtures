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
#### Directory Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Create file structure based on OLD
# setwd("V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline")
# dirsOLD <- list.dirs(path = "V:/WORK/Sockeye/Kodiak/KMA 2014-2016/Baseline/OLD", recursive = FALSE, full.names = FALSE)
# sapply(dirsOLD, dir.create)
# sapply(c("Baseline", "BAYES.baseline", "BAYES.control", "BAYES.mixture", "BAYES.output"), function(folder) {dir.create(path = paste(getwd(), "BAYES", folder, sep = "/"))})


rm(list = ls(all = TRUE))
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Mixtures")
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")