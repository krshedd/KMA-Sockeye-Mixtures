# The goal of this script is to test the KMA sockeye baseline with Eric Anderson's `rubias` package
# 1) Look at self-assignment tests
# 2) Re-analyze the mixture scenario proof tests with `rubias` and compare to `BAYES`

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Convert Baseline to `rubias` ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/KMA Commercial Harvest 2014-2016/Baseline/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get Locus Control, objects, and populations
LocusControl <- dget(file = "Objects/LocusControl103.txt")
loci89 <- dget(file = "Objects/loci89.txt")
KMA473Pops <- dget(file = "Objects/KMA473Pops.txt")
KMA473PopsGroupVec15 <- dget(file = "Objects/KMA473PopsGroupVec15.txt")
Groups15 <- dget(file = "Objects/Groups15.txt")

require(beepr)
invisible(sapply(KMA473Pops, function(silly) {assign(x = paste(silly, ".gcl", sep = ""), value = dget(file = paste0(getwd(), "/Raw genotypes/PostQCPopsloci103/", silly, ".txt")), pos = 1)} )); beep(2)
length(objects(pattern = "\\.gcl"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load gcl functions
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")
source("C:/Users/krshedd/Documents/R/GCL-R-Scripts/create_rubias_baseline.R")
kma_473pops_89loci_15groups.rubias_base <- create_rubias_baseline(sillyvec = KMA473Pops, loci = loci89, group_names = Groups15, groupvec = KMA473PopsGroupVec15)
str(kma_473pops_89loci_15groups.rubias_base, max.level = 2)

require(rubias)
require(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Baseline testing with `rubias` ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Self assignment
kma_473pops_89loci_15groups.rubias_base_sa <- self_assign(reference = kma_473pops_89loci_15groups.rubias_base, gen_start_col = 5)
str(kma_473pops_89loci_15groups.rubias_base_sa, max.level = 1)

sa_to_repu <- kma_473pops_89loci_15groups.rubias_base_sa %>% 
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))

sa_to_repu %>% 
  filter(inferred_repunit == "Black Lake") %>% 
  filter(repunit %in% c("Black Lake", "Chignik Lake")) %>% 
  ggplot(aes(x = collection, y = repu_scaled_like)) +
  geom_boxplot()

sa_to_repu %>% 
  filter(repunit == "Black Lake") %>% 
  ggplot(aes(x = inferred_repunit, y = repu_scaled_like)) +
  geom_boxplot()

kma_473pops_89loci_15groups.rubias_base_sa %>% 
  filter(repunit == "Black Lake") %>% 
  filter(inferred_repunit %in% c("Black Lake", "Chignik Lake")) %>% 
  ggplot(aes(x = inferred_collection, y = scaled_likelihood)) +
  geom_boxplot(aes(colour = inferred_repunit))

# Confusion matrix
str(sa_to_repu)
sa_to_repu %>% 
  mutate(repunit_f = factor(x = repunit, levels = Groups15)) %>% 
  mutate(inferred_repunit_f = factor(x = inferred_repunit, levels = Groups15)) %>% 
  group_by(collection, repunit_f, inferred_repunit_f) %>% 
  summarise(mean_repu_scaled_like_col = mean(repu_scaled_like)) %>% 
  ungroup() %>% 
  group_by(repunit_f, inferred_repunit_f) %>% 
  summarise(mean_repu_scaled_like = mean(mean_repu_scaled_like_col)) %>% 
  # spread(repunit, mean_repu_scaled_like)
  ggplot(aes(x = repunit_f, y = inferred_repunit_f, z = mean_repu_scaled_like)) +
  geom_tile(aes(fill = mean_repu_scaled_like)) +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Leave-one-out
kma_473pops_89loci_15groups.rubias_base_loo <- assess_reference_loo(reference = kma_473pops_89loci_15groups.rubias_base, gen_start_col = 5, reps = 50, mixsize = 200)
kma_473pops_89loci_15groups.rubias_base_loo

tmp <- kma_473pops_89loci_15groups.rubias_base_loo %>% 
  group_by(iter, repunit) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))
tmp

ggplot(tmp, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)

#~~~~~~~~~~~~~~~~~~
## Leave-one-out scenario

# Grab scenarios
fishery_scenarios <- read.table(file = "FisheryProofTestScenarios.txt", header = TRUE, sep = "\t")
arep <- sapply(colnames(fishery_scenarios)[-1], function(scn) {data.frame(repunit = Groups15, ppn = fishery_scenarios[, scn])}, simplify = FALSE)

# Test Scenarios
kma_473pops_89loci_15groups.rubias_base_loo_scn <- assess_reference_loo(reference = kma_473pops_89loci_15groups.rubias_base, gen_start_col = 5, reps = 50, mixsize = 200, alpha_repunit = arep)
dput(x = kma_473pops_89loci_15groups.rubias_base_loo_scn, file = "Objects/kma_473pops_89loci_15groups.rubias_base_loo_scn.txt")

# Summarize data at repunit level
tmp <- kma_473pops_89loci_15groups.rubias_base_loo_scn %>% 
  mutate(repunit_f = factor(x = repunit, levels = Groups15)) %>% 
  group_by(repunit_scenario, iter, repunit_f) %>% 
  summarise(true_repprop = sum(true_pi), repprop_posterior_mean = sum(post_mean_pi), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))
tmp

fishery_scenarios.chr <- colnames(fishery_scenarios)[-1]

# Plot all 7 scenarios on separate plots
sapply(fishery_scenarios.chr, function(scenario) {
  tmp.scn <- filter(.data = tmp, repunit_scenario == scenario)
  dev.new()
  print(
    ggplot(tmp.scn, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit_f)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1) +
      facet_wrap(~ repunit_f)
  )
})


ggplot(tmp, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit_f)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## MC
kma_473pops_89loci_15groups.rubias_base_mc <- assess_reference_mc(reference = kma_473pops_89loci_15groups.rubias_base, gen_start_col = 5, reps = 50, mixsize = 200)
kma_473pops_89loci_15groups.rubias_base_mc

tmp <- kma_473pops_89loci_15groups.rubias_base_mc %>% 
  group_by(iter, repunit) %>% 
  summarise(true_repprop = sum(omega), repprop_posterior_mean = sum(post_mean), repu_n = sum(n)) %>% 
  mutate(repu_n_prop = repu_n / sum(repu_n))
tmp

ggplot(tmp, aes(x = repu_n_prop, y = repprop_posterior_mean, colour = repunit)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~ repunit)