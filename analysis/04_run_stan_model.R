# This script runs the stan model

# Preamble ----------------------------------------------------------------

library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)
library(cmdstanr)
library(optparse)

source("analysis/utils.R")

# User-supplied options
option_list <- list(
  make_option(c("-v", "--variable"), 
              default = "time_culture", action ="store", type = "character", help = "Variable for which to compute effect"),
  make_option(c("-r", "--redo_stanfit"), 
              default = FALSE, action ="store", type = "logical", help = "Redo stan fitting"),
  make_option(c("-s", "--subset"), 
              default = "all", action ="store", type = "character", help = "subset of participants"),
  make_option(c("-f", "--fake_data"), 
              default = TRUE, action ="store", type = "logical", help = "whether to use fake data or not")
)

opt <- parse_args(OptionParser(option_list = option_list))

print(opt)

if (opt$subset == "all") {
  opt$subset <- NULL
}

# Load stan data ----------------------------------------------------------

stan_data <- readRDS(file = makeStanDataFilename(opt = opt))

# Compile model -----------------------------------------------------------

lc_model <- cmdstanr::cmdstan_model("analysis/stan/latent_class_model.stan")
lc_model_gen <- cmdstanr::cmdstan_model("analysis/stan/latent_class_genquant.stan")

# Run sampling ------------------------------------------------------------

if (!file.exists(makeStanOutputFilename(opt = opt)) | opt$redo_stanfit) {
  
  lc_samples <- lc_model$sample(
    data = stan_data,
    iter_warmup = 250,
    iter_sampling = 1250,
    parallel_chains = 4,
    init = .1,
    max_treedepth = 11,
    refresh = 100,
    save_warmup = FALSE
  )
  
  # Save
  lc_samples$save_object(makeStanOutputFilename(opt = opt))
  
  # Diagnostics
  lc_samples$cmdstan_diagnose()
  
} else {
  cat("-- Loading pre-computed draws \n")
  lc_samples <- readRDS(makeStanOutputFilename(opt = opt))
}

# # # Diagnostics
# lc_samples$draws("beta_rdt_sens") %>%
#   bayesplot::mcmc_intervals()
# 
# lc_samples$draws("beta_chol") %>%
#   bayesplot::mcmc_intervals()


# Generated quantities ----------------------------------------------------

genquant <- lc_model_gen$generate_quantities(lc_samples, 
                                             data = stan_data,
                                             parallel_chains = 4)

genquant$save_object(makeGenquantFilename(opt = opt))


