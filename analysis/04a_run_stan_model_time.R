# This script runs the stan model


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(optparse)

source("analysis/utils.R")

# User-supplied options
option_list <- list(
  make_option(c("-v", "--variable"), 
              default = "age_timevary", action ="store", type = "character", help = "Variable for which to compute effect")
)

opt <- parse_args(OptionParser(option_list = option_list))


# Load stan data ----------------------------------------------------------

stan_data <- readRDS(file = makeStanDataFilename(variable = opt$variable))

# Compile model -----------------------------------------------------------

lc_model <- cmdstanr::cmdstan_model("analysis/stan/latent_class_model_time.stan")
# lc_model_gen <- cmdstanr::cmdstan_model("analysis/stan/latent_class_genquant.stan")

# Run sampling ------------------------------------------------------------

lc_samples <- lc_model$sample(
  data = stan_data,
  iter_warmup = 50,
  iter_sampling = 50,
  parallel_chains = 4,
  init = .1,
  max_treedepth = 10,
  refresh = 10,
  save_warmup = TRUE,
  output_dir = "generated_data/",
  output_basename = "time_vary_4"
)


# Save
lc_samples$save_object(makeStanOutputFilename(variable = opt$variable))

# Diagnostics
lc_samples$cmdstan_diagnose()

# # Diagnostics
# lc_samples$draws("beta_sens_rdt") %>% 
#   bayesplot::mcmc_intervals()

# lc_samples$draws(c("gen_sens_rdt",
#                    "gen_spec_rdt",
#                    "gen_sens_pcr",
#                    "gen_spec_pcr",
#                    "gen_sens_culture")) %>% 
#   bayesplot::mcmc_intervals()
# 
# lc_samples$summary(c("gen_sens_rdt",
#                      "gen_spec_rdt",
#                      "gen_sens_pcr",
#                      "gen_spec_pcr",
#                      "gen_sens_culture"))

# Generated quantities ----------------------------------------------------

genquant <- lc_model_gen$generate_quantities(lc_samples, 
                                             data = stan_data,
                                             parallel_chains = 4)

genquant$save_object(makeGenquantFilename(variable = opt$variable))

# 
# get_prop <- function(v) {
#   counts <- table(v)
#   
#   # Complete with all cases
#   u_cats <- 1:8
#   all_counts <- rep(0, length(u_cats))
#   names(all_counts) <- u_cats
#   all_counts[names(counts)] <- counts
#   
#   # Compute cumulative probability of being in a risk category larger or equal
#   res <- all_counts/sum(all_counts)
#   
#   res
# }
# 
# cat_post <- genquant$summary("gen_cat", get_prop, .cores = 3)
# 
# 
# lab_data <- lab_data %>% 
#   addEpiWeek(date_col = "date_sample") %>% 
#   mutate(obs_id = row_number())
# 
# cat_long <- cat_post %>% 
#   pivot_longer(cols = as.character(1:8),
#                names_to = "category",
#                values_to = "prop") %>% 
#   mutate(catnum = as.numeric(category),
#          obs_id = str_extract(variable, "[0-9]+") %>% as.numeric()) %>% 
#   inner_join(lab_data %>% select(set) %>% mutate(obs_id = row_number()))
# 
# prob_A <- cat_long %>% 
#   filter(set == "A") %>% 
#   group_by(obs_id) %>% 
#   summarise(prop_A = sum(prop[catnum >4]))
# 
# prob_A %>% 
#   ggplot(aes(x = prop_A)) +
#   geom_histogram()
# 
# prob_B <- cat_long %>% 
#   filter(set == "B") %>% 
#   group_by(obs_id) %>% 
#   summarise(prop_B1 = sum(prop[catnum %in% c(5, 6)]),
#             prop_B2 = sum(prop[catnum %in% c(7, 8)]),
#             prop_nonB = 1 - prop_B1 - prop_B2)
# 
# prob_C <- cat_long %>% 
#   filter(set == "C") %>% 
#   group_by(obs_id) %>% 
#   summarise(prop_C1 = sum(prop[catnum %in% c(1)]),
#             prop_C2 = sum(prop[catnum %in% c(2)]),
#             prop_C3 = sum(prop[catnum %in% c(3)]),
#             prop_C4 = sum(prop[catnum %in% c(4)]),
#             # Proportion that should have had a neg RDT test
#             prop_nonC = 1 - sum(prop[catnum %in% 1:4]))
# 
# prob_C %>% 
#   ggplot(aes(x = prop_nonC)) +
#   geom_histogram()
# 
# nonC <- prob_C %>% filter(prop_nonC > .5) %>% pull(obs_id)
# 
# lab_data_C %>% 
#   inner_join(prob_C) %>% 
#   ggplot(aes(x = prop_nonC, fill = diff)) +
#   geom_histogram()
# 
# issues <- lab_data_C %>% 
#   inner_join(prob_C) %>% 
#   # filter(obs_id %in% nonC) %>% 
#   mutate(p_chol = prob_chol$mean[lab_data_C$obs_id])
# 
# ggplot(issues, aes(x = p_chol, y = prop_nonC, color = diff)) +
#   geom_point() +
#   facet_wrap(~age_cat)
# 
# pchol <- lc_samples$summary("p_chol", .cores = 4)
# 
# lab_data %>% 
#   bind_cols(pchol) %>% 
#   filter(!rdt_res, pcr_res) %>% 
#   ggplot(aes(mean)) +
#   geom_histogram()
