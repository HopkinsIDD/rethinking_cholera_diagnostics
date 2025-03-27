# This script runs simulations to show impact of test performance on RDT evaluation

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)

source("analysis/utils.R")

set.seed(12302)

# Functions ----------------------------------------------------------------

computeSimPerf <- function(sim_vec, 
                           ref = "culture") {
  pos <- sum(sim_vec[1:4])
  neg <- sum(sim_vec[5:8])
  
  if (ref == "culture") {
    true_pos <- sum(sim_vec[c(1, 3)])
    true_neg <- sum(sim_vec[c(6, 8)])
    false_pos <- sum(sim_vec[c(2, 4)])
    false_neg <- sum(sim_vec[c(5, 7)])
  } else if (ref == "pcr") {
    true_pos <- sum(sim_vec[c(1, 2)])
    true_neg <- sum(sim_vec[c(7, 8)])
    false_pos <- sum(sim_vec[c(3, 4)])
    false_neg <- sum(sim_vec[c(5, 6)])
  } else if (ref == "culture_pcr") {
    # Combined reference
    # Page, A.L., Alberti, K.P., Mondonge, V., Rauzier, J., Quilici, M.L. and Guerin, P.J., 2012. Evaluation of a rapid test for the diagnosis of cholera in the absence of a gold standard. PLoS One, 7(5), p.e37360.
    #   "Culture and PCR reference standard: a sample was considered positive if 
    #   any of the culture or PCR results were positive for the detection of 
    #   V. cholerae O1 or O139. A sample was considered negative if both culture 
    #   results were available and both negative, and PCR was also negative. 
    #   As above, specimens with only one negative culture result available were
    #   considered indeterminate and excluded from the analysis."
    true_pos <- sum(sim_vec[c(1, 2, 3)])
    true_neg <- sum(sim_vec[c(8)])
    false_pos <- sum(sim_vec[c(4)])
    false_neg <- sum(sim_vec[c(5, 6, 7)])
      
  } else {
    stop("Reference unknown: ", ref)
  }
  
  # Compute stats
  sens <- true_pos/(true_pos + false_neg)
  spec <- true_neg/(true_neg + false_pos)
  
  if (is.nan(spec)) {
    spec <- 1
  }
  
  if (is.nan(sens)) {
    sens <- 1
  }
  
  tibble(
    sens = sens,
    spec = spec
  )
}

computeProbs <- function(chol_prior,
                         sens,
                         spec) {
  
  # 1. RDT+ PCR+ CUL+
  # 2. RTD+ PCR+ CUL-
  # 3. RDT+ PCR- CUL+
  # 4. RDT+ PCR- CUL-
  # 5. RDT- PCR+ CUL+
  # 6. RDT- PCR+ CUL-
  # 7. RDT- PCR- CUL+
  # 8. RDT- PCR- CUL-
  
  p <- rep(NA, 8)
  x <- chol_prior
  
  p[1] = prod(sens) * x + prod(1-spec) * (1-x);
  p[2] = sens[1]*sens[2]*(1-sens[3]) * x + (1-spec[1])*(1-spec[2])*spec[3] * (1-x);
  p[3] = sens[1]*(1-sens[2])*sens[3] * x + (1-spec[1])*spec[2]*(1-spec[3]) * (1-x);
  p[4] = sens[1]*(1-sens[2])*(1-sens[3]) * x + (1-spec[1])*spec[2]*spec[3] * (1-x);
  p[5] = (1-sens[1])*sens[2]*sens[3] * x + spec[1]*(1-spec[2])*(1-spec[3]) * (1-x);
  p[6] = (1-sens[1])*sens[2]*(1-sens[3]) * x + spec[1]*(1-spec[2])*spec[3] * (1-x);
  p[7] = (1-sens[1])*(1-sens[2])*sens[3] * x + spec[1]*spec[2]*(1-spec[3]) * (1-x);
  p[8] = prod(1-sens) * x + prod(spec) * (1-x);
  
  p
}

simulateTests <- function(n_sim = 1,
                          n_samples,
                          chol_prior,
                          sens,
                          spec) {
  
  # Compute priors
  probs <- computeProbs(chol_prior = chol_prior,
                        sens = sens,
                        spec = spec)
  
  # Simulate multinomial
  rmultinom(n_sim, n_samples, prob = probs)
}

unpackPerf <- function(vec) {
  sens <- vec[c("rdt_sens", "pcr_sens", "culture_sens")]
  spec <- vec[c("rdt_spec", "pcr_spec", "culture_spec")]
  
  list(
    sens = sens,
    spec = spec
  )
}

# Load samples ------------------------------------------------------------

var <- "mean"

# Samples
genquant <- readRDS(str_glue("generated_data/genquant_{var}.rds"))

# Covariates
load(str_glue("generated_data/covar_data_{var}.rds"))

# Make samples
gen_perf <- str_c("gen_", getAllTests())

test_performance <- genquant$summary(gen_perf, .cores = 4)

test_performance <- test_performance %>%
  mutate(test = str_extract(variable, "rdt|pcr|culture") %>% factor(levels = c("rdt", "pcr", "culture")),
         what = str_extract(variable, "sens|spec"),
         id = str_extract(variable, "[0-9]+") %>% as.numeric()) %>%
  rowwise() %>%
  mutate(
    covar = colnames(u_covar_mats[[str_c(test, what, sep = "_")]])[max(ncol(u_covar_mats[[str_c(test, what, sep = "_")]]), 1)],
    covar_val = u_covar_mats[[str_c(test, what, sep = "_")]][id, max(ncol(u_covar_mats[[str_c(test, what, sep = "_")]]), 1)]
  ) %>% 
  ungroup() %>% 
  group_by(test, what) %>% 
  mutate(fixed = max(id) == 1) %>% 
  ungroup()

saveRDS(test_performance, file = str_glue("generated_data/generated_test_performance_{var}.rds"))


# Simulate ----------------------------------------------------------------

# Define number of samples per simulation
n_samples <- 300

# Define range of true cholera prevalence to cover
chol_priors <- seq(0.01, .99, by = .05)

sim_perf <- 
  map_df(
    c("pcr", "culture", "culture_pcr"), 
    function(z) {
      map_df(
        chol_priors, 
        function(y) {
          map_df(
            unique(test_performance$id), 
            function(x) {
              df <- filter(test_performance, id == x | fixed)
              
              # Make vector of test performances
              test_vec <- df$mean
              names(test_vec) <- df$variable %>% str_extract(str_c(getAllTests(), collapse = "|"))
              
              sens_spec <- unpackPerf(vec = test_vec)
              
              # Simulate samples
              test_sim <- simulateTests(n_sim = 1000, 
                                        n_samples = n_samples, 
                                        chol_prior = y, 
                                        sens = sens_spec$sens, 
                                        spec = sens_spec$spec)
              
              # Compute RDT performance based on reference
              rdt_perf <- apply(test_sim, 2, computeSimPerf, ref = z) %>% 
                reduce(bind_rows) %>% 
                pivot_longer(cols = everything()) %>% 
                group_by(name) %>% 
                summarise(mean = mean(value),
                          q5 = quantile(value, .05, na.rm = T),
                          q95 = quantile(value, .95, na.rm = T))
              
              rdt_perf %>% 
                mutate(chol_prior = y,
                       id = x,
                       ref = z)
            })
        })
    })

saveRDS(sim_perf, file = str_glue("generated_data/simulated_performance_{var}.rds"))

