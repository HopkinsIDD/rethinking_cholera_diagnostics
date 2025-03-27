# This script prepares the data from the seroburden repo for RDT analysis

# Preamble ----------------------------------------------------------------

library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)
library(optparse)
library(here)

source("analysis/utils.R")

# User-supplied options
option_list <- list(
  make_option(c("-v", "--variable"), 
              default = "time_culture", action ="store", type = "character", help = "Variable for which to compute effect"),
  make_option(c("-r", "--redo_data"), 
              default = FALSE, action ="store", type = "logical", help = "redo data processing"),
  make_option(c("-s", "--subset"), 
              default = "all", action ="store", type = "character", help = "subset of participants"),
  make_option(c("-f", "--fake_data"), 
              default = TRUE, action ="store", type = "logical", help = "whether to use fake data or not")
)

opt <- parse_args(OptionParser(option_list = option_list))

print(opt)

cat("---- Processing data for", opt$variable, ".\n")

if (!is.null(opt$variable) & !is.null(opt$equation)) {
  stop()
}

if (opt$subset == "all") {
  opt$subset <- NULL
}

# Load data ---------------------------------------------------------------

cat("---- A. Preparing lab data.\n")

if (opt$fake_data) {
  
  # For testing use fake data that has similar statistical properties as the true data
  lab_data <- readRDS(here("generated_data/fake_lab_data.rds"))
  
} else {
  
  if (opt$redo_data | !file.exists("generated_data/lab_data.rds")) {
    
    lab_data <- read_csv("data/serochit_clinical_2022Aug31.csv") %>% 
      select(date_sample = dttm_sample_collection, 
             date_pcr = pcr_date,
             date_culture = cult_date,
             rdt_res, 
             pcr_res, 
             cult_res, 
             age,
             sex = cat_sex,
             serotype = cat_cul_o1_serotype,
             pcr_delay,
             cult_delay,
             antibiotic_use,
             ind_ab_rx_pre_gtfcc,
             ind_ab_rx_pre_eff,
             cat_dehydration_status,
             duration_of_stay) %>% 
      # Define age categories
      mutate(age_cat = cut(age, getAgeCuts(), right = F)) %>% 
      # Test results to boolean
      # mutate(across(contains("res"), function(x) x == "Positive")) %>% 
      mutate(across(contains("date"), as.Date)) %>% 
      # Times to testing
      rowwise() %>% 
      mutate(time_to_pcr = difftime(date_pcr, date_sample, units = "days") %>% as.numeric(),
             time_to_culture = difftime(date_culture, date_sample, units = "days") %>% as.numeric()) %>% 
      mutate(rdt_period =  map_chr(date_sample, ~ defineRDTPeriod(.)))  %>% 
      mutate(time_to_pcr_cat = cut(time_to_pcr, c(seq(0, 30, by = 5), Inf)),
             time_to_culture_cat = cut(time_to_culture, c(seq(0, 30, by = 5), Inf))) %>% 
      ungroup() %>% 
      # Epidemiological period
      mutate(epi_period = map_chr(date_sample, ~ defineEpiPeriod(.))) %>% 
      mutate(rdt_period = forcats::fct_relevel(rdt_period, c("Post-batch")),
             epi_period = forcats::fct_relevel(epi_period, c("monsoon")),
             age_cat = forcats::fct_relevel(age_cat, c("[5,Inf)"))) %>% 
      mutate(cat_dehydration_status = case_when(cat_dehydration_status == "Don't know" ~ "No signs",
                                                T ~ cat_dehydration_status) %>% 
               as.factor() %>% 
               forcats::fct_relevel("Some")) 
    
    # Add climate data
    chittaghon_weather <- readRDS("generated_data/chittaghon_airport_weather_data_filled.rds")
    
    lab_data <- lab_data %>% 
      left_join(chittaghon_weather %>% 
                  # Lag by incubation period
                  select(date, temp = temp_pred),
                by = c("date_sample" = "date"))
    
    
    lab_data <- lab_data %>% 
      mutate(time_to_pcr_std = scale(time_to_pcr),
             time_to_culture_std = scale(time_to_culture),
             temp_std = scale(temp)) %>% 
      # Temp std by epi period
      group_by(epi_period) %>% 
      mutate(temp_std_period = scale(temp)) %>% 
      ungroup() %>% 
      # Fill time to pcr for data processing
      replace_na(list(time_to_pcr_std = 0,
                      time_to_culture_std = 0))
    
    # Set the week
    lab_data <- lab_data %>% 
      mutate(week = map_chr(date_sample, ~str_c(lubridate::epiyear(.), "_", lubridate::epiweek(.))))
    
    if (file.exists("data/incidence_inference_data.rdata")) {
      load("data/incidence_inference_data.rdata")
    } else {
      # Get the incidence and prob values
      prob_traj <- readRDS("../archive/cholera_seroburden/generated_data/est_prob_traj_survonly_rdt_analysis.rds")
      incid_samples <- readRDS("../archive/cholera_seroburden/generated_data/incid_prior_survonly_stan_output_rdt_analysis.rds")
      
      # TODO hack to correct age categories
      hackAgeCats <- function(df) {
        df %>% 
          filter(age_cat != "[65,Inf)") %>% 
          mutate(age_cat = case_when(age_cat == "[0,5)" ~ age_cat,
                                     TRUE ~ "[5,Inf)"))
      }
      
      chol_traj <- incid_samples$summary("I_chol", .cores = 6) %>% 
        bind_cols(prob_traj %>% select(date, age_cat)) %>% 
        hackAgeCats()
      
      awd_traj <- incid_samples$summary("I_nonchol", .cores = 6) %>% 
        bind_cols(prob_traj %>% select(date, age_cat)) %>% 
        hackAgeCats()
      
      prob_traj <- prob_traj %>% hackAgeCats()
      pop <- getSitakundatData(data_path = "../archive/cholera_seroburden/serochit_data/census_data_sex_age_sitakunda.csv",
                               age_categories = c("[0,5)", "[5,Inf)"))
      
      save(prob_traj, chol_traj, awd_traj, pop, file = "data/incidence_inference_data.rdata")
    }
    
    # Find prob of cholera for each sample
    lab_data <- lab_data %>% 
      mutate(prior_chol = NA,
             chol_incid = NA,
             awd_incid = NA)
    
    for (i in 1:nrow(lab_data)) {
      p <- prob_traj %>% 
        filter(age_cat == lab_data$age_cat[i],
               date == lab_data$date_sample[i]) %>% 
        pull(mean)
      
      if (length(p) > 0) {
        lab_data$prior_chol[i] <- p
      } else {
        lab_data$prior_chol[i] <- .1
      }
      
      chol <- chol_traj %>% 
        filter(age_cat == lab_data$age_cat[i],
               date == lab_data$date_sample[i]) %>% 
        pull(mean)
      
      nonchol <- awd_traj %>% 
        filter(age_cat == lab_data$age_cat[i],
               date == lab_data$date_sample[i]) %>% 
        pull(mean)
      
      if (length(chol) == 0) {
        lab_data$chol[i] <- mean(chol_traj$mean[chol_traj$age_cat == lab_data$age_cat[i]])/pop$pop[pop$age_cat == lab_data$age_cat[i]]
        lab_data$nonchol[i] <- mean(awd_traj$mean[awd_traj$age_cat == lab_data$age_cat[i]])/pop$pop[pop$age_cat == lab_data$age_cat[i]]
      } else {
        lab_data$chol[i] <- chol/pop$pop[pop$age_cat == lab_data$age_cat[i]]
        lab_data$nonchol[i] <- nonchol/pop$pop[pop$age_cat == lab_data$age_cat[i]]
      }
    }
    
    
    lab_data <- lab_data %>% 
      mutate(
        chol_std = scale(chol),
        nonchol_std = scale(nonchol)
      )
    
    lab_data <- lab_data %>% 
      mutate(month = map_chr(date_sample, ~ str_c(lubridate::year(.), lubridate::month(.), sep = "_")))
    
    saveRDS(lab_data, "generated_data/lab_data.rds")
    
  } else {
    lab_data <- readRDS("generated_data/lab_data.rds")
  }
}

lab_data <- lab_data %>% 
  # !! remove data not consistent with protocol
  filter(!((is.na(pcr_res) | is.na(cult_res)) & rdt_res))

# Prepare data for stan ---------------------------------------------------
cat("---- B. Preparing covariate data.\n")

if (!is.null(opt$subset)) {
  cat("Subseting data for: ", opt$subset, "\n")
  if (opt$subset == "children") {
    lab_data <- lab_data %>% 
      filter(age_cat == "[0,5)")
  } else if (opt$subset == "adults") {
    lab_data <- lab_data %>% 
      filter(age_cat == "[5,Inf)")
  } else if (opt$subset != "all") {
    stop("Did not find option to subset")
  }
}

# Change antibiotics variable if asked for
if (opt$variable == "antibioticsGTFCC") {
  lab_data <- lab_data %>% 
    mutate(antibiotic_use = case_when(ind_ab_rx_pre_gtfcc == "0" ~ 0,
                                      TRUE ~ 1))
  variable <- "antibiotics"
} else if (opt$variable == "antibioticsEff") {
  lab_data <- lab_data %>% 
    mutate(antibiotic_use = case_when(ind_ab_rx_pre_eff == "0" ~ 0,
                                      TRUE ~ 1))
  variable <- "antibiotics"
} else if (opt$variable != "antibiotics") {
  lab_data <- lab_data %>% 
    mutate(antibiotic_use = case_when(ind_ab_rx_pre_gtfcc == "0" ~ 0,
                                      TRUE ~ 1))
  variable <- opt$variable
} else {
  variable <- opt$variable
}

# # !!!! 
# lab_data <- lab_data %>% 
#   mutate(cult_res = case_when(!rdt_res & pcr_res ~ TRUE,
#                               !rdt_res & !pcr_res ~ FALSE,
#                               T ~ cult_res))


# Matrix of N x 3 of results of three tests [rdt, pcr, culture]
test_res_mat <- lab_data %>% 
  select(rdt_res, pcr_res, cult_res) %>% 
  as.matrix()

# Test results available for maringalization:
# 1: all three available
# 2: RDT and PCR
# 3: RDT only

test_cat <- apply(test_res_mat, 1, function(x) {
  n_na <- sum(is.na(x))
  case_when(n_na == 1 ~ 2,
            n_na == 2 ~ 3,
            T ~ 1)
})


# Fill NAs with FALSE
test_res_mat[is.na(test_res_mat)] <- FALSE

# Regression matrices for cholera probability
if (is.null(opt$subset)) {
  covar_mat_chol <- model.matrix(data = lab_data,
                                 ~ temp_std + age_cat + cat_dehydration_status - 1)
} else {
  covar_mat_chol <- model.matrix(data = lab_data,
                                 ~ temp_std + cat_dehydration_status - 1)
}

# !! Change the temperature to be period-scaled
lab_data <- rename(lab_data, 
                   temp_std_all = temp_std,
                   temp_std = temp_std_period)

# Regression matrices for each test
eqs <- makeEquations(variable = variable,
                     subset = opt$subset)

covar_mats <- makeCovarMats(lab_data = lab_data,
                            eqs = eqs)

if (opt$redo_data | !file.exists("generated_data/all_u_covar.rds")) {
  # Unique combinations based on dates in data
  all_u_covar <- expand.grid(
    age_cat = unique(lab_data$age_cat),
    cat_dehydration_status = unique(lab_data$cat_dehydration_status),
    # date = sort(unique(lab_data$date_sample)),
    time_to_culture_std = (c(0, 3, 7, 10, 14, 17, 21, 25) - mean(lab_data$time_to_culture, na.rm = T)) / sd(lab_data$time_to_culture, na.rm = T),
    temp_std = seq(-2, 2, by = .25),
    antibiotic_use = unique(lab_data$antibiotic_use),
    rdt_period = unique(lab_data$rdt_period),
    epi_period = unique(lab_data$epi_period)
  ) %>% 
    as_tibble()
  
  # Set references
  all_u_covar <- all_u_covar %>% 
    mutate(rdt_period = forcats::fct_relevel(rdt_period, c("Post-batch")),
           epi_period = forcats::fct_relevel(epi_period, c("monsoon")),
           age_cat = forcats::fct_relevel(age_cat, c("[5,Inf)")))  
  
  saveRDS(all_u_covar, file = "generated_data/all_u_covar.rds")
  
} else {
  all_u_covar <- readRDS("generated_data/all_u_covar.rds")
}

u_covar_dfs <- makeUCovarDFs(all_u_covar = all_u_covar,
                             eqs = eqs)

u_covar_mats <- makeUCovarMats(u_covar_dfs = u_covar_dfs,
                               eqs = eqs,
                               variable = variable)

# Save covariates
save(eqs, covar_mats, u_covar_dfs, u_covar_mats, 
     file = makeCovarDataFilename(opt = opt))

# Post-stratification matrices --------------------------------------------

# Define number of unique post-stratified strata
poststrat_matrices <- getAllPostStratMatrices(variable = variable,
                                              u_covar_mats = u_covar_mats,
                                              u_covar_dfs = u_covar_dfs,
                                              eqs = eqs,
                                              lab_data = lab_data)

U_poststrat <- max(map_dbl(poststrat_matrices, nrow))

# Fixt post-strat matrices for empty data
poststrat_matrices <- map(poststrat_matrices, function(x) {
  if(nrow(x) != U_poststrat) {
    matrix(1, nrow = U_poststrat, ncol = 1)
  } else {
    x
  }
})

# Prior infection prob ----------------------------------------------------

cat("---- C. Preparing stan data.\n")

# Test performance priors
sens_bounds <- tribble(
  ~test, ~mean, ~lo, ~hi,
  # "RDT", 0.98, 0.88, 0.99,
  "RDT", 0.7, 0.5, 0.99,
  "PCR", 0.7, 0.59, 0.9,
  "cul", 0.71, 0.59, 0.81
) %>%
  mutate(across(c("mean", "lo", "hi"), ~ logit(. - 1e-6))) %>%
  mutate(sd = (mean - lo)/2)

## tightening the culture specificity 
spec_bounds <- tribble(
  ~test, ~mean, ~lo, ~hi,
  # "RDT", 0.97, 0.89, 1,
  "RDT", 0.9, 0.6, 1,
  "PCR", 0.97, 0.8, 1,
  # "PCR", 0.97, 0.6, 1,
  "cul", 0.999, 0.997, 1
)  %>%
  mutate(across(c("mean", "lo", "hi"), ~ logit(. - 1e-6))) %>%
  mutate(sd = (mean - lo)/2)

logit_sens_prior_mu <- sens_bounds$mean
logit_sens_prior_sd <- sens_bounds$sd

# Set prior on first period of RDT batch
logit_spec_prior_mu <- spec_bounds$mean
logit_spec_prior_sd <- spec_bounds$sd
logit_spec_prior_sd[3] <- .05


# Fraction of RDT- PCR-tested ---------------------------------------------

N_month <- length(unique(lab_data$month))
frac_B <- computeFracB(lab_data)

map_month <- map_dbl(lab_data$month, ~ which(frac_B$month == .))

# Stan data ---------------------------------------------------------------

u_weeks <- lab_data %>%
  group_by(week) %>% 
  slice_min(date_sample, with_ties = F) %>% 
  arrange(date_sample) %>% 
  pull(week)

stan_data <- list(
  N = nrow(test_res_mat),    # Number of observations
  K = 3,                     # Number of test configurations
  M_chol = ncol(covar_mat_chol),    # Number of covariates for cholera infection
  y = test_res_mat,
  test_cat = test_cat,
  prior_chol = lab_data$prior_chol,    # Prior probability of being positive based on seroburden model
  X_chol = covar_mat_chol,
  logit_sens_prior_mu = logit_sens_prior_mu,
  logit_sens_prior_sd = logit_sens_prior_sd,
  logit_spec_prior_mu = logit_spec_prior_mu,
  logit_spec_prior_sd = logit_spec_prior_sd,
  sd_beta = 1,
  N_week = length(u_weeks),
  map_week = map_dbl(lab_data$week, ~ which(u_weeks == .)),
  N_month = N_month,
  frac_B = frac_B$frac_B,
  map_month = map_month
)

# Append all covariate objects
covar_stan_data <- map(
  names(eqs), 
  function(x) {
    list(
      M = ncol(covar_mats[[x]]),
      X = covar_mats[[x]],
      U = nrow(u_covar_mats[[x]]),
      X_u = u_covar_mats[[x]]) %>% 
      magrittr::set_names(str_c(names(.), x, sep = "_"))
  }) %>% 
  unlist(recursive = F)

stan_data <- append(stan_data, covar_stan_data)

# Post-stratification matrices
poststrat_stan_data <- map(
  names(poststrat_matrices), 
  function(x) {
    list(X_poststrat = poststrat_matrices[[x]]) %>% 
      magrittr::set_names(str_c(names(.), x, sep = "_"))
  }) %>% 
  unlist(recursive = F) %>% 
  append(list(U_poststrat = U_poststrat))

stan_data <- append(stan_data, poststrat_stan_data)


if(opt$variable == "age_timevary") {
  stan_data$autocorr_order <- 1 
}

saveRDS(stan_data, file = makeStanDataFilename(opt = opt))
