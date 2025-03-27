// Latent class model for RDT performance analysis

data {
  int<lower=0> N;            // Number of observations
  int<lower=1> M_rdt_sens;        // Number of covariates RDT performance
  int<lower=1> M_rdt_spec;        // Number of covariates RDT performance
  int<lower=1> M_pcr_sens;        // Number of covariates PCR performance
  int<lower=1> M_pcr_spec;        // Number of covariates PCR performance
  int<lower=1> M_culture_sens;    // Number of covariates Culture performance
  int<lower=1> M_culture_spec;    // Number of covariates Culture performance
  int<lower=1> M_chol;       // Number of covariates for the probability of having cholera
  int<lower=1> U_rdt_sens;        // Number of unique covariate combinations for RDT performance
  int<lower=1> U_rdt_spec;        // Number of unique covariate combinations for RDT performance
  int<lower=1> U_pcr_sens;        // Number of unique covariate combinations for PCR performance
  int<lower=1> U_pcr_spec;        // Number of unique covariate combinations for PCR performance
  int<lower=1> U_culture_sens;    // Number of unique covariate combinations for Culture performance
  int<lower=1> U_culture_spec;    // Number of unique covariate combinations for Culture performance
  int<lower=1> U_poststrat;       // Number of post-stratified groups
  
  array[N, 3] int<lower=0, upper=1> y;    // Test results [RDT, PCR, culture]
  array[N] int<lower=1, upper=3> test_cat;
  vector<lower=0, upper=1>[N] prior_chol;  // Prior probability of being cholera positive
  
  // Covarites
  matrix[N, M_rdt_sens] X_rdt_sens;
  matrix[N, M_pcr_sens] X_pcr_sens;
  matrix[N, M_culture_sens] X_culture_sens;
  matrix[N, M_rdt_spec] X_rdt_spec;
  matrix[N, M_pcr_spec] X_pcr_spec;
  matrix[N, M_culture_spec] X_culture_spec;
  matrix[N, M_chol] X_chol;
  
  matrix[U_rdt_sens, M_rdt_sens] X_u_rdt_sens;
  matrix[U_pcr_sens, M_pcr_sens] X_u_pcr_sens;
  matrix[U_culture_sens, M_culture_sens] X_u_culture_sens;
  matrix[U_rdt_spec, M_rdt_spec] X_u_rdt_spec;
  matrix[U_pcr_spec, M_pcr_spec] X_u_pcr_spec;
  matrix[U_culture_spec, M_culture_spec] X_u_culture_spec;
  
  // Priors for test performance intercepts
  array[3] real logit_sens_prior_mu;    // chol test sensitivity prior{RDT first half, RDT second half, PCR, culture}
  array[3] real<lower=0> logit_sens_prior_sd;   
  array[3] real logit_spec_prior_mu;    // chol test specificity prior{RDT first half, RDT second half, PCR, culture}
  array[3] real<lower=0> logit_spec_prior_sd;   
  
  real<lower=0> sd_beta;
  
  // Time mappings for outputs
  int<lower=1> N_month;
  array[N] int<lower=1,upper=N_month> map_month;
  array[N_month] real<lower=0, upper=1> frac_B;
  
  
  // Post-stratification matrices
  matrix[U_poststrat, U_rdt_sens] X_poststrat_rdt_sens;
  matrix[U_poststrat, U_rdt_spec] X_poststrat_rdt_spec;
  matrix[U_poststrat, U_pcr_sens] X_poststrat_pcr_sens;
  matrix[U_poststrat, U_pcr_spec] X_poststrat_pcr_spec;
  matrix[U_poststrat, U_culture_sens] X_poststrat_culture_sens;
  matrix[U_poststrat, U_culture_spec] X_poststrat_culture_spec;
  
}

transformed data {
  array[4, 8] int fake_data_rdt;    // Fake data to marginalize out unknwon pcr and culture tests
  array[4, N] real priors_A;        // priors for culture results when unavailable 
  array[4, N] real priors_B;        // priors for PCR and culture results when unavailable
  array[4, N] real log_priors_A;    // Log of priors
  array[4, N] real log_priors_B;
  array[3] real sens_prior_mu;      // Mean of sens prior. indices: 1:RDT first period, 2:RDT second period, 3:PCR, 4:culture
  array[3] real spec_prior_mu;      // Mean of spec prior.
  array[N, 8] int y_full;           // test result data in wide format (8 combinations of test results)
  vector[N] logit_prior_chol;       // Offset for cholera infection probability
  
  for (i in 1:N) {
    for (j in 1:8) {
      y_full[i, j] = 0;
    }
  }
  
  for (i in 1:N) {
    logit_prior_chol[i] = logit(prior_chol[i]);
  }
  
  for (i in 1:N) {
    
    if (y[i, 1] == 1 && y[i, 2] == 1 && y[i, 3] == 1) {
      // 1. RDT+ PCR+ CUL+
      y_full[i, 1] = 1;
    } else if (y[i, 1] == 1 && y[i, 2] == 1 && y[i, 3] == 0) {
      // 2. RTD+ PCR+ CUL-
      y_full[i, 2] = 1;
    } else if (y[i, 1] == 1 && y[i, 2] == 0 && y[i, 3] == 1) {
      // 3. RDT+ PCR- CUL+
      y_full[i, 3] = 1;
    } else if (y[i, 1] == 1 && y[i, 2] == 0 && y[i, 3] == 0) {
      // 4. RDT+ PCR- CUL-
      y_full[i, 4] = 1;
    } else if (y[i, 1] == 0 && y[i, 2] == 1 && y[i, 3] == 1) {
      // 5. RDT- PCR+ CUL+
      y_full[i, 5] = 1;
    } else if (y[i, 1] == 0 && y[i, 2] == 1 && y[i, 3] == 0) {
      // 6. RDT- PCR+ CUL-
      y_full[i, 6] = 1;
    } else if (y[i, 1] == 0 && y[i, 2] == 0 && y[i, 3] == 1) {
      // 7. RDT- PCR- CUL+
      y_full[i, 7] = 1;
    } else if (y[i, 1] == 0 && y[i, 2] == 0 && y[i, 3] == 0) {
      // 8. RDT- PCR- CUL-
      y_full[i, 8] = 1;
    } 
  }
  
  for (i in 1:3) {
    sens_prior_mu[i] = inv_logit(logit_sens_prior_mu[i]);
    spec_prior_mu[i] = inv_logit(logit_spec_prior_mu[i]);
  }
  
  for (i in 1:N) {
    // Prior probability of a neg RDT
    real p_rdtneg = (1-sens_prior_mu[1]) * prior_chol[i] + spec_prior_mu[1] * (1-prior_chol[i]);
    // Prior probability of a neg RDT and pos PCR
    real p_rdtneg_pcrpos = (1-sens_prior_mu[1]) * sens_prior_mu[2] * prior_chol[i] + spec_prior_mu[1] * (1-spec_prior_mu[2]) * (1-prior_chol[i]);
    // Prior probability of a neg RDT and neg PCR
    real p_rdtneg_pcrneg = (1-sens_prior_mu[1]) * (1-sens_prior_mu[2]) * prior_chol[i] + spec_prior_mu[1] * spec_prior_mu[2] * (1-prior_chol[i]);
    
    // Case A
    // 1. RDT- PCR+ CUL+
    priors_A[1, i] = sens_prior_mu[2] * sens_prior_mu[3] * (1-sens_prior_mu[1]) * prior_chol[i] +
    (1-spec_prior_mu[2]) * (1-spec_prior_mu[3]) * spec_prior_mu[1] * (1-prior_chol[i]);
    // 2. RDT- PCR+ CUL-
    priors_A[2, i] = sens_prior_mu[2] * (1-sens_prior_mu[3]) * (1-sens_prior_mu[1]) * prior_chol[i] +
    (1-spec_prior_mu[2]) * spec_prior_mu[3] * spec_prior_mu[1] * (1-prior_chol[i]);
    // 3. RDT- PCR- CUL+
    priors_A[3, i] = (1-sens_prior_mu[2]) * sens_prior_mu[3] * (1-sens_prior_mu[1]) * prior_chol[i] +
    spec_prior_mu[2] * (1-spec_prior_mu[3]) * spec_prior_mu[1] * (1-prior_chol[i]); 
    // 4. RDT- PCR- CUL-
    priors_A[4, i] = (1-sens_prior_mu[2]) * (1-sens_prior_mu[3]) * (1-sens_prior_mu[1]) * prior_chol[i] +
    spec_prior_mu[2] * spec_prior_mu[3] * spec_prior_mu[1] * (1-prior_chol[i]);
    
    
    // Case B
    // Initialize at unnormalized values
    for (k in 1:4) {
      priors_B[k, i] = priors_A[k, i];
    }
    
    // Normalize by respective values
    for (k in 1:4) {
      priors_A[k, i] = priors_A[k, i]/p_rdtneg;
      if (k <= 2) {
        priors_B[k, i] = priors_B[k, i]/p_rdtneg_pcrpos;
      } else {
        priors_B[k, i] = priors_B[k, i]/p_rdtneg_pcrneg;
      }
      
      // Pre-compute logs
      log_priors_A[k, i] = log(priors_A[k, i]);
      log_priors_B[k, i] = log(priors_B[k, i]);
    }
  }
  
  // Initialize fake data for rdt only measurements
  // We know these are only for negative samples
  for (i in 1:4) {
    for (j in 1:8) {
      fake_data_rdt[i, j] = 0;
    }
    fake_data_rdt[i, i+4] = 1;
  }
}

parameters {
  vector[M_rdt_sens] beta_rdt_sens;
  vector[M_rdt_spec] beta_rdt_spec;
  vector[M_pcr_sens] beta_pcr_sens;
  vector[M_pcr_spec] beta_pcr_spec;
  vector[M_culture_sens] beta_culture_sens;
  vector[M_culture_spec] beta_culture_spec;
  vector[M_chol] beta_chol;
}
transformed parameters {
  vector[N] p_chol = inv_logit(logit_prior_chol + X_chol * beta_chol);
}
model {
  array[8, N] real p;    // probability of each test result combination
  vector[N] sens_rdt;
  vector[N] spec_rdt;
  vector[N] sens_pcr;
  vector[N] spec_pcr;
  vector[N] sens_culture;
  vector[N] spec_culture;
  
  
  sens_rdt = inv_logit(X_rdt_sens * beta_rdt_sens);
  spec_rdt = inv_logit(X_rdt_spec * beta_rdt_spec);
  sens_pcr = inv_logit(X_pcr_sens * beta_pcr_sens);
  spec_pcr = inv_logit(X_pcr_spec * beta_pcr_spec);
  sens_culture = inv_logit(X_culture_sens * beta_culture_sens);
  spec_culture = inv_logit(X_culture_spec * beta_culture_spec);
  
  // Compute probabilities for each test result combination
  // We here use all three test type results:
  // 1. RDT+ PCR+ CUL+
  // 2. RTD+ PCR+ CUL-
  // 3. RDT+ PCR- CUL+
  // 4. RDT+ PCR- CUL-
  // 5. RDT- PCR+ CUL+
  // 6. RDT- PCR+ CUL-
  // 7. RDT- PCR- CUL+
  // 8. RDT- PCR- CUL-
  for (i in 1:N) {
    row_vector[3] sens = [sens_rdt[i], sens_pcr[i], sens_culture[i]];
    row_vector[3] spec = [spec_rdt[i], spec_pcr[i], spec_culture[i]];
    real x = p_chol[i];
    p[1, i] = prod(sens) * x + prod(1-spec) * (1-x);
    p[2, i] = sens[1]*sens[2]*(1-sens[3]) * x + (1-spec[1])*(1-spec[2])*spec[3] * (1-x);
    p[3, i] = sens[1]*(1-sens[2])*sens[3] * x + (1-spec[1])*spec[2]*(1-spec[3]) * (1-x);
    p[4, i] = sens[1]*(1-sens[2])*(1-sens[3]) * x + (1-spec[1])*spec[2]*spec[3] * (1-x);
    p[5, i] = (1-sens[1])*sens[2]*sens[3] * x + spec[1]*(1-spec[2])*(1-spec[3]) * (1-x);
    p[6, i] = (1-sens[1])*sens[2]*(1-sens[3]) * x + spec[1]*(1-spec[2])*spec[3] * (1-x);
    p[7, i] = (1-sens[1])*(1-sens[2])*sens[3] * x + spec[1]*spec[2]*(1-spec[3]) * (1-x);
    p[8, i] = prod(1-sens) * x + prod(spec) * (1-x);
  }
  
  // Multinomial log-likelihoods
  for (i in 1:N) {
    array[4] real p_rdt;
    matrix[2, 2] p_rdt_pcr;
    vector[8] p_vec = to_vector(p[, i]);
    
    if (test_cat[i] == 3) {
      // Likelihood of RDT- only observations (missing PCR and culture)
      for (k in 1:4) {
        p_rdt[k] = multinomial_lpmf(fake_data_rdt[k, ]| p_vec) + log_priors_A[k, i];
      }
      target += log_sum_exp(p_rdt);
      
    } else if (test_cat[i] == 2) {
      // Likelihood of RDT- and PCR observations (missing culture)
      {
        int k;
        if (y[i, 2] == 1) {
          k = 1;
        } else {
          k = 2;
        }
        p_rdt_pcr[k, 1] = multinomial_lpmf(fake_data_rdt[(k-1)*2+1, ]| p_vec) + log_priors_B[(k-1)*2+1, i];
        p_rdt_pcr[k, 2] = multinomial_lpmf(fake_data_rdt[(k-1)*2+2, ]| p_vec) + log_priors_B[(k-1)*2+2, i];
        target += log_sum_exp(p_rdt_pcr[k, ]);
      }
      
    } else {
      // Complete observations (RDT and PCR and culture)
      target += multinomial_lpmf(y_full[i, ]| p_vec);
    }
  }
  
  
  
  // Priors regression parameters
  beta_rdt_sens[1] ~ normal(logit_sens_prior_mu[1], logit_sens_prior_sd[1]);
  beta_rdt_spec[1] ~ normal(logit_spec_prior_mu[1], logit_spec_prior_sd[1]);
  beta_pcr_sens[1] ~ normal(logit_sens_prior_mu[2], logit_sens_prior_sd[2]);
  beta_pcr_spec[1] ~ normal(logit_spec_prior_mu[2], logit_spec_prior_sd[2]);
  beta_culture_sens[1] ~ normal(logit_sens_prior_mu[3], logit_sens_prior_sd[3]);
  beta_culture_spec[1] ~ normal(logit_spec_prior_mu[3], logit_spec_prior_sd[3]);
  
  if (M_rdt_sens > 1) {
    beta_rdt_sens[2:M_rdt_sens] ~ normal(0, sd_beta);
  }
  
  if (M_rdt_spec > 1) {
    beta_rdt_spec[2:M_rdt_spec] ~ normal(0, sd_beta);
  }
  
  if (M_pcr_sens > 1) {
    beta_pcr_sens[2:M_pcr_sens] ~ normal(0, sd_beta);
  }
  
  if (M_pcr_spec > 1) {
    beta_pcr_spec[2:M_pcr_spec] ~ normal(0, sd_beta);
  }
  
  if (M_culture_sens > 1) {
    beta_culture_sens[2:M_culture_sens] ~ normal(0, sd_beta);
  }
  
  if (M_culture_spec > 1) {
    beta_culture_spec[2:M_culture_spec] ~ normal(0, 1e-2);
  }
  
  // Prior on reg coefficients of probability of cholera
  beta_chol ~ normal(0, .5);
}

generated quantities {
  
  // Covariate effect sizes
  vector[M_rdt_sens-1] eff_size_rdt_sens;
  vector[M_rdt_spec-1] eff_size_rdt_spec;
  vector[M_pcr_sens-1] eff_size_pcr_sens;
  vector[M_pcr_spec-1] eff_size_pcr_spec;
  vector[M_culture_sens-1] eff_size_culture_sens;
  vector[M_culture_spec-1] eff_size_culture_spec;
  
  // Generated test performances
  vector[U_rdt_sens] gen_rdt_sens;
  vector[U_rdt_spec] gen_rdt_spec;
  vector[U_pcr_sens] gen_pcr_sens;
  vector[U_pcr_spec] gen_pcr_spec;
  vector[U_culture_sens] gen_culture_sens;
  vector[U_culture_spec] gen_culture_spec;
  
  // Post-stratified estimates for variable of interest
  vector[U_poststrat] gen_rdt_sens_poststrat;
  vector[U_poststrat] gen_rdt_spec_poststrat;
  vector[U_poststrat] gen_pcr_sens_poststrat;
  vector[U_poststrat] gen_pcr_spec_poststrat;
  vector[U_poststrat] gen_culture_sens_poststrat;
  vector[U_poststrat] gen_culture_spec_poststrat;
  
  // Simulated test outputs and categories
  array[N] int gen_cat;
  array[N_month] int n_A;
  array[N_month] int n_B;
  array[N_month] int n_C;
  array[N_month] int n_B_concordant;    // number of RDT- PCR- in the absence of culture (B)
  array[N_month] int n_B_rdt_FN;        // number of PCR+ among RDT- (RDT false negative) when PCR available (B)
  array[N_month] int n_C_concordant;    // number of RDT+ PCR+ culture+ (C)
  array[N_month] int n_C_rdt_FP;        // number of PCR+ and culture+ among RDT+ (RDT false positive) when all three tests available (C) 
  array[N_month] int n_C_pcr_FN;        // number of PCR- and culture+ among RDT+ (PCR false negative) when all three tests available (C) 
  array[N_month] int n_C_culture_FN;    // number of culture- and PCR+ among RDT+ (culture false negative) when all three tests available (C) 
  // array[N_month] real frac_A_pos;           // fraction of RDT+ (all samples are tested with RDT)
  // array[N_month] real frac_B_rdt_FN;        // fraction of PCR+ among RDT- (RDT false negative) when PCR available (B)
  // array[N_month] real frac_C_rdt_FP;        // fraction of PCR+ and culture+ among RDT+ (RDT false positive) when all three tests available (C) 
  // array[N_month] real frac_C_pcr_FN;        // fraction of PCR- and culture+ among RDT+ (PCR false negative) when all three tests available (C) 
  // array[N_month] real frac_C_culture_FN;    // fraction of culture- and PCR+ among RDT+ (culture false negative) when all three tests available (C) 
  
  // Log-likelihoods for model comparison
  vector[N] log_lik;
  
  
  // Generated test performances
  gen_rdt_sens = inv_logit(X_u_rdt_sens * beta_rdt_sens);
  gen_rdt_spec = inv_logit(X_u_rdt_spec * beta_rdt_spec);
  gen_pcr_sens = inv_logit(X_u_pcr_sens * beta_pcr_sens);
  gen_pcr_spec = inv_logit(X_u_pcr_spec * beta_pcr_spec);
  gen_culture_sens = inv_logit(X_u_culture_sens * beta_culture_sens);
  gen_culture_spec = inv_logit(X_u_culture_spec * beta_culture_spec);
  
  // Post-stratified estimates
  gen_rdt_sens_poststrat = X_poststrat_rdt_sens * gen_rdt_sens;
  gen_rdt_spec_poststrat = X_poststrat_rdt_spec * gen_rdt_spec;
  gen_pcr_sens_poststrat = X_poststrat_pcr_sens * gen_pcr_sens;
  gen_pcr_spec_poststrat = X_poststrat_pcr_spec * gen_pcr_spec;
  gen_culture_sens_poststrat = X_poststrat_culture_sens * gen_culture_sens;
  gen_culture_spec_poststrat = X_poststrat_culture_spec * gen_culture_spec;
  
  
  
  // Effect sizes
  if (M_rdt_sens > 1) {
    eff_size_rdt_sens = exp(beta_rdt_sens[2:M_rdt_sens]);
  }
  if (M_rdt_spec > 1) {
    eff_size_rdt_spec = exp(beta_rdt_spec[2:M_rdt_spec]);
  }
  if (M_pcr_sens > 1) {
    eff_size_pcr_sens = exp(beta_pcr_sens[2:M_pcr_sens]);
  }
  if (M_pcr_spec > 1) {
    eff_size_pcr_spec = exp(beta_pcr_spec[2:M_pcr_spec]);
  }
  if (M_culture_sens> 1) {
    eff_size_culture_sens = exp(beta_culture_sens[2:M_culture_sens]);
  }
  if (M_culture_spec > 1) {
    eff_size_culture_spec = exp(beta_culture_spec[2:M_culture_spec]);
  }
  
  { 
    
    array[8, N] real p;    // probability of each test result combination
    vector[N] sens_rdt;
    vector[N] spec_rdt;
    vector[N] sens_pcr;
    vector[N] spec_pcr;
    vector[N] sens_culture;
    vector[N] spec_culture;
    array[N, 8] int gen_samples;
    
    
    sens_rdt = inv_logit(X_rdt_sens * beta_rdt_sens);
    spec_rdt = inv_logit(X_rdt_spec * beta_rdt_spec);
    sens_pcr = inv_logit(X_pcr_sens * beta_pcr_sens);
    spec_pcr = inv_logit(X_pcr_spec * beta_pcr_spec);
    sens_culture = inv_logit(X_culture_sens * beta_culture_sens);
    spec_culture = inv_logit(X_culture_spec * beta_culture_spec);
    
    // Compute probabilities for each test result combination
    // We here use all three test type results:
    // 1. RDT+ PCR+ CUL+
    // 2. RTD+ PCR+ CUL-
    // 3. RDT+ PCR- CUL+
    // 4. RDT+ PCR- CUL-
    // 5. RDT- PCR+ CUL+
    // 6. RDT- PCR+ CUL-
    // 7. RDT- PCR- CUL+
    // 8. RDT- PCR- CUL-
    for (i in 1:N) {
      row_vector[3] sens = [sens_rdt[i], sens_pcr[i], sens_culture[i]];
      row_vector[3] spec = [spec_rdt[i], spec_pcr[i], spec_culture[i]];
      real x = p_chol[i];
      p[1, i] = prod(sens) * x + prod(1-spec) * (1-x);
      p[2, i] = sens[1]*sens[2]*(1-sens[3]) * x + (1-spec[1])*(1-spec[2])*spec[3] * (1-x);
      p[3, i] = sens[1]*(1-sens[2])*sens[3] * x + (1-spec[1])*spec[2]*(1-spec[3]) * (1-x);
      p[4, i] = sens[1]*(1-sens[2])*(1-sens[3]) * x + (1-spec[1])*spec[2]*spec[3] * (1-x);
      p[5, i] = (1-sens[1])*sens[2]*sens[3] * x + spec[1]*(1-spec[2])*(1-spec[3]) * (1-x);
      p[6, i] = (1-sens[1])*sens[2]*(1-sens[3]) * x + spec[1]*(1-spec[2])*spec[3] * (1-x);
      p[7, i] = (1-sens[1])*(1-sens[2])*sens[3] * x + spec[1]*spec[2]*(1-spec[3]) * (1-x);
      p[8, i] = prod(1-sens) * x + prod(spec) * (1-x);
      
      // Generate samples
      gen_samples[i, ] = multinomial_rng(to_vector(p[, i]), 1);
      
      // Define category
      for (j in 1:8) {
        if(gen_samples[i, j] == 1) {
          gen_cat[i] = j;
        }
      }
    }
    
    // Compute statistics for posterior retrodictive checks
    {
      // Initialize numbers in each category
      for (m in 1:N_month) {
        n_A[m] = 0;
        n_B[m] = 0;
        n_C[m] = 0;
        n_B_concordant[m] = 0;
        n_B_rdt_FN[m] = 0;
        n_C_concordant[m] = 0;
        n_C_rdt_FP[m] = 0;
        n_C_pcr_FN[m] = 0;
        n_C_culture_FN[m] = 0;
      }
      
      for (i in 1:N) {
        int m = map_month[i];
        
        if (gen_cat[i] < 5) {
          n_C[m] += 1;
          
          // RDT+ counts
          if (gen_cat[i] == 1) {
            n_C_concordant[m] += 1;
          } else if (gen_cat[i] == 2) {
            n_C_culture_FN[m] += 1;
          } else if (gen_cat[i] == 3) {
            n_C_pcr_FN[m] += 1;
          } else if (gen_cat[i] == 4) {
            n_C_rdt_FP[m] += 1;
          }
        } else if (frac_B[m] > uniform_rng(0, 1)) {
          n_B[m] += 1;
          // Simulate sampling of PCR-tested samples
          //  - categories 5-8 are RDT- susceptible of being tested by PCR
          //  - categories 5-6 are RDT- and PCR+
          if (gen_cat[i] < 7) {
            n_B_rdt_FN[m] += 1;
          } else {
            n_B_concordant[m] += 1;
          }
        } else {
          n_A[m] += 1;
        }
      }
    }
    
    // Compute log-likelihoods
    // Multinomial log-likelihoods
    for (i in 1:N) {
      array[4] real p_rdt;
      matrix[2, 2] p_rdt_pcr;
      vector[8] p_vec = to_vector(p[, i]);
      
      if (test_cat[i] == 3) {
        // Likelihood of RDT- only observations (missing PCR and culture)
        for (k in 1:4) {
          p_rdt[k] = multinomial_lpmf(fake_data_rdt[k, ]| p_vec) + log_priors_A[k, i];
        }
        
        log_lik[i] = log_sum_exp(p_rdt);
        
      } else if (test_cat[i] == 2) {
        // Likelihood of RDT- and PCR observations (missing culture)
        {
          int k;
          if (y[i, 2] == 1) {
            k = 1;
          } else {
            k = 2;
          }
          p_rdt_pcr[k, 1] = multinomial_lpmf(fake_data_rdt[(k-1)*2+1, ]| p_vec) + log_priors_B[(k-1)*2+1, i];
          p_rdt_pcr[k, 2] = multinomial_lpmf(fake_data_rdt[(k-1)*2+2, ]| p_vec) + log_priors_B[(k-1)*2+2, i];
          log_lik[i] = log_sum_exp(p_rdt_pcr[k, ]);
        }
        
      } else {
        // Complete observations (RDT and PCR and culture)
        log_lik[i] = multinomial_lpmf(y_full[i, ]| p_vec);
      }
    }
    
  }
}
