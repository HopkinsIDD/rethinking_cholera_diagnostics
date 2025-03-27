# Utility functions
addEpiWeek <- function(df, date_col = "date") {
  
  df %>% 
    mutate(week = lubridate::epiweek(!!rlang::sym(date_col)),
           year = lubridate::epiyear(!!rlang::sym(date_col)),
           epiweek = str_c(year, week, sep = "-"),
           epiweek_date = as.Date(str_c(epiweek, "-1"), format = "%Y-%U-%u"))
}

getAgeCuts <- function() {
  c(0, 5, Inf)
}

getAgeCuts2 <- function() {
  c(0, 5, 10, 15, 25, 35, 45, 55, 65, 75, 85, Inf)
}

logit <- function(x) {
  log(x/(1-x))
}

get_age_group <- function(age, 
                          num_breaks = 1,
                          break1 = 5, 
                          break2 = 65) {
  # 3 age categories
  if (num_breaks == 2) {
    grp <- ifelse(age < break1, 1,
                  ifelse(age >= break1 & age < break2, 2, 
                         3)
    )
    # 2 age categories
  } else if (num_breaks == 1) {
    grp <- ifelse(age < break1, 1, 2)
    # error if different number of categories requested
  } else {
    stop('num_breaks must be 1 or 2')
  }
  return(grp)
}
getSitakundatData <- function(data_path = "serochit_data/census_data_sex_age_sitakunda.csv",
                              age_categories) {
  read_csv(data_path) %>% 
    janitor::clean_names() %>% 
    select(age, pop = both_sexes_population) %>% 
    mutate(age = case_when(age == "80+" ~ "80-100", 
                           T ~ age)) %>% 
    mutate(al = str_extract(age, "[0-9]+(?=-)") %>% as.numeric(),
           ar = str_extract(age, "(?<=-)[0-9]+") %>% as.numeric(),
           pop_2021 = pop * 1.015^10,
           pop_2021 = case_when(age == "0-4" ~ pop_2021 * .8, 
                                T ~ pop_2021),
           age_cat = get_age_group(al, num_breaks = 1)) %>% 
    group_by(age_cat) %>% 
    summarise(pop = sum(pop_2021)) %>% 
    ungroup() %>% 
    mutate(age_cat = age_categories[age_cat])
}

getAllTests <- function() {
  str_c(rep(c("rdt", "pcr", "culture"), each = 2), rep(c("sens", "spec"), time = 3), sep = "_")
}



# Data descriptives -------------------------------------------------------

computeFracB <- function(lab_data) {
  lab_data %>% 
    filter(is.na(cult_res)) %>% 
    group_by(month) %>% 
    summarise(n_tot = n(),
              n_pcr = sum(!is.na(pcr_res)),
              frac_B = n_pcr/n_tot) %>% 
    mutate(yearnum = str_extract(month, "[0-9]{4}") %>% as.numeric(),
           monthnum = str_extract(month, "(?<=_)[0-9]+") %>% as.numeric()) %>% 
    arrange(yearnum, monthnum)
}

# Covariates --------------------------------------------------------------

makeCovarMats <- function(lab_data, 
                          eqs) {
  map(eqs, 
      function(x) {
        model.matrix(data = lab_data, x)
      }) %>% 
    magrittr::set_names(names(eqs))
}


makeUCovarDFs <- function(all_u_covar, 
                          eqs) {
  map(eqs, 
      function(x) {
        all_u_covar %>% 
          select(
            # date, 
            one_of(as.character(x) %>% 
                     str_split(" \\+ ") %>% 
                     .[[2]] %>% 
                     str_split(":") %>% 
                     unlist())) %>% 
          distinct()
      }) %>% 
    magrittr::set_names(names(eqs))
}

makeUCovarMats <- function(u_covar_dfs, 
                           eqs,
                           variable) {
  map(1:length(u_covar_dfs), 
      function(x) {
        model.matrix(data = u_covar_dfs[[x]], eqs[[x]]) %>% 
          {
            unique(.)
            # if(variable == "mean") {
            #   unique(.)
            # } else {
            #   .
            # }
          }
      }) %>% 
    magrittr::set_names(names(eqs))
  
}

defineRDTPeriod <- function(x) {
  ifelse(x < as.Date("2021-06-29"), "Pre-batch", "Post-batch")
}

defineEpiPeriod <- function(x) {
  # Align with RDT evaluation study
  case_when(lubridate::month(x) %in% c(3:6) ~ "pre-monsoon",
            lubridate::month(x) %in% c(7:8) ~ "monsoon",
            lubridate::month(x) %in% c(9:10) ~ "post-monsoon",
            T ~ "winter")
}

getEpiPeriodLimits <- function() {
  tibble(date = seq.Date(as.Date("2021-01-01"), as.Date("2022-08-30"), by = "1 days")) %>%
    mutate(season = map_chr(date, ~defineEpiPeriod(.)))  %>%
    mutate(season = factor(season, levels = c("pre-monsoon", "monsoon", "post-monsoon", "winter"))) %>%
    mutate(seas_diff = season != lag(season, 1),
           seas_diff = ifelse(is.na(seas_diff), 0, seas_diff),
           period = cumsum(seas_diff)) %>% 
    group_by(season, period) %>%
    summarise(tl = min(date),
              tr = max(date)) %>%
    ungroup()
}


# Post-stratification -----------------------------------------------------

#' getAllPostStratMatrices
#'
#' @param variable 
#' @param u_covar_mats 
#' @param u_covar_dfs 
#' @param eqs 
#' @param lab_data 
#'
#' @return
#' @export
#'
#' @examples
getAllPostStratMatrices <- function(variable,
                                    u_covar_mats,
                                    u_covar_dfs,
                                    eqs,
                                    lab_data) {
  
  map(1:length(eqs), function(i) {
    getPostStratMatrix(variable,
                       what = names(eqs)[i],
                       eq = eqs[[i]], 
                       u_covar_df = u_covar_dfs[[i]],
                       u_covar_mat = u_covar_mats[[i]],
                       lab_data  = lab_data)
  }) %>% 
    magrittr::set_names(names(eqs))
  
}

#' getPostStratMatrix
#'
#' @param variable 
#' @param what 
#' @param eq 
#' @param u_covar_df 
#' @param u_covar_mat 
#' @param lab_data 
#'
#' @return
#' @export
#'
#' @examples
getPostStratMatrix <- function(variable,
                               what,
                               eq,
                               u_covar_df,
                               u_covar_mat,
                               lab_data) {
  
  
  # Dimension along which to post-stratify
  poststrat_col <- getDAGDict()[variable] %>% as.character()
  if (is.na(poststrat_col)) {
    cat("-- No stratification for", variable, "in", what, "(not found in data)\n")
    return(matrix(1, nrow = 1, ncol = 1))
  }
  
  u_poststrat_strata <- unique(u_covar_df[[poststrat_col]]) %>% sort()
  U_poststrat <- length(u_poststrat_strata)
  
  vars <- eq %>% 
    as.character() %>% 
    .[[2]] %>% 
    str_split("\\+") %>% .[[1]] %>% 
    str_subset(poststrat_col, negate = T) %>% 
    str_subset("1", negate = T) %>% 
    str_trim()
  
  if (length(vars) == 0 & is.null(u_poststrat_strata)) {
    cat("-- No stratification for", what, "\n")
    return(matrix(1, nrow = 1, ncol = 1))
  }
  
  if (any(str_detect(vars, "std"))) {
    # Map observed numeric values to closest generated quantity
    num_vars <- str_subset(vars, "std")
    
    for (v in num_vars) {
      cat("-- Re-categorizing", v, "into unique covariate values for genquant poststratification.\n")
      ref_vars <- u_covar_df[[v]] %>% unique()
      closest_ind <- map_dbl(lab_data[[v]], function(x) {which.min(abs(ref_vars - x))})
      # Overwrite values
      lab_data[[v]] <- ref_vars[closest_ind]
    }
  }
  
  # Compute post-stratified estimates for each level.
  # For categorical variables these are computed directly from the data
  if (!str_detect(poststrat_col, "std")) {
    # Compute frequencies
    poststrat_freq <- lab_data %>% 
      group_by_at(c(poststrat_col, vars)) %>% 
      tally() %>% 
      group_by_at(poststrat_col) %>% 
      mutate(frac = n/sum(n)) %>% 
      ungroup()
    
  } else {
    # Round to allow matching to work
    u_poststrat_strata <- round(u_poststrat_strata*1e5)/1e5
    u_covar_mat[, poststrat_col] <- round(u_covar_mat[, poststrat_col]*1e5)/1e5
    
    # Compute frequencies
    poststrat_freq_simple <- lab_data %>% 
      group_by_at(vars) %>% 
      tally() %>% 
      ungroup() %>% 
      mutate(frac = n/sum(n)) %>% 
      ungroup()
    
    poststrat_freq <- map_df(u_poststrat_strata, function(x) {
      poststrat_freq_simple[[poststrat_col]] <- x
      poststrat_freq_simple
    })
  }
  
  # Make covariate matrix corresponding to frequency table
  covar_mat_poststrat_freq <- model.matrix(data = poststrat_freq, eq)
  
  poststrat_mat <- matrix(0, nrow = U_poststrat, ncol = nrow(u_covar_mat))
  counts <- rep(0, length(u_poststrat_strata))
  
  # Loop over poststrat_freq rows to match frequency tables
  for (i in 1:nrow(poststrat_freq)) {
    # What strata is the post-stratification for?
    ii <- which(u_poststrat_strata == poststrat_freq[[poststrat_col]][i])
    
    # Get index of row in u_covar_mat that matches
    j <- map_lgl(1:nrow(u_covar_mat), function(k){
      identical(covar_mat_poststrat_freq[i,],
                u_covar_mat[k, ])
    }) %>% 
      which()
    
    if (length(j) > 0) {
      counts[ii] <- counts[ii] + 1
      # Set the frequency
      poststrat_mat[ii, j] <- poststrat_freq$frac[i]
    }
  }
  
  # Check that rows sum to 1
  testthat::expect_equal(rowSums(poststrat_mat), rep(1, nrow(poststrat_mat)))
  rownames(poststrat_mat) <- u_poststrat_strata
  poststrat_mat
}

# Filenames ---------------------------------------------------------------

makeFilename <- function(prefix, opt, type = "rds") {
  suffix <- ifelse(is.null(opt$subset), "", str_c("_", opt$subset))
  suffix <- str_c(suffix, ifelse(!opt$fake_data, "", "_fakedata"))
  
  here::here(str_glue("generated_data/{prefix}_{opt$variable}{suffix}.{type}"))
}

makeCovarDataFilename <- function(opt) {
  makeFilename(prefix = "covar_data", opt = opt, type = "rdata")
}

makeStanDataFilename <- function(opt) {
  makeFilename(prefix = "stan_data", opt)
}


makeStanOutputFilename <- function(opt) {
  makeFilename(prefix = "stan_output", opt = opt)
}

makeGenquantFilename <- function(opt) {
  makeFilename(prefix = "genquant", opt)
}


# DAGs --------------------------------------------------------------------

makeEquations <- function(variable,
                          subset = NULL) {
  
  # Read dags
  dags <- readRDS("generated_data/dags.rds")
  
  if (variable == "age_timevary") {
    variable <- "age"
  }
  
  if (variable == "mean") {
    
    # Model with intercepts only
    res <- map(1:length(dags), ~ as.formula("~1")) %>% 
      magrittr::set_names(names(dags))
    
  } else if (variable == "time") {
    
    # Model with intercepts only
    res <- map(1:length(dags), function(x) {
      
      if (names(dags)[x] == "rdt_sens") {
        if (is.null(subset)) {
          as.formula("~1 + age_cat + epi_period + rdt_period")
        } else {
          as.formula("~1 + epi_period + rdt_period")
        }
      } else if (names(dags)[x] != "culture_spec") {
        if (is.null(subset)) {
          as.formula("~1 + age_cat + epi_period")
        } else {
          as.formula("~1 + epi_period")
        }
      } else {
        as.formula("~1")
      }
    }) %>% 
      magrittr::set_names(names(dags))
    
  } else if (variable %in% names(getDAGDict())) {
    
    # Get conditioning set
    adjust_set <- map(
      names(dags), 
      function(x) {
        dag <- dags[[x]]
        
        if (variable == "batch" & !(x %in% c("culture_spec"))) {
          dagitty::adjustedNodes(dag) <- c("period", "age")
        }
        
        if (variable == "temp" & (x %in% c("rdt_sens"))) {
          dagitty::adjustedNodes(dag) <- c("batch")
        }
        
        if ((variable %in% c("epi_period", "age")) & x %in% c("rdt_sens", "pcr_sens", "culture_sens")) {
          dagitty::adjustedNodes(dag) <- c("antibiotics")
        }
        
        if (x %in% c("rdt_sens")) {
          if (!("batch" %in% dagitty::adjustedNodes(dag))) {
            dagitty::adjustedNodes(dag) <- c(dagitty::adjustedNodes(dag), "batch")
          }
        }
        
        if (x %in% c("culture_sens")) {
          if (!("time_culture" %in% dagitty::adjustedNodes(dag)) & variable != "time_culture") {
            dagitty::adjustedNodes(dag) <- c(dagitty::adjustedNodes(dag), "time_culture")
          }
        }
        
        res <- try(dagitty::adjustmentSets(dag, 
                                           outcome = x, 
                                           exposure = c(variable), 
                                           effect = "total")) %>% 
          unlist()
        
        if(length(dagitty::adjustedNodes(dag)) > 0) {
          res <- c(res, dagitty::adjustedNodes(dag)) %>% unique()
        }
        
        res
      })
    
    if (!is.null(subset)) {
      adjust_set <- map(adjust_set, function(x){
        x <- str_subset(x, "age", negate = TRUE)
        if (length(x) == 0 || x == "") {
          return(NULL)
        } else {
          return(x)
        }
      })
    }
    
    res <- map(1:length(dags), 
               function(x) {
                 # Get all ancestors to check if variable is present
                 var_ancestors <- dagitty::ancestors(dags[[x]], names(dags)[x])
                 
                 if (!(variable %in% var_ancestors)) {
                   # If variable not among ancestors return simple intercept model
                   return(as.formula("~ 1"))
                 } else {
                   other_terms <- " + 1"
                   if (!is.null(adjust_set[[x]])) {
                     other_terms <- str_c(other_terms, " + ", str_c(parseAdjustSet(adjust_set[[x]]), collapse = " + "))
                   }
                   return(as.formula(str_glue("~ {getDAGDict()[variable]} {other_terms}")))
                 }
               }) %>% 
      magrittr::set_names(names(dags))
    
    
    # if (variable == "batch") {
    #   res$rdt_sens <- ~  age_cat + rdt_period  + epi_period
    # }
  } else {
    stop("variable ", variable, " is not among availble controls.")
  }
  
  res
}

getDAGDict <- function() {
  c("antibiotics" = "antibiotic_use",
    "period" = "epi_period",
    "batch" = "rdt_period",
    "temp" = "temp_std",
    "age" = "age_cat",
    "virbio" = "chol_std",
    "severity" = "cat_dehydration_status",
    "time_culture" = "time_to_culture_std")
}

parseAdjustSet <- function(x) {
  res <- unlist(x)
  
  map_chr(res, ~ getDAGDict()[.])
}


# Postprocessing ----------------------------------------------------------

#' getEffectSizes
#'
#' @param genquant 
#' @param covar_mats 
#' @param effects 
#'
#' @return
#' @export
#'
#' @examples
getEffectSizes <- function(genquant,
                           covar_mats,
                           effects) {
  
  genquant$summary(effects) %>%
    mutate(test = str_extract(variable, "rdt|pcr|culture") %>% factor(levels = c("rdt", "pcr", "culture")),
           what = str_extract(variable, "sens|spec"),
           id = str_extract(variable, "[0-9]+") %>% as.numeric()) %>%
    rowwise() %>%
    mutate(
      covar = colnames(covar_mats[[str_c(test, what, sep = "_")]])[id+1]
    ) 
}

#' getEffectsForVar
#'
#' @param var 
#'
#' @return
#' @export
#'
#' @examples
getEffectsForVar <- function(var) {
  
  if (var == "time_culture") {
    effects <- "eff_size_culture_sens"
  } else if (var == "batch") {
    effects <- c("eff_size_rdt_sens")
  } else {
    effects <- c(
      "eff_size_rdt_sens",
      "eff_size_pcr_sens",
      "eff_size_culture_sens",
      "eff_size_rdt_spec",
      "eff_size_pcr_spec")
  }
  effects
}


# Compare RDR against PCR and culture
computeSensSpec <- function(df, 
                            target_col, 
                            ref_col,
                            wide = TRUE) {
  res <- df %>% 
    mutate(target := !!rlang::sym(target_col),
           ref := !!rlang::sym(ref_col)) %>% 
    filter(!is.na(ref), !is.na(target)) %>% 
    summarise(
      tot_pos = sum(ref),
      false_neg = sum(!target & ref),
      tot_neg = sum(!ref),
      false_pos = sum(target & !ref)
    ) %>% 
    mutate(sens = 1 - false_neg/tot_pos,
           sens_lo = Hmisc::binconf(tot_pos - false_neg, tot_pos)[2],
           sens_hi = Hmisc::binconf(tot_pos - false_neg, tot_pos)[3],
           spec = 1 - false_pos/tot_neg,
           spec_lo = Hmisc::binconf(tot_neg - false_pos, tot_neg)[2],
           spec_hi = Hmisc::binconf(tot_neg - false_pos, tot_neg)[3])
  
  if (!wide) {
    res <- res  %>% 
      pivot_longer(cols = c(contains("sens"), contains("spec"))) %>% 
      mutate(what = str_extract(name, "sens|spec"),
             bound = str_extract(name, "lo|hi")) %>% 
      replace_na(list(bound = "mean")) %>% 
      select(-name) %>% 
      pivot_wider(values_from = "value",
                  names_from = "bound")
  }
  
  res
}


# PPV/NPV -----------------------------------------------------------------

computePPV <- function(p, sens, spec) {
  (sens * p)/(sens * p + (1-spec) * (1 - p))
}

computeNPV <- function(p, sens, spec) {
  (spec * (1-p))/(spec * (1-p) + (1-sens) *p)
}



# Figures -----------------------------------------------------------------

computeProps <- function(df, col_var) {
  df %>% 
    summarise(
      n = n(),
      n_pos = sum(!!rlang::sym(col_var)),
      frac_pos = n_pos/n) %>% 
    rowwise() %>% 
    mutate(lo = Hmisc::binconf(n_pos, n)[2],
           hi = Hmisc::binconf(n_pos, n)[3]) %>% 
    ungroup()
}


plotProps <- function(df, x_col) {
  pd <- position_dodge(width = .2)
  df %>% 
    ggplot(aes(x = !!rlang::sym(x_col), y = frac_pos)) +
    geom_point(position = pd, aes(size = n)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = pd, width = 0, alpha = .6) +
    theme_bw()
}

plotPropsColour <- function(df, x_col, colour_col) {
  pd <- position_dodge(width = .2)
  df %>% 
    ggplot(aes(x = !!rlang::sym(x_col), y = frac_pos, colour = !!rlang::sym(colour_col))) +
    geom_point(position = pd, aes(size = n)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = pd, width = 0, alpha = .6) +
    theme_bw()
}



colors_tests <- function() {
  test_colors <- c("RDT" = "#800ABF", 
                   "PCR" = "#DE6B0D", 
                   "culture" = "#308584", 
                   "fake" = "white",
                   "composite" = "black")
  test_colors
}

colors_seasons <- function() {
  c("pre-monsoon" =  "#BDA86E",
    "monsoon" = "#A11B1B",
    "post-monsoon" = "#334F1B",
    "winter" = "#6D82C7")
}


#' plot_pp
#' Plot the posterior retrodictive checks
#' 
#' @param df 
#' @param simple 
#'
#' @return
#' @export
#'
#' @examples
plot_pp <- function(df, 
                    simple = FALSE) { 
  
  p <- ggplot(df, aes(x = month_date)) +
    geom_point(aes(y = mean), alpha = .7) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0, alpha = .5) +
    geom_point(aes(y = n), col = "red", alpha = .7) +
    theme_bw() +
    facet_grid(label ~ ., scales = "free") +
    labs(x = "date", y = "Number of \nsamples") +
    theme(strip.text = element_text(size = 6.5))
  
  if (simple) {
    p +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())
  } else {
    p
  }
}
# Post-processing ---------------------------------------------------------

#' custom_summaries
#' Custom summaries to get the 95% CrI
#'
#' @return
#' @export
#'
#' @examples
custom_summaries <- function() {
  
  c(
    "mean", "median", "custom_quantile2",
    posterior::default_convergence_measures(),
    posterior::default_mcse_measures()
  )
}

#' cri_interval
#' The Credible interval to report in summaries
#' @return
#' @export
#'
#' @examples
cri_interval <- function() {
  c(0.025, 0.975)
}

#' custom_quantile2
#' quantile functoin with custom cri
#'
#' @param x 
#' @param cri 
#'
#' @return
#' @export
#'
#' @examples
custom_quantile2 <- function(x, cri = cri_interval()) {
  posterior::quantile2(x, probs = cri)
}
