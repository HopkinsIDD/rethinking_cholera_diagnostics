# DAG analysis for covariates
# https://docs.google.com/presentation/d/1ESkJrq-faQveq42JOw3CBwg7w-uJEoN0pFYrJnGnsg8/edit

# Preamble ----------------------------------------------------------------

library(dagitty)
library(tidyverse)

# Functions ---------------------------------------------------------------

makeDAG <- function(what) {
  
  if (what == "rdt_sens") {
    what_in_arrows <- str_glue(
      "batch -> {what}
      shedding -> {what}
      microbiome -> {what}
      temp -> {what}"
    )
  } else if (what == "rdt_spec") {
    what_in_arrows <- str_glue(
      "microbiome -> {what}"
    )
  } else if (what == "pcr_sens") {
    what_in_arrows <- str_glue(
      "shedding -> {what}
      microbiome -> {what}"
    )
  } else if (what == "pcr_spec") {
    what_in_arrows <- str_glue(
      "microbiome -> {what}"
    )
  } else if (what == "culture_sens") {
    what_in_arrows <- str_glue(
      "shedding -> {what}
    time_culture -> {what}
    microbiome -> {what}
      temp -> {what}"
    )
  } else if (what == "culture_spec") {
    what_in_arrows <- ""
  } else {
    stop("Please provide valid what: ", what)
  }
  
  
  dag_str <- str_glue(
    "dag {{
    period [exposure]
    temp [exposure]
    age [exposure]
    batch [exposure]
    microbiome [unobserved]
    vibrio [exposure]
    inf_dose [unobserved]
    load [unobserved]
    severity [exposure]
    shedding [unobserved]
    antibiotics [exposure]
    time_culture [exposure]

    {what} [outcome]
    
    period -> temp
    period -> microbiome
    period -> vibrio
    
    temp -> microbiome
    temp -> load
    temp -> vibrio
    age -> inf_dose
    age -> severity
    age -> microbiome
    vibrio -> inf_dose -> load -> shedding
    
    microbiome -> load
    antibiotics -> load
    antibiotics -> microbiome
    load -> severity
    load -> shedding
    age -> antibiotics

    {what_in_arrows}
    }}")
  
  out_dag <- dagitty(dag_str)
  
  outcomes(out_dag) <- what
  latents(out_dag) <- c("shedding", "inf_dose", "load", "microbiome")
  
  coords <- list(
    x = c(what = 6,
          period = 3,
          temp = 3,
          age = 1,
          vibrio = 2,
          inf_dose = 3,
          severity = 4,
          load = 4,
          shedding = 5,
          microbiome = 5,
          batch = 6,
          antibiotics = 1,
          time_culture = 5
    ),
    y = 5 - c(what = 3,
              period = 1,
              temp = 2,
              age = 3,
              vibrio = 4,
              inf_dose = 4,
              severity = 5,
              load = 4,
              shedding = 4,
              microbiome = 3,
              batch = 1,
              antibiotics = 5,
              time_culture = 1
    )
  )
  
  names(coords$x)[1] <- what
  names(coords$y)[1] <- what
  dagitty::coordinates(out_dag) <- coords
  
  out_dag
}


makeCombDag <- function() {
  dag_str <- str_glue(
    "dag {{
    period [exposure]
    temp [exposure]
    age [exposure]
    batch [exposure]
    microbiome [unobserved]
    vibrio [exposure]
    inf_dose [unobserved]
    load [unobserved]
    severity [exposure]
    shedding [unobserved]
    antibiotics [exposure]
    time_culture [exposure]

    pos [outcome]
    
    period -> temp
    period -> microbiome
    period -> vibrio
    
    temp -> microbiome
    temp -> load
    temp -> vibrio
    age -> inf_dose
    age -> severity
    age -> microbiome
    vibrio -> inf_dose -> load -> shedding
    
    microbiome -> load
    antibiotics -> load
    antibiotics -> microbiome
    load -> severity
    load -> shedding
    age -> antibiotics

    {what_in_arrows}
    }}")
}

# DAGs for test performance -----------------------------------------------

# Loop over tests
tests <- str_c(rep(c("rdt", "pcr", "culture"), each = 2), 
               rep(c("sens", "spec"), times = 3),
               sep = "_")

dags <- map(tests, ~ makeDAG(.)) %>% magrittr::set_names(tests)

saveRDS(dags, "generated_data/dags.rds")

dag_plots <- map(1:length(tests), function(x) {
  dags[[x]] %>% 
    ggdag::tidy_dagitty() %>% 
    mutate(latent = name %in% c("shedding", "inf_dose", "load", "microbiome")) %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +   
    ggdag::geom_dag_point(aes(color = latent)) +
    ggdag::geom_dag_edges(edge_colour = "gray") +
    ggdag::geom_dag_text(col = "black") +
    ggdag::theme_dag() +
    ggtitle(tests[x]) +
    theme(panel.border = element_rect(color = "black", fill = NA),
          legend.position = c(.1, .8),
          legend.key.size = unit(x = .3, units = "in"))
})

all_dags <- cowplot::plot_grid(plotlist = dag_plots,
                               nrow = 3, ncol = 2,
                               align = "tblr") +
  theme(panel.background = element_rect(fill = "white"))

ggsave(all_dags, filename = "figures/all_dags.png", width = 16, height = 14, dpi = 300)


map(names(dags), ~ try(dagitty::adjustmentSets(dags[[.]], 
                                               outcome = ., 
                                               exposure = c(""), 
                                               effect = "total")))

map(names(dags), ~ try(dagitty::impliedConditionalIndependencies(dags[[.]])))

# DAG for cholera positivity ----------------------------------------------
dag_chol <- dagitty::dagitty( 
  "dag {
    age [exposure]
    cholera [exposure]
    severity_prior [exposure]
    severity_hosp [exposure]
    shedding_prior [unobserved]
    shedding_hosp [unobserved]
    antibiotics [exposure]
    hospitalization [exposure]
    
    
    result [outcome]
    
    age -> cholera
    age -> antibiotics

    cholera -> shedding_prior
    shedding_prior -> severity_prior -> antibiotics -> shedding_hosp
    shedding_prior -> shedding_hosp
    shedding_hosp -> severity_hosp
    severity_prior -> hospitalization
    cholera -> shedding_hosp
    
    cholera -> result
    shedding_hosp -> result

    }")


coords <- list(
  x = c(
    age = 1,
    cholera = 1,
    shedding_prior = 2,
    shedding_hosp = 3,
    severity_prior = 2,
    severity_hosp = 3,
    antibiotics = 1.5,
    hospitalization = 2.5,
    result = 4
  ),
  y = 5 - c(
    age = 3,
    cholera = 1,
    shedding_prior = 2,
    shedding_hosp = 2,
    severity_prior = 3,
    severity_hosp = 3,
    antibiotics = 4,
    hospitalization = 4,
    result = 2
  )
)

latents(dag_chol) <- c("cholera", "severity_prior", "shedding_prior", "shedding_hosp")
adjustedNodes(dag_chol) <- "hospitalization"

dagitty::coordinates(dag_chol) <- coords

rethinking::drawdag(dag_chol)

dagitty::impliedConditionalIndependencies(dag_chol)

dagitty::adjustmentSets(dag_chol, exposure = "severity_hosp")

