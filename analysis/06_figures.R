# This script makes figures and tables for manuscript


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(here)
library(ggvenn)
library(cowplot)

source("analysis/utils.R")

# Functions ---------------------------------------------------------------

computeTestPerf <- function(df, 
                            target_col, 
                            ref_col,
                            what = "sens") {
  df %>% 
    mutate(target := !!rlang::sym(target_col),
           ref := !!rlang::sym(ref_col)) %>% 
    filter(!is.na(ref), !is.na(target)) %>% 
    { 
      x <- .
      if (what == "sens") {
        summarise(
          x,
          tot = sum(ref),
          errors = sum(!target & ref)
        )  
      } else {
        summarise(
          x,
          tot = sum(!ref),
          errors = sum(target & !ref)
        ) 
      }
    } %>%
    rowwise() %>% 
    mutate(perf = 1 - errors/tot,
           perf_lo = Hmisc::binconf(tot - errors, tot)[2],
           perf_hi = Hmisc::binconf(tot - errors, tot)[3],
           what = what) %>% 
    ungroup()
}

date_limits <- c("2021-01-01", "2022-08-30") %>% as.Date()



# Load data ---------------------------------------------------------------


# Load lab data and add definitions for plot
lab_data <- readRDS(here("generated_data/lab_data.rds")) %>% 
  addEpiWeek(date_col = "date_sample") %>% 
  mutate(obs_id = row_number()) %>% 
  # Define RDT and PCR results for plot
  mutate(rdt_pcr = str_c("RDT", ifelse(rdt_res, "+", "-"), " & PCR", ifelse(pcr_res, "+", "-")),
         rdt_pcr = ifelse(is.na(rdt_pcr), "RDT- & no PCR", rdt_pcr),
         rdt_pcr = factor(rdt_pcr, levels = c("RDT+ & PCR+",
                                              "RDT+ & PCR-",
                                              "RDT- & PCR+",
                                              "RDT- & PCR-",
                                              "RDT- & no PCR") %>% 
                            rev()),
         rdt_pcr_cult = case_when(is.na(cult_res) ~ str_c(as.character(rdt_pcr), " & no culture"),
                                  T ~ str_c(as.character(rdt_pcr), " & culture", ifelse(cult_res, "+", "-"))),
         rdt_pcr_cult = factor(rdt_pcr_cult, levels = c("RDT+ & PCR+ & culture+",
                                                        "RDT+ & PCR- & culture+",
                                                        "RDT+ & PCR+ & culture-",
                                                        "RDT+ & PCR- & culture-",
                                                        "RDT- & PCR+ & no culture",
                                                        "RDT- & PCR- & no culture",
                                                        "RDT- & no PCR & no culture") %>%
                                 rev())
  ) %>% 
  mutate(epi_period = factor(epi_period, levels = c("pre-monsoon", "monsoon", "post-monsoon", "winter")),
         age_cat = factor(age_cat, levels = c("[0,5)", "[5,Inf)")))

# Weather data for plot 
chittaghon_weather <- readRDS("generated_data/chittaghon_airport_weather_data_filled.rds")

# Table 1: summary table --------------------------------------------------

table1_data <- lab_data %>% 
  mutate(across(contains("_res"), function(x) {ifelse(x, "positive", "negative")}),
         antibiotic_use = map_chr(antibiotic_use, ~ifelse(., "Yes", "No")),
         age_cat = factor(age_cat,
                          levels = c("[0,5)", "[5,Inf)"),
                          labels = c("0-4", "5+"))) %>% 
  select(age_cat, 
         Sex = sex, 
         `Dehydration status` = cat_dehydration_status, 
         `Antibiotic use` = antibiotic_use, 
         `Time to culture` = time_to_culture, 
         `RDT result` = rdt_res, 
         `PCR result` = pcr_res, 
         `Culture result` = cult_res) %>% 
  gtsummary::tbl_summary(by = "age_cat") 

table1_data %>% 
  gtsummary::as_gt() %>%  
  gt::gtsave(filename = "figures/table_1.docx")



# Figure 1: time series and distribution ----------------------------------

# Make long format for ploting
long_weatherdat <- chittaghon_weather %>%
  select(date, temp, relhum, prcp) %>% 
  pivot_longer(cols = c("temp"),
               names_to = "var") %>% 
  # bind IMERG percipitation estiamtes
  bind_rows(
    readRDS("generated_data/imerg_precipitation_sitakunda.rds") %>% 
      as_tibble() %>% 
      select(date, value = imerg) %>% 
      mutate(var = "precip")
  ) %>% 
  addEpiWeek() %>% 
  group_by(epiweek_date, var) %>% 
  summarise(mean_val = mean(value),
            sum_val = sum(value)) %>% 
  ungroup() %>% 
  mutate(value = case_when(var == "precip" ~ sum_val,
                           var == "temp" ~ mean_val,
                           TRUE ~ NA)) %>% 
  rename(date = epiweek_date) %>% 
  filter(date >= min(lab_data$date_sample), date <= date_limits[2])

weather_dict <- c("season" = "season",
                  "temp" = "mean\ntemperature [C]",
                  "precip" = "weekly rainfall\n[mm]")

addSeasonRibbons <- function() {
  geom_rect(data = getEpiPeriodLimits(), 
            aes(xmin = tl, xmax = tr, ymin = 0, ymax = Inf, fill = season),
            inherit.aes = F, 
            alpha = .10) 
}

p_season <- tibble(var = "season", value = NA) %>% 
  mutate(var = factor(var, levels = c("season", "precip", "temp"))) %>% 
  ggplot(aes(x = date, y = value)) +
  geom_linerange(data = getEpiPeriodLimits() %>% 
                   mutate(var = "season"), 
                 aes(xmin = tl, xmax = tr, y = 1, lty = season, color = season),
                 inherit.aes = F)  +
  geom_vline(data = getEpiPeriodLimits(), 
             aes(xintercept = tr, lty = season, color = season),
             inherit.aes = F, lwd = .3)  +
  geom_label(data = getEpiPeriodLimits() %>% 
               rowwise() %>% 
               mutate(date_mid = mean.Date(c(tl, tr)),
                      var = "season") %>% 
               ungroup() %>% 
               mutate(season2 = case_when(season == "post-monsoon" ~ "post-\nmonsoon",
                                          TRUE ~ season)), 
             aes(x = date_mid, y = 1, label = season2, color = season),
             inherit.aes = F, size = 3) +
  facet_grid(var ~ ., scales = "free", labeller = labeller(var = weather_dict),
             switch = "y") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  guides(lty = "none", color = "none") +
  scale_color_manual(values = colors_seasons())

p_clim <- long_weatherdat %>% 
  mutate(var = factor(var, levels = c("precip", "temp"))) %>% 
  ggplot(aes(x = date, y = value)) +
  geom_vline(data = getEpiPeriodLimits(), 
             aes(xintercept = tr, lty = season, color = season),
             inherit.aes = F, lwd = .3)  +
  geom_line(data = filter(long_weatherdat, var == "temp"), color = "red") +
  geom_bar(data = filter(long_weatherdat, var == "precip"), fill = "darkblue",
           stat = "identity") +
  facet_grid(var ~ ., scales = "free", labeller = labeller(var = weather_dict),
             switch = "y") +
  theme_bw() +
  scale_x_date(limits = date_limits, date_breaks = "3 months") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.justification = "top") +
  scale_color_manual(values = colors_seasons())


p_data2 <- lab_data %>%
  ggplot(aes(x = epiweek_date))  +
  geom_vline(data = getEpiPeriodLimits(), 
             aes(xintercept = tr, lty = season, color = season),
             inherit.aes = F, lwd = .3) +
  geom_bar(aes(fill = rdt_pcr_cult,
               x = epiweek_date), 
           inherit.aes = F,
           width = 6) +
  scale_fill_manual(values = rev(c("#800ABF", "#CF8FF7", "#308584", "#40B8B4", "#857E7E", "#DE6B0D", "#F7C08C"))) +
  guides(fill = guide_legend("Test results")) +
  theme_bw() +
  labs(y = "Weekly count", x = NULL)  +
  theme_bw() +
  scale_x_date(limits = date_limits, date_breaks = "3 months") +
  facet_grid(age_cat~., switch = "y") +
  scale_color_manual(values = colors_seasons())


p_fig1AB <- plot_grid(
  p_season,
  p_clim,
  p_data2 +
    guides(fill = "none", lty = "none", color = "none"),
  ncol = 1,
  labels = c("a", NA, "b"),
  align = "v",
  axis = "lr",
  rel_heights = c(.6, 2, 4)
)

p_fig1AB


lab_pcronly <- str_glue("PCR-only\n(RDT-)")
lab_culture <-  str_glue("PCR & Culture\n(RDT+)")

modifyLabLabels <- function(df) {
  df %>% 
    mutate(rdt_pcr_cult = as.character(rdt_pcr_cult),
           rdt_pcr_cult = str_c("<", 
                                str_extract(rdt_pcr_cult, "(?<=RDT)[\\+,\\-]"),
                                ",",
                                ifelse(str_detect(rdt_pcr_cult, "no PCR"), "NA",
                                       str_extract(rdt_pcr_cult, "(?<=PCR)[\\+,\\-]")),
                                ",",
                                ifelse(str_detect(rdt_pcr_cult, "no culture"), "NA", 
                                       str_extract(rdt_pcr_cult, "(?<=culture)[\\+,\\-]")),
                                ">"),
           rdt_pcr_cult = factor(rdt_pcr_cult, 
                                 levels = c("<-,NA,NA>",
                                            "<-,-,NA>",
                                            "<-,+,NA>",
                                            "<+,-,->",
                                            "<+,+,->",
                                            "<+,-,+>",
                                            "<+,+,+>")))
}

p_data_rdt_pcr_cult <- lab_data %>% 
  filter(!is.na(cult_res)) %>% 
  count(age_cat, cult_res, rdt_pcr_cult) %>% 
  group_by(age_cat) %>% 
  mutate(frac = n/sum(n)) %>% 
  ungroup() %>% 
  mutate(what = lab_culture) %>% 
  bind_rows(
    lab_data %>% 
      filter(!is.na(pcr_res), is.na(cult_res)) %>% 
      count(age_cat, rdt_pcr_cult) %>% 
      group_by(age_cat) %>% 
      mutate(frac = n/sum(n)) %>% 
      ungroup() %>% 
      mutate(what = lab_pcronly)
  ) %>% 
  modifyLabLabels() %>% 
  mutate(what = factor(what, levels = c(lab_pcronly, lab_culture))) %>% 
  ggplot(aes(x = factor(1), y = frac)) +
  geom_bar(aes(fill = rdt_pcr_cult), 
           width = .5,
           stat = "identity") +
  scale_fill_manual(values = rev(c("#800ABF", "#CF8FF7", "#308584", "#40B8B4", "#857E7E", "#DE6B0D"))) +
  guides(fill = guide_legend("Test results\n<RDT,PCR,culture>")) +
  theme_bw() +
  facet_grid(age_cat ~ what, switch = "y") +
  labs(y = "proportion") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()#,
        # panel.grid = element_blank()
  )


p_data_rdt_pcr_cult

p_fig1 <- ggdraw(p_fig1AB +
                   theme(plot.margin = unit(c(1, 15, 1, 1), units = "lines")))  +
  draw_plot(p_data_rdt_pcr_cult +
              theme(
                strip.text.y = element_blank(),
                panel.border = element_blank(),
                plot.margin = unit(c(1, 1, 0, 0), units = "lines"),
                axis.line = element_line(),
                plot.background = element_blank()
              ),
            x = .68, y = .045, width = .31, height = .615) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave(p_fig1, filename = "figures/figure_1_new.png", 
       width = 13.5, height = 8, dpi = 300)



# Table 2: RDT estimates --------------------------------------------------
genquant_mean <- readRDS(makeGenquantFilename(variable = "mean"))
genquant_age <- readRDS(makeGenquantFilename(variable = "age"))

perf <- map_df(c("sens", "spec"), function(what) {
  bind_rows(
    genquant_age$summary(str_glue("gen_rdt_{what}_poststrat")) %>%
      mutate(age_cat = c("0-4", "5+")),
    genquant_mean$summary(str_glue("gen_rdt_{what}_poststrat")) %>%
      mutate(age_cat = "all")
  ) %>%
    mutate(what = what)
}) %>%
  mutate(age_cat = factor(age_cat, levels = c("all", "0-4", "5+")))

# Compute un-adjusted sens and spec against PCR
unadjusted_perf <- lab_data %>% 
  bind_rows(lab_data %>% mutate(age_cat = "overall")) %>% 
  group_by(age_cat) %>% 
  group_modify(function(x, y){
    computeSensSpec(df = x,
                    target = "rdt_res",
                    ref = "pcr_res",
                    wide = FALSE)
  }) %>%
  mutate(age_cat = factor(age_cat,
                          levels = c("overall", "[0,5)", "[5,Inf)"),
                          labels = c("all", "0-4", "5+"))) 



table2_data <- inner_join(
  unadjusted_perf %>% 
    mutate(total = case_when(what == "sens" ~ tot_pos,
                             T ~ tot_neg),
           errors = case_when(what == "sens" ~ false_neg,
                              T ~ false_pos)) %>% 
    mutate(text_raw = str_c(
      formatC(mean*100, format = "f", digits = 1),
      " (",
      formatC(lo*100, format = "f", digits = 1),
      "-",
      formatC(hi*100, format = "f", digits = 1),
      ")"
    )) %>% 
    select(age_cat, what, total, errors, text_raw),
  perf  %>% 
    mutate(text_adj = str_c(
      formatC(mean*100, format = "f", digits = 1),
      " (",
      formatC(q5*100, format = "f", digits = 1),
      "-",
      formatC(q95*100, format = "f", digits = 1),
      ")"
    )) %>% 
    select(age_cat, what, text_adj)
) %>% 
  ungroup()


table2_data %>% 
  filter(age_cat == "all") %>% 
  select(-age_cat) %>% 
  magrittr::set_colnames(c("Metric",
                           "N reference",
                           "N discordant",
                           "Unadjusted estimate against PCR",
                           "Adjusted estimate")) %>% 
  flextable::as_flextable(max_row = Inf) %>%
  flextable::set_table_properties(width = 1, layout = "autofit") %>%
  flextable::save_as_docx(path = "figures/table_2.docx")


# Figure 2: Post-stratification ------------------------------------------------
gq_files <- dir(here("generated_data/"), pattern = "genquant", full.names = T) %>%
  # str_subset("mean", negate = T) %>%
  # !! Change once re-run
  str_subset("time\\.", negate = T) %>% 
  str_subset("batch|temp", negate = T) %>% 
  str_subset("genquant(.)*rdata", negate = T)

test_dict <- c(
  "rdt" = "RDT",
  "pcr" = "PCR",
  "culture" = "culture"
)


# Post-stratified estimates
poststrat <- map(gq_files, function(x) {
  genquant <- readRDS(x)
  
  suffix <- str_extract(x, '(?<=genquant_)(.)*(?=\\.rds)')
  age_cat <- str_extract(suffix, "adults|children")
  if (!is.na(age_cat)) {
    var <- str_remove(suffix, str_c("_", age_cat))
  } else {
    var <- suffix
    age_cat <- "all"
  }
  
  load(here(str_glue("generated_data/covar_data_{suffix}.rdata")))
  stan_data <- readRDS(here(str_glue("generated_data/stan_data_{suffix}.rds")))
  
  map_df(c("rdt", "pcr", "culture"), function(test) {
    map_df(c("sens", "spec"), function(what) {
      if (!(test == "culture" & what == "spec")) {
        cat(x, test, what, "\n")
        ps <- genquant$summary(str_glue("gen_{test}_{what}_poststrat")) %>% 
          mutate(test = test,
                 what = what,
                 covariate = var)
        
        strata <- rownames(stan_data[[str_glue("X_poststrat_{test}_{what}")]])
        if (length(strata) > 0) {
          # Fix strata if numeric
          if (var == "time_culture") {
            strata <- as.numeric(strata) * sd(lab_data$time_to_culture, na.rm = T) + 
              mean(lab_data$time_to_culture, na.rm = T)
            strata <- as.factor(round(strata)) %>% fct_reorder(strata)
          } else if (var == "temp") {
            strata <- as.numeric(strata)
            strata <- as.factor(round(strata)) %>% fct_reorder(strata)
          } 
          ps$strata <- strata
        } else if (var == "mean") {
          ps$strata <- "overall"
          ps$covariate <- "overall"
        }
        
        ps %>% 
          mutate(what = factor(what, 
                               levels = c("sens", "spec"),
                               labels = c("sensitivity", "specificity")),
                 test = factor(test_dict[test], levels = test_dict),
                 age_cat = age_cat)
      }
    })
  })
}) 

names(poststrat) <- map_chr(gq_files, ~ str_extract(., '(?<=genquant_)(.)*(?=\\.rds)')) %>% 
  str_replace_all("mean", "overall")

poststrat  %>% 
  saveRDS(file = "generated_data/poststratified_estimates_new.rds")

label_dict <- c(
  "overall" = "overall",
  "age" = "age",
  "period" = "season",
  "antibiotics" = "antibiotics",
  "time_culture" = "time to culture"
)

ps_figures <- map(
  c("overall", "age", "period", "antibioticsGTFCC", "time_culture"), 
  function(x) {
    
    ps <- poststrat[[x]] %>% 
      mutate(strata = case_when(strata == "[0,5)" ~ "0-4",
                                strata == "[5,Inf)" ~ "5+",
                                # strata == "[65,Inf)" ~ "65+",
                                strata == "0" & str_detect(covariate, "antibiotics")  ~ "No",
                                strata == "1" & str_detect(covariate, "antibiotics")  ~ "Yes",
                                T ~ strata))
    
    if (x == "antibioticsGTFCC") {
      ps$covariate <- "antibiotic use"
    }
    
    # Add fake datapoint for specicity
    ps_fake <- tibble(mean = c(.5, 1), 
                      test = "fake", 
                      what = "specificity", 
                      covariate =  ps$covariate[1],
                      strata = ps$strata[1]) %>% 
      bind_rows(tibble(mean = c(.35, 1), 
                       test = "fake", 
                       what = "sensitivity", 
                       covariate =  ps$covariate[1],
                       strata = ps$strata[1]))
    
    
    if (x == "period") {
      ps <- mutate(ps, strata = factor(strata, levels = c("pre-monsoon", "monsoon", "post-monsoon", "winter")))
    }
    
    
    if (ps$covariate[1] %in% c("time_culture")) {
      xlab <- "days"
    } else {
      xlab <- NULL
    }
    
    if (ps$covariate[1] %in% c("overall")) {
      ylab <- "Post-stratified estimate"
    } else {
      ylab <- NULL
    }
    
    pd <- position_dodge(width = .3)
    
    if (ps$covariate[1] != c("time_culture")) {
      p <- ps %>% 
        filter(!is.na(strata)) %>% 
        ggplot(aes(x = strata, y = mean, color = test)) +
        geom_errorbar(position = pd, width = 0, aes(ymin = q5, ymax = q95)) +
        geom_point(position = pd, aes(pch = test))   +
        geom_point(data = ps_fake, alpha = 0)   +
        guides(color = "none", shape = "none") +
        theme_bw()
    } else {
      p <- ps %>% 
        filter(!is.na(strata)) %>%
        mutate(strata = as.numeric(as.character(strata))) %>%  
        ggplot(aes(x = strata, y = mean)) +
        geom_ribbon(aes(ymin = q5, ymax = q95, fill = test), alpha = .2, show.legend = TRUE) +
        geom_point(data = ps_fake %>% filter(what == "sensitivity"), 
                   aes(x = 10), alpha = 0)   +
        geom_line(aes(color = test), show.legend = TRUE) +
        geom_rug(inherit.aes = F,
                 data = lab_data %>% 
                   filter(time_to_culture > 0),
                 aes(x = time_to_culture, y = .5),
                 position = "jitter",
                 sides = "b",
                 alpha = .3) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.direction = "horizontal")
    }
    
    p <- p + 
      facet_grid(what ~ covariate, scales = "free", space = "free_x", drop = FALSE,
                 labeller = labeller(covariate = label_dict),
                 switch = "y") +
      labs(x = xlab, y = ylab) +
      # scale_y_continuous(limits = c(.3, 1)) +
      scale_color_manual(values = colors_tests(), drop = F, limits = c("RDT", "PCR", "culture")) +
      scale_fill_manual(values = colors_tests(), drop = F, limits = c("RDT", "PCR", "culture"))
    
    if (!(ps$covariate[1] %in% c("overall", "time_culture"))) {
      p <- p + theme(strip.text.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
    }
    
    if (x == "period") {
      p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    }
    
    p +
      theme(panel.grid.minor.y = element_blank())
    
  })


p_poststrat <- plot_grid(
  plotlist = ps_figures,
  nrow = 1,
  rel_widths = map_dbl(
    c("overall", "age", "period", "antibiotics", "time_culture"), 
    function(x) {
      ps <- poststrat[[x]] %>% filter(!is.na(strata))
      length(unique(ps$strata))*ifelse(x == "time_culture", .5, 1)
    }) + c(1.5, rep(0, 3), 3),
  align = "h",
  axis = "tb"
)

ggsave(p_poststrat, filename = "figures/figure_2_new.png", width = 11, height = 5)


# Figure 3: PPV/NPV -------------------------------------------------------

prior_chol_median <- lab_data %>% 
  group_by(age_cat, epi_period) %>% 
  summarise(prior_chol = median(prior_chol)) %>%
  mutate(age_cat = factor(age_cat,
                          levels = c("overall", "[0,5)", "[5,Inf)"),
                          labels = c("all", "0-4", "5+"))) 

# PPV and NPV
sens_spec_gen <- genquant_age$draws("gen_rdt_sens_poststrat") %>% 
  posterior::as_draws() %>% 
  posterior::as_draws_df() %>% 
  as_tibble() %>% 
  pivot_longer(cols = contains("post")) %>% 
  mutate(age_cat = c("0-4", "5+")[as.numeric(str_extract(name, "[1-2]"))]) %>% 
  select(-name) %>% 
  rename(sens = value) %>% 
  inner_join(
    genquant_age$draws("gen_rdt_spec_poststrat") %>% 
      posterior::as_draws() %>% 
      posterior::as_draws_df() %>% 
      as_tibble() %>% 
      pivot_longer(cols = contains("post")) %>% 
      mutate(age_cat = c("0-4", "5+")[as.numeric(str_extract(name, "[1-2]"))]) %>% 
      select(-name) %>% 
      rename(spec = value)
  ) 

# Compute PPV and NPV at mean prevalences
pv_mean <- sens_spec_gen %>% 
  inner_join(prior_chol_median) %>% 
  rowwise() %>% 
  mutate(ppv = computePPV(p = prior_chol, sens = sens, spec = spec),
         npv = computeNPV(p = prior_chol, sens = sens, spec = spec)) %>% 
  pivot_longer(cols = c("ppv", "npv")) %>% 
  group_by(age_cat, epi_period, name) %>% 
  summarise(prior_chol = mean(prior_chol),
            mean = mean(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975)) %>% 
  ungroup() %>% 
  mutate(age_cat = factor(age_cat,
                          labels = c("0-4", "5+")),
         epi_period = factor(epi_period, levels = c("pre-monsoon", "monsoon", "post-monsoon", "winter"))) 


# To text format
pv_mean %>% 
  mutate(
    prior_chol =  formatC(prior_chol*100, format = "f", digits = 1),
    text = str_c(
      formatC(mean*100, format = "f", digits = 1),
      " (",
      formatC(lo*100, format = "f", digits = 1),
      "-",
      formatC(hi*100, format = "f", digits = 1),
      ")"
    )) %>% 
  select(age_cat, epi_period, name, prior_chol, text) %>% 
  pivot_wider(names_from = "name",
              values_from = "text")

# Figure with varying levels of prevalence
# Weekly prior of cholera by age category
prior_chol_all <- lab_data %>% 
  # summarise(prior_chol = mean(prior_chol)) %>%
  mutate(age_cat = factor(age_cat,
                          levels = c("overall", "[0,5)", "[5,Inf)"),
                          labels = c("all", "0-4", "5+")),
         epi_period = factor(epi_period, levels = c("pre-monsoon", "monsoon", "post-monsoon", "winter"))) 

# Compute mean line
pv_pred <- sens_spec_gen %>% 
  inner_join(
    expand.grid(
      age_cat = unique(sens_spec_gen$age_cat),
      prior_chol = seq(1e-7,  max(lab_data$prior_chol), by = .01)
    )
  ) %>% 
  mutate(ppv = computePPV(p = prior_chol, sens = sens, spec = spec),
         npv = computeNPV(p = prior_chol, sens = sens, spec = spec)) %>% 
  pivot_longer(cols = c("ppv", "npv")) %>% 
  group_by(age_cat, name, prior_chol) %>% 
  summarise(prior_chol = mean(prior_chol),
            mean_value = mean(value),
            q025 = quantile(value, 0.025),
            q975 = quantile(value, 0.975)) %>% 
  ungroup() 

w <- .5
# Compute IQRs for prevalence
prior_chol_iqr_age <- prior_chol_all %>% 
  group_by(epi_period, age_cat) %>% 
  summarise(lo = quantile(prior_chol, .25),
            hi = quantile(prior_chol, .75),
            median = median(prior_chol))

p_fig3 <- pv_pred %>% 
  ggplot(aes(x = prior_chol)) +
  geom_rect(data = prior_chol_iqr_age, inherit.aes = F, 
            aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = epi_period),
            alpha = .15) +
  geom_vline(data = prior_chol_iqr_age, 
             aes(xintercept = median, lty = epi_period),
             lwd = .7, color = "white", alpha = .3) +
  geom_vline(data = prior_chol_iqr_age, 
             aes(xintercept = median, color = epi_period, lty = epi_period),
             lwd = .2) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .3) +
  geom_line(aes(y = mean_value), col = "black") +
  geom_rug(data = prior_chol_all, alpha = .1, aes(color = epi_period)) +
  geom_errorbar(data = pv_mean, aes(ymin = lo, ymax = hi, color = epi_period), 
                # col = "red", 
                width = .01) +
  geom_point(data = pv_mean, aes(y = mean, pch = epi_period, color = epi_period), #col = "red", #pch = 21, 
             fill = "white", size = 2.5) +
  facet_grid(name ~ age_cat, scales = "free",
             labeller = labeller(name = c("npv" = "Negative Predictive value",
                                          "ppv" = "Positive Predictive value"))) +
  theme_bw() +
  labs(x = "Cholera prevalence among AWD cases",
       y = "Negative/Positive predictive value",
       color = "Season",
       linetype = "Season",
       shape = "Season",
       fill = "Season") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_color_manual(values = colors_seasons()) +
  scale_fill_manual(values = colors_seasons())


ggsave(p_fig3, filename = "figures/figure_3.png", width = 10, height = 6, dpi = 300)

# Save for use in manuscript statistics report
save(pv_mean, sens_spec_gen, pv_pred, file = "generated_data/ppv_npv_data.rdata")

# Figure 4: simulations ---------------------------------------------------

sim_data <- readRDS("generated_data/simulated_performance_mean.rds") %>% 
  rename(what = name) %>% 
  mutate(ref = factor(ref, 
                      levels = c("culture_pcr", "pcr", "culture"),
                      labels = c("composite", "PCR", "culture")),
         what = ifelse(what == "sens", "sensitivity", "specificity"))

test_performance <- readRDS("generated_data/generated_test_performance_mean.rds")


# Load data from systematic review
sysreview <- readxl::read_xlsx("data/journal.pmed.1004286.s003.xlsx") %>% 
  janitor::clean_names() %>% 
  mutate(prop_positive = number_of_cases_positive/number_of_cases_tested,
         in_outbreak = is.na(outbreak_start))

pd <- position_dodge(width = .3)

# Compute IQRs for prevalence
prior_chol_iqr <- prior_chol_all %>% 
  group_by(epi_period) %>% 
  summarise(lo = quantile(prior_chol, .25),
            hi = quantile(prior_chol, .75),
            median = median(prior_chol))

p_fig4_v2 <- sim_data %>% 
  ggplot(aes(x = chol_prior, y = mean)) +
  geom_rect(data = prior_chol_iqr, inherit.aes = F, 
            aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = epi_period),
            alpha = .15) +
  geom_vline(data = prior_chol_iqr, inherit.aes = F, 
             aes(xintercept = median, lty = epi_period),
             lwd = .7, color = "white", alpha = .3) +
  geom_vline(data = prior_chol_iqr, inherit.aes = F, 
             aes(xintercept = median, color = epi_period, lty = epi_period),
             lwd = .2) +
  geom_hline(data = test_performance %>% 
               mutate(what = ifelse(what == "sens", "sensitivity", "specificity")) %>% 
               filter(test == "rdt"), aes(yintercept = mean),
             lty = 2, lwd = .3, col = "red") +
  # geom_line(position = pd, lwd = .3, alpha = .5) +
  geom_errorbar(aes(ymin = q5, ymax = q95), 
                width = 0, alpha = 1,
                lwd = .3, position = pd) +
  geom_point(size = 1) +
  facet_grid(what ~ ref, scales = "free") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = colors_seasons()) +
  scale_fill_manual(values = colors_seasons()) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  # scale_color_manual(values = c("#84846C", "#594444", "#42515C")) +
  labs(x = "True cholera prevalence among AWD", 
       y = "Simulated RDT performance estimate") +
  guides(fill = guide_legend("Season"),
         color = guide_legend("Season"),
         lty = guide_legend("Season")) +
  theme(legend.position = "bottom")

ggsave(p_fig4_v2, filename = "figures/figure_4.png", width = 8, height = 5.3, dpi = 300)



# Supfig S1: venns by age cat ----------------------------------------

p_venn_culture_age <- map(levels(lab_data$age_cat), function(x){
  dat <- lab_data %>% 
    filter(!is.na(cult_res), age_cat == x)
  
  dat %>% 
    select(`RDT+` = rdt_res, `PCR+` = pcr_res, `Culture+` = cult_res) %>% 
    ggvenn(text_size = 3.5, set_name_size = 5) +
    ggtitle(str_glue("Age: {x} Culture and PCR-tested sample\n(RDT+ only, N = {dat %>% nrow()})"))
}) %>% 
  plot_grid(
    plotlist = .,
    nrow = 1
  )

p_venn_pcr_age <- map(levels(lab_data$age_cat), function(x){
  dat <- lab_data %>% 
    filter(!is.na(pcr_res), age_cat == x)
  
  dat %>% 
    select(`RDT+` = rdt_res, `PCR+` = pcr_res) %>% 
    ggvenn(text_size = 3.5, set_name_size = 5) +
    ggtitle(str_glue("Age: {x} PCR-tested samples\n(N = {dat %>% nrow()})"))
}) %>% 
  plot_grid(
    plotlist = .,
    nrow = 1
  )

# Sampling diagram
# https://docs.google.com/presentation/d/1uLd6mQ0ujb1HQ3TyzbltasfukeTGlTIVBhusucnw8jQ/edit?usp=sharing

p_sfig_venns <- plot_grid(
  ggdraw() + draw_image("figures/serochit_sampling.png"),
  plot_grid(
    p_venn_pcr_age +
      theme(plot.margin = unit(c(1, 0, 2, 1), units = "lines")),
    p_venn_culture_age +
      theme(plot.margin = unit(c(1, 0, 2, 1), units = "lines")),
    ncol = 1,
    rel_heights = c(.7, 1),
    labels = c("b", "c")
  ),
  nrow = 1,
  labels = c("a", NA),
  rel_widths = c(1, 3)
)+
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave(p_sfig_venns, filename = "figures/sfigure_venns.png", width = 15, height = 10, dpi = 300)


# Supfig S2: effect sizes -------------------------------------------------


eff_sizes <-  dir(here("generated_data/"), pattern = "genquant", full.names = T) %>%
  # str_subset("mean", negate = T) %>%
  # !! Change once re-run
  str_subset("time\\.", negate = T) %>% 
  str_subset("temp", negate = T) %>% 
  str_subset("genquant(.)*rdata", negate = T) %>%
  # str_subset("adult", negate = F) %>%
  str_subset("child|adult", negate = T) %>%
  str_subset("batch|age|GTFCC|period|time") %>% 
  map_df(function(x) {
    genquant <- readRDS(x)
    
    var <- str_extract(x, '(?<=genquant_)(.)*(?=\\.rds)')
    if (var == "mean") {
      return(NULL)
    }
    
    load(here(str_glue("generated_data/covar_data_{var}.rdata")))
    
    if (var == "antibioticsGTFCC") {
      var <- "antibiotics"
    }
    
    getEffectSizes(genquant = genquant,
                   covar_mats = covar_mats,
                   effects = getEffectsForVar(var))  %>%
      filter(str_detect(covar, getDAGDict()[var]))
  })


# Save for later use
saveRDS(eff_sizes, file = "generated_data/effect_sizes.rds")

pd <- position_dodge(width = .3)
p_eff <- eff_sizes %>%
  mutate(covar_class = case_when(str_detect(covar, "age") ~ "age\n(ref:5-64)",
                                 str_detect(covar, "rdt") ~ "RDT batch < Jun 30th\n(ref:>Jun 30th)",
                                 str_detect(covar, "period") ~ "season\n(ref:mosoon)",
                                 str_detect(covar, "time_to_culture") ~ "time to\nculture\n(one sd from mean)",
                                 str_detect(covar, "antibio") ~ "antibiotics\n(ref:No)",
                                 T ~ "other"),
         covar_class = factor(covar_class, levels = c("age\n(ref:5-64)",
                                                      "season\n(ref:mosoon)", 
                                                      "antibiotics\n(ref:No)",
                                                      "time to\nculture\n(one sd from mean)",
                                                      "RDT batch < Jun 30th\n(ref:>Jun 30th)")), 
         covar = str_remove(covar, "age_cat|epi_period"),
         covar = case_when(covar == "[0,5)" ~ "0-4",
                           # covar == "[65,Inf)" ~ "65+",
                           covar == "antibiotic_use" ~ "Yes",
                           covar == "time_to_culture_std" ~ "scaled time",
                           covar == "rdt_periodPre-batch" ~ "RDT batch",
                           T ~ covar),
         covar = factor(covar, 
                        levels = c("0-4", "pre-monsoon", "post-monsoon", "winter", "Yes", "scaled time", "RDT batch")),
         what = factor(what, 
                       levels = c("sens", "spec"),
                       labels = c("sensitivity", "specificity")),
         test = test_dict[test],
         test = factor(test, levels = test_dict)) %>% 
  ggplot(aes(x = covar, y = mean, color = test)) +
  geom_hline(aes(yintercept = 1), lty = 2, color = "darkgray") +
  geom_point(position = pd, aes(pch = test)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), position = pd, width = 0) +
  facet_grid(what ~ covar_class, switch = "y", scales = "free_x", space = "free_x") +
  theme_bw() +
  scale_y_log10() +
  labs(x = NULL, y = 'Odds ratio', color = "Test", shape = "Test") +
  theme(
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.border = element_blank(),
    axis.line.x = element_line()
  ) +
  scale_color_manual(values = colors_tests()) +
  guides(color = "none", shape = "none")

ggsave(p_eff, filename = "figures/sfigure_effect_sizes.png", width = 11, height = 5)



# Supfig S3: poststart by age ----------------------------------------------


ps_figures_age <- map(
  c("period", "antibioticsGTFCC", "time_culture"), 
  function(x) {
    
    ps <- poststrat[str_detect(names(poststrat), x)] %>% 
      bind_rows() %>% 
      mutate(age_cat = factor(age_cat, levels = c("all", "children", "adults")),
             test = factor(test, levels = c("RDT", "PCR", "culture"))) %>% 
      mutate(strata = case_when(strata == "[0,5)" ~ "0-4",
                                strata == "[5,Inf)" ~ "5+",
                                strata == "0" & str_detect(covariate, "antibiotics") ~ "No",
                                strata == "1" & str_detect(covariate, "antibiotics") ~ "Yes",
                                T ~ strata))
    
    if (x == "antibioticsGTFCC") {
      ps$covariate <- "antibiotic use"
    }
    
    # Add fake datapoint for specicity
    ps_fake <- tibble(mean = c(.5, 1), 
                      test = "culture", 
                      what = "specificity", 
                      covariate =  ps$covariate[1],
                      strata = ps$strata[1],
                      age_cat = "adults") %>% 
      bind_rows(tibble(mean = c(.35, 1), 
                       test = "culture", 
                       what = "sensitivity", 
                       covariate =  ps$covariate[1],
                       strata = ps$strata[1],
                       age_cat = "adults")) %>% 
      mutate(test = factor(test, levels = c("RDT", "PCR", "culture")))
    
    
    if (x == "period") {
      ps <- mutate(ps, strata = factor(strata, levels = c("pre-monsoon", "monsoon", "post-monsoon", "winter")))
    }
    
    if (ps$covariate[1] %in% c("time_culture")) {
      xlab <- "days"
    } else {
      xlab <- NULL
    }
    
    if (ps$covariate[1] %in% c("overall")) {
      ylab <- "Post-stratified estimate"
    } else {
      ylab <- NULL
    }
    
    pd <- position_dodge(width = .3)
    
    if (ps$covariate[1] != c("time_culture")) {
      p <- ps %>% 
        filter(!is.na(strata)) %>% 
        ggplot(aes(x = strata, y = mean, color = test, group = age_cat)) +
        geom_errorbar(position = pd, width = 0, aes(ymin = q5, ymax = q95), lwd = .3) +
        geom_point(position = pd, aes(pch = age_cat))   +
        geom_point(data = ps_fake, alpha = 0)   +
        guides(color = "none", shape = "none") +
        theme_bw()   +
        ggh4x::facet_nested(what ~ covariate + test, scales = "free", space = "free_x", drop = FALSE,
                            labeller = labeller(covariate = label_dict),
                            switch = "y") 
    } else {
      p <- ps %>% 
        filter(!is.na(strata)) %>%
        mutate(strata = as.numeric(as.character(strata))) %>%  
        ggplot(aes(x = strata, y = mean, group = age_cat)) +
        geom_ribbon(aes(ymin = q5, ymax = q95, fill = test), alpha = .1, show.legend = TRUE) +
        geom_point(data = ps_fake %>% filter(what == "sensitivity"), 
                   aes(x = 10), alpha = 0)   +
        geom_line(aes(color = test, lty = age_cat), show.legend = TRUE) +
        geom_rug(inherit.aes = F,
                 data = lab_data %>% 
                   filter(time_to_culture > 0),
                 aes(x = time_to_culture, y = .5),
                 position = "jitter",
                 sides = "b",
                 alpha = .3) +
        theme_bw() +
        theme(legend.position = "right",
              legend.direction = "vertical") +
        guides(linetype = guide_legend("age"))  +
        facet_grid(what ~ covariate, scales = "free", space = "free_x", drop = FALSE,
                   labeller = labeller(covariate = label_dict),
                   switch = "y") 
    }
    
    p <- p +
      labs(x = xlab, y = ylab) +
      # scale_y_continuous(limits = c(.3, 1)) +
      scale_color_manual(values = colors_tests(), drop = F, limits = c("RDT", "PCR", "culture")) +
      scale_fill_manual(values = colors_tests(), drop = F, limits = c("RDT", "PCR", "culture"))
    
    if (!(ps$covariate[1] %in% c("period", "time_culture"))) {
      p <- p + theme(strip.text.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
    }
    
    if (x == "period") {
      p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      p <- p + guides(shape = guide_legend("age"))# +
      # theme(legend.position = c(.14, .25))
    }
    
    p 
    
  })

p_poststrat_age <- plot_grid(
  plotlist = ps_figures_age,
  nrow = 1,
  rel_widths = c(2, 1.4, 2),
  align = "h",
  axis = "tb"
) +
  theme(plot.margin = unit(c(0, 1, 0, 1.5), units = "lines"),
        plot.background = element_rect(fill = "white", color = "white"))

ggsave(p_poststrat_age, filename = "figures/supfig_poststrat_age.png", 
       width = 14, height = 5)


# Supfig S4: antibiotic definitions ---------------------------------------


ps_ab <- poststrat[str_detect(names(poststrat), "antibiotics")] %>% 
  bind_rows() %>% 
  mutate(age_cat = factor(age_cat, levels = c("all", "children", "adults")),
         test = factor(test, levels = c("RDT", "PCR", "culture"))) %>% 
  mutate(strata = case_when(strata == "[0,5)" ~ "0-4",
                            strata == "[5,Inf)" ~ "5+",
                            strata == "0" & str_detect(covariate, "antibiotics") ~ "No",
                            strata == "1" & str_detect(covariate, "antibiotics") ~ "Yes",
                            T ~ strata))

# Add fake datapoint for specicity
ps_fake <- tibble(mean = c(.5, 1), 
                  test = "culture", 
                  what = "specificity", 
                  covariate =  "antibiotics",
                  strata = ps_ab$strata[1],
                  age_cat = "adults") %>% 
  bind_rows(tibble(mean = c(.35, 1), 
                   test = "culture", 
                   what = "sensitivity", 
                   covariate =  "antibiotics",
                   strata = ps_ab$strata[1],
                   age_cat = "adults"))  %>% 
  mutate(age_cat = factor(age_cat, levels = c("all", "children", "adults")),
         test = factor(test, levels = c("RDT", "PCR", "culture")))

ab_dict <- c("antibiotics" = "any antibiotics",
             "antibioticsGTFCC" = "GTFFC recommendation",
             "antibioticsEff" = "efficient against V. cholerae")

pd <- position_dodge(width = .3)
p_poststrat_ab <- ps_ab %>% 
  filter(!is.na(strata)) %>% 
  mutate(covariate = ab_dict[covariate],
         covarite = factor(covariate, levels = ab_dict)) %>% 
  ggplot(aes(x = strata, y = mean, color = covariate, group = covariate)) +
  geom_errorbar(position = pd, width = 0, aes(ymin = q5, ymax = q95), lwd = .3) +
  geom_point(position = pd, aes(pch = age_cat))   +
  # geom_point(data = ps_fake, alpha = 0)   +
  ggh4x::facet_nested(what + test  ~  age_cat,
                      switch = "y", drop = T) +
  labs(x = "antibiotic use", y = "test performance",
       color = "definition of\nantibiotic use",
       shape = "age") +
  theme_bw()

ggsave(p_poststrat_ab, filename = "figures/supfig_poststrat_antibiotics.png", 
       width = 8, height = 7)


# Supfig S5: NPV/PPV for culture ------------------------------------------
# This is an ugly copy/paste, could make into function

# PPV and NPV
sens_spec_gen_culture <- genquant_age$draws("gen_culture_sens_poststrat") %>% 
  posterior::as_draws() %>% 
  posterior::as_draws_df() %>% 
  as_tibble() %>% 
  pivot_longer(cols = contains("post")) %>% 
  mutate(age_cat = c("0-4", "5+")[as.numeric(str_extract(name, "[1-2]"))]) %>% 
  select(-name) %>% 
  rename(sens = value) %>% 
  inner_join(
    genquant_age$draws("gen_culture_spec_poststrat") %>% 
      posterior::as_draws() %>% 
      posterior::as_draws_df() %>% 
      as_tibble() %>% 
      pivot_longer(cols = contains("post")) %>% 
      mutate(age_cat = c("0-4", "5+")[as.numeric(str_extract(name, "[1-2]"))]) %>% 
      select(-name) %>% 
      rename(spec = value)
  ) 

# Compute PPV and NPV at mean prevalences
pv_mean_culture <- sens_spec_gen_culture %>% 
  inner_join(prior_chol_median) %>% 
  rowwise() %>% 
  mutate(ppv = computePPV(p = prior_chol, sens = sens, spec = spec),
         npv = computeNPV(p = prior_chol, sens = sens, spec = spec)) %>% 
  pivot_longer(cols = c("ppv", "npv")) %>% 
  group_by(age_cat, epi_period, name) %>% 
  summarise(prior_chol = mean(prior_chol),
            mean = mean(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975)) %>% 
  ungroup() %>% 
  mutate(age_cat = factor(age_cat,
                          labels = c("0-4", "5+")),
         epi_period = factor(epi_period, levels = c("pre-monsoon", "monsoon", "post-monsoon", "winter"))) 

pv_pred_culture <- sens_spec_gen_culture %>% 
  inner_join(
    expand.grid(
      age_cat = unique(sens_spec_gen_culture$age_cat),
      prior_chol = seq(1e-7,  max(lab_data$prior_chol), by = .01)
    )
  ) %>% 
  mutate(ppv = computePPV(p = prior_chol, sens = sens, spec = spec),
         npv = computeNPV(p = prior_chol, sens = sens, spec = spec)) %>% 
  pivot_longer(cols = c("ppv", "npv")) %>% 
  group_by(age_cat, name, prior_chol) %>% 
  summarise(prior_chol = mean(prior_chol),
            mean_value = mean(value),
            q025 = quantile(value, 0.025),
            q975 = quantile(value, 0.975)) %>% 
  ungroup() 


p_supfig_pv_culture <- pv_pred_culture %>% 
  ggplot(aes(x = prior_chol)) +
  geom_rect(data = prior_chol_iqr_age, inherit.aes = F, 
            aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = epi_period),
            alpha = .15) +
  geom_vline(data = prior_chol_iqr_age, 
             aes(xintercept = median, lty = epi_period),
             lwd = .7, color = "white", alpha = .3) +
  geom_vline(data = prior_chol_iqr_age, 
             aes(xintercept = median, color = epi_period, lty = epi_period),
             lwd = .2) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .3) +
  geom_line(aes(y = mean_value), col = "black") +
  geom_rug(data = prior_chol_all, alpha = .1, aes(color = epi_period)) +
  geom_errorbar(data = pv_mean_culture, aes(ymin = lo, ymax = hi, color = epi_period), 
                # col = "red", 
                width = .01) +
  geom_point(data = pv_mean_culture, aes(y = mean, pch = epi_period, color = epi_period), #col = "red", #pch = 21, 
             fill = "white", size = 2.5) +
  facet_grid(name ~ age_cat, scales = "free",
             labeller = labeller(name = c("npv" = "Negative Predictive value",
                                          "ppv" = "Positive Predictive value"))) +
  theme_bw() +
  labs(x = "Cholera prevalence among AWD cases",
       y = "Negative/Positive predictive value",
       color = "Season",
       linetype = "Season",
       shape = "Season",
       fill = "Season") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_color_manual(values = colors_seasons()) +
  scale_fill_manual(values = colors_seasons())


ggsave(p_supfig_pv_culture, filename = "figures/supfig_npv_ppv_culture.png", 
       width = 10, height = 6, dpi = 300)


# Supfig S6: Posterior retrodictive checks --------------------------------

lab_data <- readRDS("generated_data/lab_data.rds") %>%
  mutate(set = case_when(is.na(cult_res) & is.na(pcr_res) ~ "A",
                         is.na(cult_res) & !is.na(pcr_res) ~ "B",
                         T ~ "C"))

# Compute fraction for which PCR is available
frac_B <- computeFracB(lab_data) %>%
  mutate(month_date = str_c(yearnum, monthnum, "01", sep = "-") %>% as.Date())

# Generated samples for model with mean and all age cats
genquant <- readRDS("generated_data/genquant_age.rds")

# Compare counts in each sampling case
# A: only RdT is available
# B: RDT and PCR available
# C: RDT, PCR and culture available

counts_A <- lab_data %>%
  filter(set == "A") %>%
  # Which test differs from the others?
  mutate(diff = case_when(T ~ "n_A"),
         label = "<RDT:-, PCR:NA, cult:NA>") %>%
  count(month, diff, label) %>%
  inner_join(frac_B, .)

counts_B <- lab_data %>%
  filter(set == "B") %>%
  # Which test differs from the others?
  mutate(diff = case_when(!rdt_res & pcr_res ~ "rdt_FN",
                          T ~ "concordant"),
         label = case_when(!rdt_res & pcr_res ~ "<RDT:-, PCR:+, cult:NA>",
                           T ~ "<RDT:-, PCR:-, cult:NA>")) %>%
  count(month, diff, label) %>%
  inner_join(frac_B, .)

counts_C <- lab_data %>%
  filter(set == "C") %>%
  # Which test differs from the others?
  mutate(diff = case_when(!pcr_res & !cult_res ~ "rdt_FP",
                          !pcr_res & cult_res ~ "pcr_FN",
                          pcr_res & !cult_res ~ "culture_FN",
                          T ~ "concordant"),
         label = case_when(!pcr_res & !cult_res ~ "<RDT:+, PCR:-, cult:->",
                           !pcr_res & cult_res ~ "<RDT:+, PCR:-, cult:+>",
                           pcr_res & !cult_res ~ "<RDT:+, PCR:+, cult:->",
                           T ~ "<RDT:+, PCR:+, cult:+>"))  %>%
  count(month, diff, label) %>%
  inner_join(frac_B %>% select(month, month_date, yearnum, monthnum), .)


sim_counts_A <- genquant$summary(c("n_A"), custom_summaries()) %>%
  mutate(month = str_extract(variable, "[0-9]+") %>% as.numeric(),
         month_date = frac_B$month_date[month],
         month = frac_B$month[month],
         diff = str_extract(variable, "(.)*(?=\\[)")) %>%
  left_join(counts_A) %>%
  replace_na(list(n = 0))

sim_counts_B <- genquant$summary(c("n_B_concordant", "n_B_rdt_FN"), custom_summaries()) %>%
  mutate(month = str_extract(variable, "[0-9]+") %>% as.numeric(),
         month_date = frac_B$month_date[month],
         month = frac_B$month[month],
         diff = str_extract(variable, "(.)*(?=\\[)") %>% str_remove("n_B_")) %>%
  left_join(counts_B) %>%
  replace_na(list(n = 0)) %>% 
  mutate(label = case_when(diff == "concordant" ~ "<RDT:-, PCR:-, cult:NA>",
                           T ~ "<RDT:-, PCR:+, cult:NA>"))

sim_counts_C <- genquant$summary(c("n_C_concordant", "n_C_rdt_FP",
                                   "n_C_pcr_FN", "n_C_culture_FN"), custom_summaries()) %>%
  mutate(month = str_extract(variable, "[0-9]+") %>% as.numeric(),
         month_date = frac_B$month_date[month],
         month = frac_B$month[month],
         diff = str_extract(variable, "(.)*(?=\\[)") %>% str_remove("n_C_")) %>%
  left_join(counts_C) %>%
  replace_na(list(n = 0)) %>% 
  mutate(label = case_when(diff == "concordant" ~ "<RDT:+, PCR:+, cult:+>",
                           diff == "rdt_FP" ~ "<RDT:+, PCR:-, cult:->",
                           diff == "pcr_FN" ~ "<RDT:+, PCR:-, cult:+>",
                           T ~ "<RDT:-, PCR:+, cult:NA>"))

p_pp <- cowplot::plot_grid(
  plot_pp(sim_counts_A, simple = TRUE),
  plot_pp(sim_counts_B, simple = TRUE),
  plot_pp(sim_counts_C),
  ncol = 1,
  rel_heights = c(1, 2, 4),
  align = "v",
  axis = "lr",
  labels = "auto"
)

ggsave(p_pp, filename = "figures/supfig_posterior_checks.png",
       width = 10, 
       height = 9)
