pacman::p_load(
  haven, tidyverse, labelled, janitor, lubridate,
  survey, brms, lme4, srvyr, loo, bayesplot, posterior,
  tidybayes, ggdist, patchwork, marginaleffects, gt
)

#Install cmdstar if you do not have for bayesian analysis 
#remotes::install_github("stan-dev/cmdstanr")
library(cmdstanr)
#install_cmdstan(cores = 2)

# Read data
data <- read_dta("2015 STEPS Data.dta") 

# Data preparation
data1 <- data |>
  select(
    PID,
    "Sex" = C1,
    "Age" = C3,
    "Marital_status" = C7,
    "Residence" = X2a,
    "Education" = C5,
    "htn_diagnosed" = H2a,
    "htn_diagnosed_12m" = H2b,
    "htn_treated_2w" = H3,
    "systolic_bp_1" = M4a,
    "diastolic_bp_1" = M4b,
    "systolic_bp_2" = M5a,
    "diastolic_bp_2" = M5b,
    "systolic_bp_3" = M6a,
    "diastolic_bp_3" = M6b,
    "diabetes_diagnosed" = H7a,
    "salt_added_at_table" = D5,
    "salt_added_during_cooking" = D6,
    "processed_high_salt_food_freq" = D7,
    "diabetes_med_2w" = H8,
    "insulin_use" = H9,
    "tobaco_current" = T1,
    "smokeless_current" = T6,
    "vigorous_work_days" = P2, 
    "vigorous_work_min" = P3b,
    "moderate_work_days" = P5, 
    "moderate_work_min" = P6b,
    "transport_days" = P8, 
    "transport_min" = P9b,
    "vigorous_rec_days" = P11, 
    "vigorous_rec_min" = P12b,
    "moderate_rec_days" = P14,
    "moderate_rec_min" = P15b,
    "type_of_oil" = D12,
    "eat_outside" = D13, 
    "day_of_interview" = I4a,
    "month_of_interview" = I4b,
    "year_of_interview" = I4c,
    "day_of_birth" = C2a,
    "month_of_birth" = C2b,
    "year_of_birth" = C2c,
    "height" = M11,
    "weight" = M12,
    "alcohol_ever" = A1,
    "salt_added_freq" = D5,
    "cluster" = I1,
    "household_id" = I2,
    starts_with("A"),
    wealth, psu, stratum, wstep2) |>
  mutate(
    # Hypertension
    systolic  = (systolic_bp_2 + systolic_bp_3) / 2,
    diastolic = (diastolic_bp_2 + diastolic_bp_3) / 2,
    hypertension = case_when(
      systolic >= 140 | diastolic >= 90 | 
        htn_diagnosed == 1 | 
        htn_treated_2w == 1 ~ 1, 
      .default = 0),
    
    # For sensitivity analysis
    htn_measured_only = ifelse(systolic >= 140 | diastolic >= 90, 1, 0),
    htn_diagnosed_only = ifelse(htn_diagnosed == 1, 1, 0),
    htn_treated_only = ifelse(htn_treated_2w == 1, 1, 0),
    
    # Age
    birth_date = make_date(year_of_birth, month_of_birth, day_of_birth),
    interview_date = make_date(year_of_interview, month_of_interview, day_of_interview),
    age_computed = floor(time_length(interval(birth_date, interview_date), "years")),
    Age_s = case_when(
      !is.na(Age) ~ Age,
      !is.na(age_computed) & age_computed >= 18 ~ age_computed,
      TRUE ~ NA_real_),
    Age_s = cut(Age_s, 
                breaks = c(18, 25, 30, 35, 40, 45, 50, 60, 70),
                labels = c("18-24", "25-29", "30-34", "35-39", 
                           "40-44", "45-49", "50-59", "60-70"),
                right = FALSE),
    
    # BMI
    bmi = weight / (height/100)^2,
    bmi_cat = case_when(
      bmi < 18.5 ~ "Underweight",
      bmi >= 18.5 & bmi < 25 ~ "Normal",
      bmi >= 25 & bmi < 30 ~ "Overweight",
      bmi >= 30 ~ "Obese",
      TRUE ~ NA_character_
    ),
    bmi_cat = factor(bmi_cat, levels = c("Normal",
                                         "Underweight",
                                         "Overweight", 
                                         "Obese")),
    
    # Diabetes
    has_diabetes = case_when(
      diabetes_diagnosed == 1 | diabetes_med_2w == 1 |
        insulin_use == 1 ~ "Yes",
      TRUE ~ "No"
    ),
    
    # Salt intake
    high_salt_intake = case_when(
      salt_added_at_table %in% c(1,2) | 
        salt_added_during_cooking %in% c(1,2) |
        processed_high_salt_food_freq %in% c(1,2) ~ "Yes",
      TRUE ~ "No"
    ),
    
    # Tobacco
    current_smoker = case_when(
      tobaco_current == 1 | smokeless_current == 1 ~ "Yes",
      TRUE ~ "No"
    ),
    
    # Alcohol
    weekly_alcohol_consumption = rowSums(
      across(A10a:A10g, ~replace_na(., 0))) +
      coalesce(A12a, 0) + coalesce(A12b, 0) + 
      coalesce(A12c, 0) + coalesce(A12d, 0) + coalesce(A12e, 0),
    
    heavy_episodic_drinking = ifelse(A8 >= 1, 1, 0),
    heavy_episodic_drinking = replace_na(heavy_episodic_drinking, 0),
    
    high_weekly_intake = case_when(
      Sex == "Male" & weekly_alcohol_consumption >= 14 ~ 1,
      Sex == "Female" & weekly_alcohol_consumption >= 7  ~ 1,
      TRUE ~ 0
    ),
    
    high_alcohol_intake = case_when(
      heavy_episodic_drinking == 1 | high_weekly_intake == 1 ~ "Yes",
      TRUE ~ "No"
    ),
    
    # Fat intake
    type_of_oil = ifelse(type_of_oil == 77, NA, type_of_oil),
    eat_outside = ifelse(eat_outside == 77, NA, eat_outside),
    high_saturated_fat_oil = ifelse(type_of_oil %in% c(3,4,6,7), 1, 0),
    high_outside_meals = ifelse(eat_outside >= 3, 1, 0),
    high_fat_intake = case_when(
      high_saturated_fat_oil == 1 | high_outside_meals == 1 ~ "Yes",
      TRUE ~ "No"
    ),
    
    # Physical Activity
    across(c(vigorous_work_days, vigorous_work_min,
             moderate_work_days, moderate_work_min,
             transport_days, transport_min,
             vigorous_rec_days, vigorous_rec_min,
             moderate_rec_days, moderate_rec_min),
           ~ replace_na(.x, 0)),
    vigorous_minutes =
      (vigorous_work_days * vigorous_work_min) +
      (vigorous_rec_days * vigorous_rec_min),
    moderate_minutes =
      (moderate_work_days * moderate_work_min) +
      (moderate_rec_days * moderate_rec_min) +
      (transport_days * transport_min),
    total_met_min =
      (vigorous_minutes * 8) +
      (moderate_minutes * 4),
    insufficient_activity = case_when(
      vigorous_minutes >= 75 ~ "No",
      moderate_minutes >= 150 ~ "No",
      (moderate_minutes + 2 * vigorous_minutes) >= 150 ~ "No",
      TRUE ~ "Yes"
    ),
    
    #Marital status
    Marital_status = ifelse(Marital_status == "", NA, Marital_status),
    
    #Residence
    Residence = factor(Residence, labels = c("Rural", "Urban"))
  ) |>
  filter(!is.na(Marital_status)) |> 
  filter(wstep2 > 0) |> 
  mutate(across(where(is.character), factor)) |>
  select(PID, Age_s, Sex, Education, wealth, Residence, 
         Marital_status, bmi_cat, high_salt_intake,
         current_smoker, has_diabetes, high_alcohol_intake,
         insufficient_activity, high_fat_intake, 
         household_id, cluster, psu, stratum, wstep2,
         hypertension, htn_measured_only, htn_diagnosed_only, htn_treated_only)

#Survey design
data1_svy <- svydesign(
  ids = ~psu,
  strata = ~stratum,
  weights = ~wstep2,
  data = data1,
  nest = TRUE
)

# ----------------------------------------------------
# Objective 1
# ----------------------------------------------------
#To estimate population-averaged RRs for individual hypertension risk factors using
#logistic and Bayesian models

# Logistic Regression
model_freq_logistic <- svyglm(
  hypertension ~ Age_s + Sex + bmi_cat + has_diabetes + 
    current_smoker + high_alcohol_intake + high_fat_intake + 
    high_salt_intake + insufficient_activity +
    Education + wealth + Residence,
  design = data1_svy,
  family = quasibinomial(link = "logit")
)

# Modified Poisson
model_freq_poisson <- svyglm(
  hypertension ~ Age_s + Sex + bmi_cat + has_diabetes + 
    current_smoker + high_alcohol_intake + high_fat_intake + 
    high_salt_intake  + insufficient_activity +
    Education + wealth + Residence,
  design = data1_svy,
  family = quasipoisson(link = "log")
)

# BAYESIAN MODEL
data_unweighted <- data1 |>
  drop_na(hypertension, Age_s, Sex, bmi_cat, has_diabetes,
          current_smoker, high_alcohol_intake, high_fat_intake,
          high_salt_intake, insufficient_activity, Education, 
          wealth, Residence)

data_bayes <- data_unweighted |>
  mutate(
    weight_norm = wstep2 * (n() / sum(wstep2, na.rm = TRUE)), # Normalize weights
    weight_norm = ifelse(is.na(weight_norm), 1, weight_norm))

prior_single <- c(
  prior(normal(0, 2.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"))

model_bayes_single <- brm(
  hypertension | weights(weight_norm) ~ 
    Age_s + Sex + bmi_cat + has_diabetes + 
    current_smoker + high_alcohol_intake + high_fat_intake + 
    high_salt_intake + insufficient_activity +
    Education + wealth + Residence,
  data = data_bayes,
  family = bernoulli(link = "logit"),
  prior = prior_single,
  chains = 4, iter = 4000, warmup = 2000,
  seed = 123,
  backend = "cmdstanr",
  threads = threading(2)
)

# Marginal Risk Ratio using G computation
# Logistic Regression
freq_logistic_RRs <- avg_comparisons(
  model_freq_logistic,
  comparison = "ratio",
  conf_level = 0.95,
  wts = model_freq_logistic$weights
) |>
  as.data.frame() |>
  mutate(Model = "Freq Logistic (Marginal RR)")

# Poisson model
poisson_RRs <- data.frame(
  term = names(coef(model_freq_poisson))[-1],
  estimate = exp(coef(model_freq_poisson))[-1],
  conf.low = exp(confint(model_freq_poisson))[-1, 1],
  conf.high = exp(confint(model_freq_poisson))[-1, 2]
) |> 
  mutate(Model = "Modified Poisson (RR)")

# Bayesian Model
bayes_RRs <- avg_comparisons(
  model_bayes_single,
  comparison = "ratio", 
  conf_level = 0.95,
  wts = data_bayes$weight_norm,
  ndraws = 1000
) |> 
  as.data.frame() |> 
  select(term, contrast, estimate, conf.low, conf.high) |> 
  mutate(Model = "Bayesian Logistic (Marginal RR)")

clean_poisson_terms <- \(df) {
  df |>
    mutate(
      Variable = case_when(
        str_detect(term, "Age_s") ~ "Age_s",
        str_detect(term, "Sex") ~ "Sex",
        str_detect(term, "bmi_cat") ~ "bmi_cat",
        str_detect(term, "Education") ~ "Education",
        str_detect(term, "wealth") ~ "wealth",
        str_detect(term, "Residence") ~ "Residence",
        str_detect(term, "diabetes") ~ "has_diabetes",
        str_detect(term, "smoker") ~ "current_smoker",
        str_detect(term, "alcohol") ~ "high_alcohol_intake",
        str_detect(term, "fat") ~ "high_fat_intake",
        str_detect(term, "salt") ~ "high_salt_intake",
        TRUE ~ term
      ),
      Level = str_remove(term, Variable)
    )
}

gcomp_clean <- bind_rows(bayes_RRs, freq_logistic_RRs) |>
  mutate(
    Variable = term,
    Level = str_extract(contrast, "(?<=mean\\().+?(?=\\))") 
  )

poisson_clean <- poisson_RRs |>
  clean_poisson_terms()

table_data <- bind_rows(gcomp_clean, poisson_clean) |>
  mutate(
    Estimate_CI = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
    Model_Short = case_when(
      str_detect(Model, "Bayesian") ~ "Bayesian RR",
      str_detect(Model, "Freq Logistic") ~ "Logistic RR",
      str_detect(Model, "Poisson") ~ "Poisson RR"
    )
  ) |>
  select(Variable, Level, Model_Short, Estimate_CI)

table_objective_1 <- table_data |>
  pivot_wider(
    names_from = Model_Short,
    values_from = Estimate_CI
  ) |>
  group_by(Variable) |>
  gt() |>
  tab_header(
    title = "Comparative Analysis of Hypertension Risk Factors",
    subtitle = "Comparison of Risk Ratios (95% CI) across Frequentist and Bayesian Frameworks"
  ) |>
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_row_groups()
  ) |>
  tab_footnote(
    footnote = "Estimates derived using Design-Weighted G-Computation.",
    locations = cells_column_labels(columns = contains("Logistic"))
  ) |>
  tab_footnote(
    footnote = "Estimates derived using Pseudo-Posterior Weighted G-Computation.",
    locations = cells_column_labels(columns = contains("Bayesian"))
  ) |>
  sub_missing(missing_text = "-") 

table_objective_1

# -----------------------------------------------------------
# Objective 2
# -----------------------------------------------------------
# To estimate cluster-specific RRs for household and community influences using
# multilevel logistic and Bayesian models

#Individuals and household - Frequentist method
model_freq_2level <- glmer(
  hypertension ~ Age_s + Sex + bmi_cat + has_diabetes + 
    current_smoker + high_alcohol_intake + high_fat_intake + 
    high_salt_intake + Education + wealth + Residence + insufficient_activity +
    (1 | household_id),
  data = data_unweighted,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa", 
                         optCtrl = list(maxfun = 2e5))
)

summary(model_freq_2level)

# ICC at household level
variance_components <- as.data.frame(VarCorr(model_freq_2level))
sigma2_hh <- variance_components$vcov[1]
icc_household <- sigma2_hh / (sigma2_hh + pi^2/3)
icc_household

#Individual, hoursehold and PSU - Frequentist method
model_freq_3level <- glmer(
  hypertension ~ Age_s + Sex + bmi_cat + has_diabetes + 
    current_smoker + high_alcohol_intake + high_fat_intake + 
    high_salt_intake + insufficient_activity +
    Education + wealth + Residence + 
    (1 | psu/household_id),
  data = data_unweighted,
  family = binomial(),
  control = glmerControl(optimizer = "bobyqa", 
                         optCtrl = list(maxfun = 2e5))
)

summary(model_freq_3level)

# Calculate ICCs
vc_3level <- as.data.frame(VarCorr(model_freq_3level))
sigma2_psu <- vc_3level$vcov[1]
sigma2_hh_3level <- vc_3level$vcov[2]
total_var <- sigma2_psu + sigma2_hh_3level + pi^2/3

icc_psu <- sigma2_psu / total_var
icc_hh_within_psu <- sigma2_hh_3level / total_var

icc_psu
icc_hh_within_psu

# Predict including random effects for each individual
cluster_specific_preds_freq <- data_unweighted |>
  mutate(
    pred_cluster_specific = predict(model_freq_3level, type = "response"), # with random effects 
    pred_population_avg = predict(model_freq_3level, re.form = NA, type = "response") # without random effects
  )

# Summary by household
household_risks_freq <- cluster_specific_preds_freq |>
  group_by(household_id, psu) |>
  summarise(
    n_members = n(),
    mean_cluster_risk = mean(pred_cluster_specific),
    mean_pop_risk = mean(pred_population_avg),
    difference = mean_cluster_risk - mean_pop_risk,
    .groups = "drop"
  )
summary(household_risks_freq$mean_cluster_risk)

# Summary by PSU
psu_risks_freq <- cluster_specific_preds_freq |>
  group_by(psu) |>
  summarise(
    n_households = n_distinct(household_id),
    n_individuals = n(),
    mean_cluster_risk = mean(pred_cluster_specific),
    mean_pop_risk = mean(pred_population_avg),
    difference = mean_cluster_risk - mean_pop_risk,
    .groups = "drop"
  )
summary(psu_risks_freq$mean_cluster_risk)

#Bayesian
# Individual and household
prior_2level <- c(
  prior(normal(0, 2.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(cauchy(0, 1), class = "sd")  
)

model_bayes_2level <- brm(
  hypertension | weights(weight_norm) ~ 
    Age_s + Sex + bmi_cat + has_diabetes + 
    current_smoker + high_alcohol_intake + high_fat_intake + 
    high_salt_intake + insufficient_activity + 
    Education + wealth + Residence +
    (1 | household_id),
  data = data_bayes,
  family = bernoulli(),
  prior = prior_2level,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  seed = 123,
  backend = "cmdstanr"
)
summary(model_bayes_2level)

#ICC - Bayesian
bayes_2level_icc <- model_bayes_2level |>
  spread_draws(sd_household_id__Intercept) |>
  mutate(
    sigma2_hh = sd_household_id__Intercept^2,
    icc = sigma2_hh / (sigma2_hh + pi^2/3)
  ) |>
  median_qi(icc, .width = 0.95)
bayes_2level_icc

#Individual, household and psu
prior_3level <- c(
  prior(normal(0, 2.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(cauchy(0, 1), class = "sd", group = "psu"),
  prior(cauchy(0, 1), class = "sd", group = "psu:household_id")
)

model_bayes_3level <- brm(
  hypertension | weights(weight_norm) ~ 
    Age_s + Sex + bmi_cat + has_diabetes + 
    current_smoker + high_alcohol_intake + high_fat_intake + 
    high_salt_intake + insufficient_activity +
    Education + wealth + Residence +
    (1 | psu/household_id),
  data = data_bayes,
  family = bernoulli(),
  prior = prior_3level,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  control = list(adapt_delta = 0.95),
  seed = 123,
  backend = "cmdstanr"
)
summary(model_bayes_3level)

# Cluster-specific predictions with random effects
cluster_preds_bayes <- posterior_epred(
  model_bayes_3level,
  newdata = data_bayes,
  ndraws = 500,
  re_formula = NULL  
)

# Population-averaged predictions without random effects
pop_preds_bayes <- posterior_epred(
  model_bayes_3level,
  newdata = data_bayes,
  ndraws = 500,
  re_formula = NA  
)

# Summary at household level 
household_risks_bayes <- data_bayes |>
  mutate(
    median_cluster_risk = apply(cluster_preds_bayes, 2, median),
    lower_cluster_risk = apply(cluster_preds_bayes, 2, quantile, 0.025),
    upper_cluster_risk = apply(cluster_preds_bayes, 2, quantile, 0.975),
    median_pop_risk = apply(pop_preds_bayes, 2, median)
  ) |>
  group_by(household_id, psu) |>
  summarise(
    n_members = n(),
    mean_cluster_risk = mean(median_cluster_risk),
    lower_cluster_risk = mean(lower_cluster_risk),
    upper_cluster_risk = mean(upper_cluster_risk),
    mean_pop_risk = mean(median_pop_risk),
    .groups = "drop"
  )
summary(household_risks_bayes$mean_cluster_risk)

# Cluster specific relative risk
calc_cluster_specific_rr <- \(model, data, exposure_var, 
                              reference_level, 
                              is_bayesian = FALSE) {
  
  exposure_levels <- unique(data[[exposure_var]])
  exposure_levels <- c(reference_level, 
                       setdiff(exposure_levels, reference_level))
  
  cluster_rr_results <- list()
  
  for (cluster in unique(data$psu)) {
    cluster_data <- data |> filter(psu == cluster)
    
    if (nrow(cluster_data) < 5) next  
    
    cluster_risks <- list()
    
    for (level in exposure_levels) {
      cluster_data_cf <- cluster_data
      cluster_data_cf[[exposure_var]] <- level
      
      if (is_bayesian) {
        preds <- posterior_epred(model, newdata = cluster_data_cf, 
                                 ndraws = 500, re_formula = NULL)
        cluster_risk <- mean(colMeans(preds))
      } else {
        preds <- predict(model, newdata = cluster_data_cf, type = "response")
        cluster_risk <- mean(preds)
      }
      
      cluster_risks[[as.character(level)]] <- cluster_risk
    }
    
    ref_risk <- cluster_risks[[1]]
    cluster_rrs <- sapply(cluster_risks, function(x) x / ref_risk)
    
    cluster_rr_results[[as.character(cluster)]] <- tibble(
      psu = cluster,
      exposure_level = names(cluster_rrs),
      cluster_rr = cluster_rrs
    )
  }
  
  return(bind_rows(cluster_rr_results))
}

# Cluster specific RRs for diabetes 
cluster_rr_diabetes_bayes <- calc_cluster_specific_rr(
  model_bayes_3level,
  data_bayes,
  exposure_var = "has_diabetes",
  reference_level = "No",
  is_bayesian = TRUE
)

# Summary of cluster-specific variation
cluster_rr_summary <- cluster_rr_diabetes_bayes |>
  filter(exposure_level == "Yes") |> 
  summarise(
    n_clusters = n(),
    mean_rr = mean(cluster_rr),
    median_rr = median(cluster_rr),
    sd_rr = sd(cluster_rr),
    min_rr = min(cluster_rr),
    max_rr = max(cluster_rr),
    q25 = quantile(cluster_rr, 0.25),
    q75 = quantile(cluster_rr, 0.75)
  )
cluster_rr_summary

#Population average vs cluster specific average
pop_avg_diabetes_rr <- bayes_RRs |>
  filter(str_detect(term, "diabetes")) |>
  pull(estimate)

# Cluster-specific median
cluster_median_rr <- cluster_rr_summary$median_rr

comparison_table <- tibble(
  Estimand = c("Population-Averaged RR", "Median Cluster-Specific RR", 
               "Difference", "% Difference"),
  Diabetes = c(
    sprintf("%.2f", pop_avg_diabetes_rr),
    sprintf("%.2f", cluster_median_rr),
    sprintf("%.2f", pop_avg_diabetes_rr - cluster_median_rr),
    sprintf("%.1f%%", 
            (pop_avg_diabetes_rr - cluster_median_rr) / pop_avg_diabetes_rr * 100)
  )
)
comparison_table

# Cluster-specific comparisons(only diabetes)
cluster_comparisons_bayes <- avg_comparisons(
  model_bayes_3level,
  variables = list(has_diabetes = c("No", "Yes")),
  comparison = "ratio",
  by = "psu", 
  wts = data_bayes$weight_norm,
  ndraws = 500
) |>
  as.data.frame()

# Summary of variation in cluster-specific effects
cluster_variation_summary <- cluster_comparisons_bayes |>
  summarise(
    n_clusters = n(),
    mean_rr = mean(estimate),
    median_rr = median(estimate),
    sd_rr = sd(estimate),
    cv = sd_rr / mean_rr, 
    iqr = IQR(estimate)
  )
cluster_variation_summary

bayes_3level_vars <- model_bayes_3level |>
  spread_draws(sd_psu__Intercept, `sd_psu:household_id__Intercept`) |>
  mutate(
    var_psu = sd_psu__Intercept^2,
    var_hh = `sd_psu:household_id__Intercept`^2,
    var_resid = 3.29, 
    icc_psu = var_psu / (var_psu + var_hh + var_resid),
    icc_hh = var_hh / (var_psu + var_hh + var_resid)
  ) |>
  median_qi(icc_psu, icc_hh, .width = 0.95)

table_objective_2 <- tibble(
  Parameter = c(
    "ICC (PSU level) - Frequentist",
    "ICC (Household within PSU) - Frequentist",
    "ICC (PSU level) - Bayesian",
    "ICC (Household within PSU) - Bayesian",
    "",
    "Population-Avg RR (Diabetes)",
    "Median Cluster-Specific RR (Diabetes)",
    "Range of Cluster RRs",
    "SD of Cluster RRs",
    "Coefficient of Variation"
  ),
  Estimate = c(
    sprintf("%.4f", icc_psu),
    sprintf("%.4f", icc_hh_within_psu),
    sprintf("%.4f [%.4f, %.4f]", 
            bayes_2level_icc$icc,
            bayes_2level_icc$.lower,
            bayes_2level_icc$.upper),
    sprintf("%.4f [%.4f, %.4f]", 
            bayes_3level_vars$icc_psu, 
            bayes_3level_vars$icc_psu.lower, 
            bayes_3level_vars$icc_psu.upper),  # 
    "",
    sprintf("%.2f", pop_avg_diabetes_rr),
    sprintf("%.2f", cluster_rr_summary$median_rr),
    sprintf("%.2f - %.2f", 
            cluster_rr_summary$min_rr, 
            cluster_rr_summary$max_rr),
    sprintf("%.2f", cluster_rr_summary$sd_rr),
    sprintf("%.2f", cluster_variation_summary$cv)
  )
) |>
  gt() |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = Parameter == "")
  ) |>
  tab_footnote(
    footnote = "ICC quantifies the proportion of variance at each level",
    locations = cells_body(columns = Parameter, rows = 1:4)
  )

table_objective_2


# -----------------------------------------------------------
# Objective 3
# -----------------------------------------------------------
# To evaluate the impact of outcome misclassification of hypertension status through
# probabilistic bias analysis.

misclass_components <- data1 |>
  summarise(
    n = n(),
    composite_htn = sum(hypertension == 1, na.rm = TRUE),
    measured_only = sum(htn_measured_only == 1, na.rm = TRUE),
    diagnosed_only = sum(htn_diagnosed_only == 1, na.rm = TRUE),
    treated_only = sum(htn_treated_only == 1, na.rm = TRUE)
  ) |>
  mutate(
    composite_prev = composite_htn / n * 100,
    measured_prev = measured_only / n * 100,
    diagnosed_prev = diagnosed_only / n * 100,
    treated_prev = treated_only / n * 100
  )
misclass_components

cross_tab <- data1 |>
  mutate(
    detection_method = case_when(
      htn_measured_only == 1 & htn_diagnosed_only == 0 & htn_treated_only == 0 ~ "Measured only",
      htn_diagnosed_only == 1 & htn_measured_only == 0 ~ "Diagnosed only",
      htn_treated_only == 1 & htn_measured_only == 0 ~ "Treated only",
      htn_measured_only == 1 & htn_diagnosed_only == 1 ~ "Measured + Diagnosed",
      htn_measured_only == 1 & htn_treated_only == 1 ~ "Measured + Treated",
      hypertension == 1 ~ "Multiple components",
      TRUE ~ "No hypertension"
    )
  ) |>
  count(detection_method) |>
  mutate(pct = n / sum(n) * 100)
cross_tab

# 3.2 Sensitivity Analysis
outcomes_to_test <- c("hypertension", "htn_measured_only", 
                      "htn_diagnosed_only")

sensitivity_results <- list()

for (outcome in outcomes_to_test) {
  cat("Fitting model with outcome:", outcome, "\n")
  
  data_sens <- data_bayes |>
    mutate(outcome_var = .data[[outcome]])
  
  model_sens <- brm(
    outcome_var | weights(weight_norm) ~ 
      Age_s + Sex + bmi_cat + has_diabetes + 
      current_smoker + high_alcohol_intake + high_fat_intake + 
      high_salt_intake + Education + wealth + 
      Residence + insufficient_activity +
      (1 | psu/household_id),
    data = data_sens,
    family = bernoulli(),
    prior = prior_3level,
    chains = 2,
    iter = 2000,
    warmup = 1000,
    seed = 123,
    backend = "cmdstanr",
    refresh = 0,
    file = paste0("model_sens_", outcome)
  )
  
  diabetes_coef <- fixef(model_sens)["has_diabetesYes", ]
  
  sensitivity_results[[outcome]] <- tibble(
    outcome_definition = outcome,
    diabetes_OR = exp(diabetes_coef["Estimate"]),
    diabetes_OR_lower = exp(diabetes_coef["Q2.5"]),
    diabetes_OR_upper = exp(diabetes_coef["Q97.5"])
  )
}

sensitivity_comparison <- bind_rows(sensitivity_results) |>
  mutate(
    outcome_label = case_when(
      outcome_definition == "hypertension" ~ "Composite (Primary)",
      outcome_definition == "htn_measured_only" ~ "Measured BP Only",
      outcome_definition == "htn_diagnosed_only" ~ "Prior Diagnosis Only"
    ),
    OR_CI = sprintf("%.2f (%.2f, %.2f)", 
                    diabetes_OR, diabetes_OR_lower, diabetes_OR_upper)
  )
sensitivity_comparison

# Probabilistic Bias Analysis
# Stan model for measurement error correction
stan_misclass_code <- "
data {
  int<lower=0> N;                      
  int<lower=0> K;                      
  matrix[N, K] X;                      
  array[N] int<lower=0,upper=1> y_obs; 
  vector[N] weights;                   
  
  // Prior parameters for sensitivity/specificity
  real<lower=0> sens_alpha;
  real<lower=0> sens_beta;
  real<lower=0> spec_alpha;
  real<lower=0> spec_beta;
}

parameters {
  vector[K] beta;                      // regression coefficients
  real<lower=0,upper=1> sensitivity;   // P(y_obs=1 | y_true=1)
  real<lower=0,upper=1> specificity;   // P(y_obs=0 | y_true=0)
}

transformed parameters {
  vector[N] pi_true;
  
  // Linear predictor (Log-Odds of True Disease)
  pi_true = inv_logit(X * beta); 
}

model {
  // Priors
  beta ~ normal(0, 1);                 
  sensitivity ~ beta(sens_alpha, sens_beta);
  specificity ~ beta(spec_alpha, spec_beta);
  
  // Likelihood accounting for misclassification AND WEIGHTS
  for (n in 1:N) {
    real p_obs_1;
    
    // Law of total probability: P(Obs=1) = P(True=1)*Se + P(True=0)*(1-Sp)
    p_obs_1 = pi_true[n] * sensitivity + (1 - pi_true[n]) * (1 - specificity);
    
    if (y_obs[n] == 1) {
      // Weighting the log-likelihood
      target += weights[n] * log(p_obs_1); 
    } else {
      target += weights[n] * log(1 - p_obs_1);
    }
  }
}

generated quantities {
  vector[N] y_true_prob;               // Predicted probability of TRUE disease
  real prev_observed;
  real prev_true;
  
  // Calculate prevalence (weighted)
  prev_observed = dot_product(to_vector(y_obs), weights) / sum(weights);
  prev_true = dot_product(pi_true, weights) / sum(weights);
  
  // Calculate posterior probability of true disease for each person
  for (n in 1:N) {
    real p_obs_1 = pi_true[n] * sensitivity + (1 - pi_true[n]) * (1 - specificity);
    
    if (y_obs[n] == 1) {
      y_true_prob[n] = (pi_true[n] * sensitivity) / p_obs_1;
    } else {
      y_true_prob[n] = (pi_true[n] * (1 - sensitivity)) / (1 - p_obs_1);
    }
  }
}
"
X_matrix <- model.matrix(
  ~ Age_s + Sex + bmi_cat + has_diabetes + 
    current_smoker + high_alcohol_intake + high_fat_intake + 
    high_salt_intake + Education + wealth + Residence + insufficient_activity,
  data = data_bayes
)

# Literature-informed priors for sensitivity/specificity
# Based on Kislaya et al. (2019) and Goncalves et al. (2018)
# Composite definition: Sensitivity ~80-90%, Specificity ~85-95%

stan_data <- list(
  N = nrow(X_matrix),
  K = ncol(X_matrix),
  X = X_matrix,
  y_obs = data_bayes$hypertension,
  weights = data_bayes$weight_norm,
  sens_alpha = 85,      
  sens_beta = 15,
  spec_alpha = 90,      
  spec_beta = 10
)

# Fitting Stan model
stan_file <- write_stan_file(stan_misclass_code)
model_compiled <- cmdstan_model(stan_file)

fit_misclass <- model_compiled$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  seed = 123,
  refresh = 500,
  adapt_delta = 0.95
)

# Checking convergence
fit_misclass$cmdstan_diagnose()

draws_misclass <- fit_misclass$draws(format = "df")

# Sensitivity and specificity posteriors
sens_spec_summary <- fit_misclass$summary(
  variables = c("sensitivity", "specificity", "prev_observed", "prev_true")
)
sens_spec_summary

# Sensitivity/specificity posteriors
sens_draws <- draws_misclass$sensitivity
spec_draws <- draws_misclass$specificity

p1 <- ggplot(data.frame(sensitivity = sens_draws), aes(x = sensitivity)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = mean(sens_draws), color = "red", linetype = "dashed") +
  labs(title = "Posterior: Sensitivity", x = "Sensitivity", y = "Density") +
  theme_minimal()

p2 <- ggplot(data.frame(specificity = spec_draws), aes(x = specificity)) +
  geom_density(fill = "darkgreen", alpha = 0.6) +
  geom_vline(xintercept = mean(spec_draws), color = "red", linetype = "dashed") +
  labs(title = "Posterior: Specificity", x = "Specificity", y = "Density") +
  theme_minimal()

p1 + p2

# Bias-corrected vs naive estimates
beta_names <- colnames(X_matrix)
bias_corrected_coefs <- fit_misclass$summary(variables = "beta")

naive_coefs <- fixef(model_bayes_3level)

comparison_bias <- tibble(
  Parameter = beta_names,
  Naive_Estimate = naive_coefs[, "Estimate"],
  BiasCorr_Estimate = bias_corrected_coefs$mean,
  Difference = BiasCorr_Estimate - Naive_Estimate,
  Pct_Change = (Difference / Naive_Estimate) * 100
)
comparison_bias

# Using only key parameters
key_params <- comparison_bias |>
  filter(str_detect(Parameter, "diabetes|bmi_cat|Age_s")) |>
  mutate(
    Naive_OR = exp(Naive_Estimate),
    BiasCorr_OR = exp(BiasCorr_Estimate),
    OR_Difference = BiasCorr_OR - Naive_OR
  ) |>
  select(Parameter, Naive_OR, BiasCorr_OR, OR_Difference)
key_params

# Grid of sensitivity/specificity values
sensitivity_grid <- expand_grid(
  sensitivity = seq(0.70, 0.95, by = 0.05),
  specificity = seq(0.80, 0.98, by = 0.03)
)

sens_analysis_results <- list()

for (i in 1:nrow(sensitivity_grid)) {
  if (i %% 10 == 0) cat("  Scenario", i, "of", nrow(sensitivity_grid), "\n")
  
  sens_val <- sensitivity_grid$sensitivity[i]
  spec_val <- sensitivity_grid$specificity[i]
  
  sens_alpha <- sens_val * 100
  sens_beta <- (1 - sens_val) * 100
  spec_alpha <- spec_val * 100
  spec_beta <- (1 - spec_val) * 100
  
  stan_data_sens <- stan_data
  stan_data_sens$sens_alpha <- sens_alpha
  stan_data_sens$sens_beta <- sens_beta
  stan_data_sens$spec_alpha <- spec_alpha
  stan_data_sens$spec_beta <- spec_beta
  
  fit_sens <- model_compiled$sample(
    data = stan_data_sens,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 123 + i,
    refresh = 0,
    show_messages = FALSE
  )
  
  # Extracting diabetes coefficient 
  diabetes_idx <- which(beta_names == "has_diabetesYes")
  diabetes_beta <- fit_sens$summary(variables = paste0("beta[", diabetes_idx, "]"))
  
  sens_analysis_results[[i]] <- tibble(
    sensitivity = sens_val,
    specificity = spec_val,
    diabetes_OR = exp(diabetes_beta$mean),
    diabetes_OR_lower = exp(diabetes_beta$q5),
    diabetes_OR_upper = exp(diabetes_beta$q95)
  )
}

sens_results_df <- bind_rows(sens_analysis_results)

# Summary statistics across scenarios
sens_summary <- sens_results_df |>
  summarise(
    mean_OR = mean(diabetes_OR),
    sd_OR = sd(diabetes_OR),
    min_OR = min(diabetes_OR),
    max_OR = max(diabetes_OR),
    range_OR = max_OR - min_OR
  )
sens_summary

table_objective_3 <- tibble(
  Analysis = c(
    "Composite Definition Prevalence",
    "Measured BP Only Prevalence",
    "Prior Diagnosis Only Prevalence",
    "",
    "Estimated Sensitivity (Mean [95% CrI])",
    "Estimated Specificity (Mean [95% CrI])",
    "",
    "Diabetes OR - Naive",
    "Diabetes OR - Bias-Corrected",
    "Absolute Difference",
    "",
    "Sensitivity Analysis: OR Range",
    "Sensitivity Analysis: SD"
  ),
  Value = c(
    sprintf("%.1f%%", misclass_components$composite_prev),
    sprintf("%.1f%%", misclass_components$measured_prev),
    sprintf("%.1f%%", misclass_components$diagnosed_prev),
    "",
    sprintf("%.3f [%.3f, %.3f]",
            sens_spec_summary$mean[1],
            sens_spec_summary$q5[1],
            sens_spec_summary$q95[1]),
    sprintf("%.3f [%.3f, %.3f]",
            sens_spec_summary$mean[2],
            sens_spec_summary$q5[2],
            sens_spec_summary$q95[2]),
    "",
    sprintf("%.2f", key_params$Naive_OR[1]),
    sprintf("%.2f", key_params$BiasCorr_OR[1]),
    sprintf("%.2f", key_params$OR_Difference[1]),
    "",
    sprintf("%.2f - %.2f", sens_summary$min_OR, sens_summary$max_OR),
    sprintf("%.2f", sens_summary$sd_OR)
  )
) |>
  gt() |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = Analysis == "")
  )

table_objective_3


# -----------------------------------------------------------
# Objective 4
# -----------------------------------------------------------
# To compare the performance metrics of Bayesian models with logistic models using
# predictive fit

#Cluster-Aware K-Fold Cross-Validation
# PSU-level folds
set.seed(123)
psu_ids <- unique(data_bayes$psu)
n_folds <- 10

# Stratifying on outcome prevalence
psu_prevalence <- data_bayes |>
  group_by(psu) |>
  summarise(prev = mean(hypertension)) |>
  mutate(prev_quartile = ntile(prev, 4))

fold_assignment <- psu_prevalence |>
  group_by(prev_quartile) |>
  mutate(fold = sample(rep(1:n_folds, length.out = n()))) |>
  ungroup() |>
  select(psu, fold)

data_cv <- data_bayes |>
  left_join(fold_assignment, by = "psu")

calc_cv_metrics <- \(predictions, observed, weights = NULL) {
  
  if (is.null(weights)) weights <- rep(1, length(observed))
  
  # Brier score
  brier <- weighted.mean((predictions - observed)^2, weights)
  
  # Log loss
  log_loss <- -weighted.mean(
    observed * log(predictions + 1e-10) + 
      (1 - observed) * log(1 - predictions + 1e-10),
    weights
  )
  
  # Calibration
  calib_data <- tibble(
    pred = predictions,
    obs = observed,
    weight = weights,
    pred_decile = ntile(pred, 10)
  ) |>
    group_by(pred_decile) |>
    summarise(
      mean_pred = weighted.mean(pred, weight),
      mean_obs = weighted.mean(obs, weight),
      n = n()
    ) |>
    filter(n >= 5)
  
  if (nrow(calib_data) > 2) {
    calib_model <- lm(mean_obs ~ mean_pred, data = calib_data)
    calib_intercept <- coef(calib_model)[1]
    calib_slope <- coef(calib_model)[2]
  } else {
    calib_intercept <- NA
    calib_slope <- NA
  }
  
  # AUC (using trapezoidal rule)
  roc_data <- tibble(
    pred = predictions,
    obs = observed
  ) |>
    arrange(desc(pred)) |>
    mutate(
      tpr = cumsum(obs) / sum(obs),
      fpr = cumsum(1 - obs) / sum(1 - obs)
    )
  
  auc <- sum(diff(roc_data$fpr) * (head(roc_data$tpr, -1) + tail(roc_data$tpr, -1)) / 2)
  
  return(list(
    brier = brier,
    log_loss = log_loss,
    calib_intercept = calib_intercept,
    calib_slope = calib_slope,
    auc = auc,
    calibration_data = calib_data
  ))
}

# CV results
models_to_compare <- c("Freq_Logistic", "Freq_Poisson", "Bayes_Single", 
                       "Bayes_2Level", "Bayes_3Level")

cv_results_all <- list()

for (model_name in models_to_compare) {
  cat("Cross-validating:", model_name, "\n")
  
  cv_fold_results <- list()
  
  for (fold_i in 1:n_folds) {
    cat("  Fold", fold_i, "/", n_folds, "\n")
    
    # Split data
    train_data <- data_cv |> filter(fold != fold_i)
    test_data <- data_cv |> filter(fold == fold_i)
    
    if (model_name == "Freq_Logistic") {
      train_svy <- svydesign(
        ids = ~psu, strata = ~stratum, weights = ~wstep2,
        data = train_data, nest = TRUE
      )
      
      model_fold <- svyglm(
        hypertension ~ Age_s + Sex + bmi_cat + has_diabetes + 
          current_smoker + high_alcohol_intake + high_fat_intake + 
          high_salt_intake + Education + wealth + Residence +
          insufficient_activity,
        design = train_svy ,
        family = quasibinomial()
      )
      
      preds <- predict(model_fold, newdata = test_data, type = "response")
      
    } else if (model_name == "Freq_Poisson") {
      train_svy <- svydesign(
        ids = ~psu, strata = ~stratum, weights = ~wstep2,
        data = train_data, nest = TRUE
      )
      
      model_fold <- svyglm(
        hypertension ~ Age_s + Sex + bmi_cat + has_diabetes + 
          current_smoker + high_alcohol_intake + high_fat_intake + 
          high_salt_intake + Education + wealth + Residence +
          insufficient_activity,
        design = train_svy,
        family = quasipoisson()
      )
      
      preds <- predict(model_fold, newdata = test_data, type = "response")
      preds <- pmin(preds, 1) 
      
    } else if (model_name == "Bayes_Single") {
      model_fold <- brm(
        hypertension | weights(weight_norm) ~ 
          Age_s + Sex + bmi_cat + has_diabetes + 
          current_smoker + high_alcohol_intake + high_fat_intake + 
          high_salt_intake + Education + wealth + Residence + 
          insufficient_activity,
        data = train_data,
        family = bernoulli(),
        prior = c(prior(normal(0, 2.5), class = "Intercept"),
                  prior(normal(0, 1), class = "b")),
        chains = 2, iter = 1500, warmup = 750,
        seed = 123 + fold_i,
        refresh = 0,
        backend = "cmdstanr"
      )
      
      preds <- colMeans(posterior_epred(model_fold, newdata = test_data, 
                                        ndraws = 500))
      
    } else if (model_name == "Bayes_2Level") {
      model_fold <- brm(
        hypertension | weights(weight_norm) ~ 
          Age_s + Sex + bmi_cat + has_diabetes + 
          current_smoker + high_alcohol_intake + high_fat_intake + 
          high_salt_intake + Education + wealth + Residence +
          insufficient_activity +
          (1 | household_id),
        data = train_data,
        family = bernoulli(),
        prior = prior_2level,
        chains = 2, iter = 1500, warmup = 750,
        seed = 123 + fold_i,
        refresh = 0,
        backend = "cmdstanr"
      )
      
      preds <- colMeans(posterior_epred(model_fold, newdata = test_data,
                                        ndraws = 500, allow_new_levels = TRUE))
      
    } else if (model_name == "Bayes_3Level") {
      model_fold <- brm(
        hypertension | weights(weight_norm) ~ 
          Age_s + Sex + bmi_cat + has_diabetes + 
          current_smoker + high_alcohol_intake + high_fat_intake + 
          high_salt_intake + Education + wealth + Residence +
          insufficient_activity +
          (1 | psu/household_id),
        data = train_data,
        family = bernoulli(),
        prior = prior_3level,
        chains = 2, iter = 1500, warmup = 750,
        control = list(adapt_delta = 0.95),
        seed = 123 + fold_i,
        refresh = 0,
        backend = "cmdstanr"
      )
      
      preds <- colMeans(posterior_epred(model_fold, newdata = test_data,
                                        ndraws = 500, allow_new_levels = TRUE))
    }
    
    # Metrics
    metrics <- calc_cv_metrics(preds, test_data$hypertension, 
                               test_data$weight_norm)
    
    cv_fold_results[[fold_i]] <- tibble(
      model = model_name,
      fold = fold_i,
      n_train = nrow(train_data),
      n_test = nrow(test_data),
      brier = metrics$brier,
      log_loss = metrics$log_loss,
      calib_intercept = metrics$calib_intercept,
      calib_slope = metrics$calib_slope,
      auc = metrics$auc
    )
  }
  cv_results_all[[model_name]] <- bind_rows(cv_fold_results)
}

table_objective_4 <- cv_results_all |>
  bind_rows() |> 
  group_by(model) |>
  summarise(
    Brier = sprintf("%.4f (%.4f)", mean(brier), sd(brier)/sqrt(n())),
    LogLoss = sprintf("%.3f (%.3f)", mean(log_loss), sd(log_loss)/sqrt(n())),
    AUC = sprintf("%.3f (%.3f)", mean(auc), sd(auc)/sqrt(n())),
    Calib_Slope = sprintf("%.2f (%.2f)", mean(calib_slope, na.rm=TRUE), sd(calib_slope, na.rm=TRUE)/sqrt(n()))
  ) |>
  arrange(desc(model)) |> 
  gt() |>
  cols_label(
    model = "Model Type",
    Brier = "Brier Score",
    LogLoss = "Log Loss ",
    AUC = "AUC",
    Calib_Slope = "Calibration Slope"
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )
table_objective_4
