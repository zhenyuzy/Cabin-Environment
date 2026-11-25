library(lavaan)
library(dplyr)

# Sensitivity analysis copied from code_IJHCS.R.
# The original model is kept unchanged. Three sensitivity models separately
# adjust for elapsed time, categorical humidity, and measured relative humidity.

data <- read.csv("./data_IJHCS.csv")

data <- data %>%
  filter(
    !ID %in% c(
      "A17", "A28", "A210", "A212", "A37",
      "B11", "B19", "B21", "B28", "B29", "B211"
    ),
    Elapsed_time > 0
  )

data$Cold_Discomfort <- ifelse(data$TSV < 0, -data$TSV, 0)
data$Heat_Discomfort <- ifelse(data$TSV > 0, data$TSV, 0)
data$Cold_Discomfort <- data$Cold_Discomfort / 3
data$Heat_Discomfort <- data$Heat_Discomfort / 3

# Temperature dummy: reference = 2 (neutral)
data$Temp_Cool <- ifelse(data$Temperature == 1, 1, 0)
data$Temp_Warm <- ifelse(data$Temperature == 3, 1, 0)
data$Temp_Neutral <- ifelse(data$Temperature == 2, 1, 0)

# CO2 dummy: reference = 0 (low)
data$CO2_High <- ifelse(data$CO2 == 1, 1, 0)

# Takeover type dummy: reference = 0
data$takeover_type_1 <- ifelse(data$takeover_type == 1, 1, 0)
data$takeover_type_2 <- ifelse(data$takeover_type == 2, 1, 0)

# Participant-level mean relative humidity measured during each experiment.
# A/B indicate low/high CO2; 1/2/3 indicate the temperature condition.
env_data_dir <- "/Users/zhenyu/Documents/CabinEXP1DATA/EnvData"
env_groups <- c("A1", "A2", "A3", "B1", "B2", "B3")

env_humidity <- bind_rows(lapply(env_groups, function(group_name) {
  env_file <- file.path(env_data_dir, paste0(group_name, "_Mean_Env_Data.csv"))
  env_data <- read.csv(env_file, check.names = FALSE)
  data.frame(
    ID = as.character(env_data$ID),
    Env_group = group_name,
    Humidity_actual = as.numeric(env_data[["Mean_RH (%)"]])
  )
}))

if (anyDuplicated(env_humidity$ID)) {
  stop("Duplicate participant IDs were found in the environmental files.")
}

group_humidity_means <- env_humidity %>%
  group_by(Env_group) %>%
  summarise(Humidity_group_mean = mean(Humidity_actual, na.rm = TRUE), .groups = "drop")

data$Env_group <- paste0(ifelse(data$CO2 == 0, "A", "B"), data$Temperature)
data <- data %>%
  left_join(env_humidity, by = c("ID", "Env_group")) %>%
  left_join(group_humidity_means, by = "Env_group") %>%
  mutate(
    Humidity_imputed = is.na(Humidity_actual),
    Humidity_RH = coalesce(Humidity_actual, Humidity_group_mean)
  )

if (anyNA(data$Humidity_RH)) {
  stop("Humidity remains missing after within-group mean imputation.")
}

# Scaled covariates improve numerical conditioning. Because the comparison uses
# std.all, this does not change their standardized coefficients.
data$Elapsed_time_cov <- as.numeric(scale(data$Elapsed_time))
humidity_category_map <- c(
  "10" = 2, "11" = 1, "20" = 2,
  "21" = 3, "30" = 2, "31" = 3
)
condition_code <- paste0(data$Temperature, data$CO2)
data$Humidity_category <- unname(humidity_category_map[condition_code])

if (anyNA(data$Humidity_category)) {
  stop("Some Temperature x CO2 groups do not have a humidity category mapping.")
}

data$Humidity_category_cov <- as.numeric(scale(data$Humidity_category))
data$Humidity_actual_cov <- as.numeric(scale(data$Humidity_RH))

vars_to_process <- c(
  "HR_mean", "HRV_SDNN", "HRV_RMSSD", "HRV_LFHF", "HRV_LF", "HRV_HF",
  "HRV_pNN20", "HRV_LFn", "HRV_HFn", "SCL_mean", "SCR_amplitude_mean",
  "SCR_frequency", "RR_mean", "RR_amplitude_mean", "Gamma_relative",
  "Beta_relative", "Alpha_relative", "Theta_relative", "Delta_relative",
  "temp_mean", "EMG_Activation_N", "EMG_Amplitude_Mean", "TSV", "TP", "PAQ",
  "KSS", "TOR_time", "reaction_time", "MinTTC", "SDOffsetfromLaneCenter",
  "MeanLongitudinalSpeed", "MeanLatitudinalSpeed", "SDLongitudinalSpeed",
  "SDLatitudinalSpeed", "SDLongitudinalAcceleration",
  "MinLongitudinalAcceleration", "MeanLongitudinalAcceleration",
  "MaxLongitudinalAcceleration", "SDLatitudinalAcceleration",
  "MeanOffsetfromLaneCenter", "SDSteeringAngleVelocity",
  "RateofSteeringAngleReversal", "MaxBrakeDepth", "saccade_frequency",
  "saccade_peak_velocity", "saccade_mean_velocity", "mean_saccade_duration",
  "mean_fixation_duration", "fixation_frequency", "mean_pupil_diameter",
  "x_std", "y_std", "ANNI", "mean_blink_duration", "blink_frequency"
)

replace_outliers <- function(x) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- qnt[2] - qnt[1]
  if (iqr == 0) return(x)
  lower <- qnt[1] - 3 * iqr
  upper <- qnt[2] + 3 * iqr
  x[x < lower] <- lower
  x[x > upper] <- upper
  x
}

data[vars_to_process] <- lapply(data[vars_to_process], replace_outliers)
data[vars_to_process] <- scale(data[vars_to_process])

measurement_model <- '
  ecg =~ HR_mean
        + HRV_pNN20
  eda =~ SCL_mean
        + SCR_amplitude_mean
  resp =~ RR_mean
  driving =~ SDOffsetfromLaneCenter
            + MaxSteeringAngleVelocity
  visual =~ saccade_frequency
            + mean_fixation_duration
            + fixation_frequency
            + ANNI
'

structural_model <- '
  ecg              ~ Temp_Cool + Temp_Warm
  eda              ~ Temp_Warm
  resp             ~ CO2_High
  EMG_Activation_N ~ Temp_Cool + Temp_Warm
  temp_mean        ~ Temp_Cool + Temp_Warm

  Cold_Discomfort  ~ Temp_Cool + Temp_Warm + resp + temp_mean
  Heat_Discomfort  ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + temp_mean
  PAQ              ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + temp_mean
  KSS              ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + Cold_Discomfort + Heat_Discomfort

  TOR_time         ~ CO2_High + Cold_Discomfort + eda + PAQ + KSS + takeover_type_1 + takeover_type_2
  driving          ~ Temp_Warm + CO2_High + Heat_Discomfort + EMG_Activation_N + PAQ + KSS + TOR_time + takeover_type_1 + takeover_type_2
  visual           ~ Temp_Cool + Temp_Warm + CO2_High + Heat_Discomfort + resp + PAQ + KSS + TOR_time + takeover_type_1 + takeover_type_2
'

# All three sensitivity models use the same adjustment targets: physiological
# outcomes, thermal discomfort, perceived air quality (PAQ), and fatigue (KSS).
elapsed_adjustment <- '
  ecg              ~ Elapsed_time_cov
  eda              ~ Elapsed_time_cov
  resp             ~ Elapsed_time_cov
  EMG_Activation_N ~ Elapsed_time_cov
  temp_mean        ~ Elapsed_time_cov
  Cold_Discomfort  ~ Elapsed_time_cov
  Heat_Discomfort  ~ Elapsed_time_cov
  PAQ              ~ Elapsed_time_cov
  KSS              ~ Elapsed_time_cov
'

humidity_category_adjustment <- '
  ecg              ~ Humidity_category_cov
  eda              ~ Humidity_category_cov
  resp             ~ Humidity_category_cov
  EMG_Activation_N ~ Humidity_category_cov
  temp_mean        ~ Humidity_category_cov
  Cold_Discomfort  ~ Humidity_category_cov
  Heat_Discomfort  ~ Humidity_category_cov
  PAQ              ~ Humidity_category_cov
  KSS              ~ Humidity_category_cov
'

humidity_actual_adjustment <- '
  ecg              ~ Humidity_actual_cov
  eda              ~ Humidity_actual_cov
  resp             ~ Humidity_actual_cov
  EMG_Activation_N ~ Humidity_actual_cov
  temp_mean        ~ Humidity_actual_cov
  Cold_Discomfort  ~ Humidity_actual_cov
  Heat_Discomfort  ~ Humidity_actual_cov
  PAQ              ~ Humidity_actual_cov
  KSS              ~ Humidity_actual_cov
'

model_original <- paste(measurement_model, structural_model, sep = "\n")
model_elapsed <- paste(model_original, elapsed_adjustment, sep = "\n")
model_humidity_category <- paste(model_original, humidity_category_adjustment, sep = "\n")
model_humidity_actual <- paste(model_original, humidity_actual_adjustment, sep = "\n")

fit_sem <- function(model_syntax) {
  sem(
    model_syntax,
    data = data,
    missing = "ml",
    cluster = "ID",
    ncpus = 1
  )
}

fit_original <- fit_sem(model_original)
fit_elapsed <- fit_sem(model_elapsed)
fit_humidity_category <- fit_sem(model_humidity_category)
fit_humidity_actual <- fit_sem(model_humidity_actual)

fit_list <- list(
  Original = fit_original,
  Elapsed_time_adjusted = fit_elapsed,
  Humidity_category_adjusted = fit_humidity_category,
  Humidity_actual_adjusted = fit_humidity_actual
)

extract_structural <- function(fit, model_name) {
  parameterEstimates(fit, standardized = TRUE) %>%
    filter(op == "~") %>%
    transmute(
      outcome = lhs,
      predictor = rhs,
      path = paste(lhs, "~", rhs),
      model = model_name,
      estimate = est,
      se = se,
      z = z,
      p_value = pvalue,
      std_coefficient = std.all
    )
}

all_structural_paths <- bind_rows(
  Map(extract_structural, fit_list, names(fit_list))
)

original_significant <- all_structural_paths %>%
  filter(model == "Original", p_value < 0.05) %>%
  select(outcome, predictor, path, original_std = std_coefficient, original_p = p_value)

original_all_paths <- all_structural_paths %>%
  filter(model == "Original") %>%
  select(outcome, predictor, path, original_std = std_coefficient, original_p = p_value)

elapsed_values <- all_structural_paths %>%
  filter(model == "Elapsed_time_adjusted") %>%
  select(path, elapsed_std = std_coefficient, elapsed_p = p_value)

humidity_category_values <- all_structural_paths %>%
  filter(model == "Humidity_category_adjusted") %>%
  select(
    path,
    humidity_category_std = std_coefficient,
    humidity_category_p = p_value
  )

humidity_actual_values <- all_structural_paths %>%
  filter(model == "Humidity_actual_adjusted") %>%
  select(
    path,
    humidity_actual_std = std_coefficient,
    humidity_actual_p = p_value
  )

significant_path_comparison <- original_significant %>%
  left_join(elapsed_values, by = "path") %>%
  left_join(humidity_category_values, by = "path") %>%
  left_join(humidity_actual_values, by = "path") %>%
  mutate(
    elapsed_change = elapsed_std - original_std,
    elapsed_percent_change = 100 * elapsed_change / abs(original_std),
    humidity_category_change = humidity_category_std - original_std,
    humidity_category_percent_change = 100 * humidity_category_change / abs(original_std),
    humidity_actual_change = humidity_actual_std - original_std,
    humidity_actual_percent_change = 100 * humidity_actual_change / abs(original_std),
    elapsed_significant = elapsed_p < 0.05,
    humidity_category_significant = humidity_category_p < 0.05,
    humidity_actual_significant = humidity_actual_p < 0.05
  ) %>%
  arrange(outcome, predictor)

all_path_comparison <- original_all_paths %>%
  left_join(elapsed_values, by = "path") %>%
  left_join(humidity_category_values, by = "path") %>%
  left_join(humidity_actual_values, by = "path") %>%
  mutate(
    elapsed_change = elapsed_std - original_std,
    humidity_category_change = humidity_category_std - original_std,
    humidity_actual_change = humidity_actual_std - original_std,
    original_significant = original_p < 0.05,
    elapsed_significant = elapsed_p < 0.05,
    humidity_category_significant = humidity_category_p < 0.05,
    humidity_actual_significant = humidity_actual_p < 0.05
  ) %>%
  arrange(outcome, predictor)

get_fit_measure <- function(fit, measure_name) {
  values <- fitMeasures(fit)
  if (measure_name %in% names(values)) {
    return(unname(values[[measure_name]]))
  }
  NA_real_
}

extract_fit_results <- function(fit, model_name) {
  fit_options <- lavInspect(fit, "options")
  data.frame(
    model = model_name,
    converged = lavInspect(fit, "converged"),
    estimator = fit_options$estimator,
    standard_errors = fit_options$se,
    test = paste(fit_options$test, collapse = "; "),
    N_observations = as.numeric(lavInspect(fit, "nobs")),
    N_clusters = n_distinct(data$ID),
    N_free_parameters = as.numeric(lavInspect(fit, "npar")),
    chisq = get_fit_measure(fit, "chisq"),
    df = get_fit_measure(fit, "df"),
    chisq_p = get_fit_measure(fit, "pvalue"),
    RMSEA = get_fit_measure(fit, "rmsea"),
    RMSEA_90CI_lower = get_fit_measure(fit, "rmsea.ci.lower"),
    RMSEA_90CI_upper = get_fit_measure(fit, "rmsea.ci.upper"),
    CFI = get_fit_measure(fit, "cfi"),
    TLI = get_fit_measure(fit, "tli"),
    SRMR = get_fit_measure(fit, "srmr"),
    NFI = get_fit_measure(fit, "nfi"),
    GFI = get_fit_measure(fit, "gfi"),
    AGFI = get_fit_measure(fit, "agfi"),
    RMR = get_fit_measure(fit, "rmr"),
    AIC = get_fit_measure(fit, "aic"),
    BIC = get_fit_measure(fit, "bic"),
    chisq_scaled = get_fit_measure(fit, "chisq.scaled"),
    df_scaled = get_fit_measure(fit, "df.scaled"),
    chisq_p_scaled = get_fit_measure(fit, "pvalue.scaled"),
    RMSEA_robust = get_fit_measure(fit, "rmsea.robust"),
    RMSEA_90CI_lower_robust = get_fit_measure(fit, "rmsea.ci.lower.robust"),
    RMSEA_90CI_upper_robust = get_fit_measure(fit, "rmsea.ci.upper.robust"),
    CFI_robust = get_fit_measure(fit, "cfi.robust"),
    TLI_robust = get_fit_measure(fit, "tli.robust"),
    stringsAsFactors = FALSE
  )
}

fit_comparison <- bind_rows(Map(extract_fit_results, fit_list, names(fit_list)))

reviewer_fit_table <- fit_comparison %>%
  transmute(
    model,
    model_role = ifelse(model == "Original", "Full specified model", "Sensitivity model"),
    chisq_df = sprintf("%.2f (%d)", chisq, as.integer(df)),
    chisq_p,
    N = N_observations,
    N_clusters,
    GFI,
    AGFI,
    RMSEA,
    RMSEA_90CI_lower,
    RMSEA_90CI_upper,
    RMSEA_90CI = sprintf(
      "%.3f [%.3f, %.3f]", RMSEA, RMSEA_90CI_lower, RMSEA_90CI_upper
    ),
    CFI,
    TLI,
    SRMR
  )

all_parameter_estimates <- bind_rows(lapply(names(fit_list), function(model_name) {
  parameterEstimates(
    fit_list[[model_name]],
    ci = TRUE,
    level = 0.95,
    standardized = TRUE
  ) %>%
    mutate(model = model_name, .before = 1)
}))

all_r_squared <- bind_rows(lapply(names(fit_list), function(model_name) {
  r2_values <- lavInspect(fit_list[[model_name]], "r2")
  data.frame(
    model = model_name,
    outcome = names(r2_values),
    r_squared = as.numeric(r2_values),
    row.names = NULL
  )
}))

covariate_effects <- all_structural_paths %>%
  filter(predictor %in% c(
    "Elapsed_time_cov", "Humidity_category_cov", "Humidity_actual_cov"
  ))

write.csv(
  significant_path_comparison,
  "IJHCS_significant_path_sensitivity_comparison.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  all_path_comparison,
  "IJHCS_all_path_sensitivity_comparison.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  fit_comparison,
  "IJHCS_sensitivity_fit_comparison.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  reviewer_fit_table,
  "IJHCS_reviewer_fit_reporting_table.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  all_parameter_estimates,
  "IJHCS_all_models_parameter_estimates.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  all_r_squared,
  "IJHCS_all_models_r_squared.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  covariate_effects,
  "IJHCS_sensitivity_covariate_effects.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  all_structural_paths,
  "IJHCS_sensitivity_all_structural_paths.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

humidity_merge_audit <- data %>%
  distinct(
    ID, Env_group, Humidity_actual, Humidity_group_mean,
    Humidity_imputed, Humidity_RH
  ) %>%
  arrange(Env_group, ID)

write.csv(
  humidity_merge_audit,
  "IJHCS_subject_humidity_merge_audit.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

model_summary_text <- unlist(lapply(names(fit_list), function(model_name) {
  c(
    paste0("========== ", model_name, " =========="),
    capture.output(summary(
      fit_list[[model_name]],
      fit.measures = TRUE,
      standardized = TRUE,
      rsquare = TRUE,
      ci = TRUE
    )),
    ""
  )
}))
writeLines(model_summary_text, "IJHCS_all_models_complete_summaries.txt")

model_syntax_text <- c(
  "========== Original (full specified model) ==========", model_original, "",
  "========== Elapsed time sensitivity ==========", model_elapsed, "",
  "========== Categorical humidity sensitivity ==========", model_humidity_category, "",
  "========== Actual humidity sensitivity ==========", model_humidity_actual
)
writeLines(model_syntax_text, "IJHCS_all_models_syntax.txt")
saveRDS(fit_list, "IJHCS_all_fitted_models.rds")

cat("\nSample and condition summary\n")
cat("Rows:", nrow(data), " Unique participants:", n_distinct(data$ID), "\n")
print(table(data$Elapsed_time, useNA = "ifany"))
print(humidity_merge_audit %>% count(Env_group, Humidity_imputed))
cat("Imputed participant IDs:",
    paste(humidity_merge_audit$ID[humidity_merge_audit$Humidity_imputed], collapse = ", "),
    "\n")
cat("Humidity RH range:", paste(range(data$Humidity_RH), collapse = " to "), "\n")

cat("\nModel convergence\n")
print(vapply(fit_list, lavInspect, logical(1), what = "converged"))

cat("\nFit comparison\n")
print(fit_comparison, digits = 4)

cat("\nReviewer-requested fit reporting table\n")
print(reviewer_fit_table, digits = 4)

cat("\nCovariate effects\n")
print(covariate_effects, digits = 4)

cat("\nOriginally significant structural paths and sensitivity comparison\n")
print(significant_path_comparison, digits = 4)

