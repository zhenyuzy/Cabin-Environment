library(lavaan)
library(dplyr)

# Untrimmed structural model after measurement-model screening.
# Data preparation is copied from code_IJHCS.R.

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

model <- '
  # measurement model retained after measurement-model screening
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

  # untrimmed structural model
  ecg              ~ Temp_Cool + Temp_Warm + CO2_High
  eda              ~ Temp_Cool + Temp_Warm + CO2_High
  resp             ~ Temp_Cool + Temp_Warm + CO2_High
  EMG_Activation_N ~ Temp_Cool + Temp_Warm + CO2_High
  temp_mean        ~ Temp_Cool + Temp_Warm + CO2_High

  Cold_Discomfort  ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + EMG_Activation_N + temp_mean
  Heat_Discomfort  ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + EMG_Activation_N + temp_mean
  PAQ              ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + EMG_Activation_N + temp_mean
  KSS              ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + EMG_Activation_N + temp_mean + PAQ + Cold_Discomfort + Heat_Discomfort

  TOR_time         ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + EMG_Activation_N + temp_mean + PAQ + Heat_Discomfort + Cold_Discomfort + KSS + takeover_type_1 + takeover_type_2
  driving          ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + EMG_Activation_N + temp_mean + PAQ + Heat_Discomfort + Cold_Discomfort + KSS + TOR_time + takeover_type_1 + takeover_type_2
  visual           ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + EMG_Activation_N + temp_mean + PAQ + Heat_Discomfort + Cold_Discomfort + KSS + TOR_time + takeover_type_1 + takeover_type_2
'

fit <- sem(
  model,
  data = data,
  missing = "ml",
  cluster = "ID",
  ncpus = 1
)

get_measure <- function(values, name) {
  if (name %in% names(values)) unname(values[[name]]) else NA_real_
}

fit_values <- fitMeasures(fit)
fit_options <- lavInspect(fit, "options")

fit_results <- data.frame(
  model = "Untrimmed_after_measurement_screening",
  converged = lavInspect(fit, "converged"),
  post_check = lavInspect(fit, "post.check"),
  estimator = fit_options$estimator,
  standard_errors = fit_options$se,
  test = paste(fit_options$test, collapse = "; "),
  N_observations = as.numeric(lavInspect(fit, "nobs")),
  N_clusters = n_distinct(data$ID),
  N_free_parameters = as.numeric(lavInspect(fit, "npar")),
  chisq = get_measure(fit_values, "chisq"),
  df = get_measure(fit_values, "df"),
  chisq_df_ratio = get_measure(fit_values, "chisq") / get_measure(fit_values, "df"),
  chisq_p = get_measure(fit_values, "pvalue"),
  GFI = get_measure(fit_values, "gfi"),
  AGFI = get_measure(fit_values, "agfi"),
  CFI = get_measure(fit_values, "cfi"),
  TLI = get_measure(fit_values, "tli"),
  SRMR = get_measure(fit_values, "srmr"),
  RMSEA = get_measure(fit_values, "rmsea"),
  RMSEA_90CI_lower = get_measure(fit_values, "rmsea.ci.lower"),
  RMSEA_90CI_upper = get_measure(fit_values, "rmsea.ci.upper"),
  NFI = get_measure(fit_values, "nfi"),
  RMR = get_measure(fit_values, "rmr"),
  AIC = get_measure(fit_values, "aic"),
  BIC = get_measure(fit_values, "bic"),
  chisq_scaled = get_measure(fit_values, "chisq.scaled"),
  df_scaled = get_measure(fit_values, "df.scaled"),
  chisq_df_ratio_scaled = get_measure(fit_values, "chisq.scaled") / get_measure(fit_values, "df.scaled"),
  chisq_p_scaled = get_measure(fit_values, "pvalue.scaled"),
  CFI_robust = get_measure(fit_values, "cfi.robust"),
  TLI_robust = get_measure(fit_values, "tli.robust"),
  RMSEA_robust = get_measure(fit_values, "rmsea.robust"),
  RMSEA_90CI_lower_robust = get_measure(fit_values, "rmsea.ci.lower.robust"),
  RMSEA_90CI_upper_robust = get_measure(fit_values, "rmsea.ci.upper.robust"),
  stringsAsFactors = FALSE
)

reviewer_fit_table <- fit_results %>%
  transmute(
    model,
    chisq_df = sprintf("%.2f (%d)", chisq, as.integer(df)),
    chisq_df_ratio,
    chisq_p,
    N = N_observations,
    N_clusters,
    GFI,
    AGFI,
    CFI,
    TLI,
    SRMR,
    RMSEA,
    RMSEA_90CI_lower,
    RMSEA_90CI_upper,
    RMSEA_90CI = sprintf("%.3f [%.3f, %.3f]", RMSEA, RMSEA_90CI_lower, RMSEA_90CI_upper)
  )

parameter_results <- parameterEstimates(
  fit,
  ci = TRUE,
  level = 0.95,
  standardized = TRUE
)

structural_paths <- parameter_results %>%
  filter(op == "~") %>%
  mutate(
    path = paste(lhs, "~", rhs),
    significant_05 = pvalue < 0.05,
    .before = 1
  )

measurement_loadings <- parameter_results %>%
  filter(op == "=~")

r2_values <- lavInspect(fit, "r2")
r_squared <- data.frame(
  outcome = names(r2_values),
  r_squared = as.numeric(r2_values),
  row.names = NULL
)

valid_latents <- measurement_loadings %>%
  count(lhs, name = "n_indicators") %>%
  filter(n_indicators >= 2) %>%
  pull(lhs)

cr_ave_table <- measurement_loadings %>%
  filter(lhs %in% valid_latents) %>%
  group_by(lhs) %>%
  summarise(
    CR = (sum(std.all))^2 / ((sum(std.all))^2 + sum(1 - std.all^2)),
    AVE = mean(std.all^2),
    .groups = "drop"
  )

# Comparison with the structurally trimmed original model fitted in the
# sensitivity-analysis script, when that output is available.
trimmed_fit_file <- "IJHCS_sensitivity_fit_comparison.csv"
if (file.exists(trimmed_fit_file)) {
  trimmed_fit <- read.csv(trimmed_fit_file) %>%
    filter(model == "Original") %>%
    transmute(
      model = "Trimmed_structural_model",
      converged,
      post_check = NA,
      estimator,
      standard_errors,
      test,
      N_observations,
      N_clusters,
      N_free_parameters,
      chisq,
      df,
      chisq_df_ratio = chisq / df,
      chisq_p,
      GFI,
      AGFI,
      CFI,
      TLI,
      SRMR,
      RMSEA,
      RMSEA_90CI_lower,
      RMSEA_90CI_upper,
      NFI,
      RMR,
      AIC,
      BIC,
      chisq_scaled,
      df_scaled,
      chisq_df_ratio_scaled = chisq_scaled / df_scaled,
      chisq_p_scaled,
      CFI_robust,
      TLI_robust,
      RMSEA_robust,
      RMSEA_90CI_lower_robust,
      RMSEA_90CI_upper_robust
    )
  untrimmed_trimmed_comparison <- bind_rows(fit_results, trimmed_fit)
} else {
  untrimmed_trimmed_comparison <- fit_results
}

trimmed_paths_file <- "IJHCS_sensitivity_all_structural_paths.csv"
if (file.exists(trimmed_paths_file)) {
  trimmed_paths <- read.csv(trimmed_paths_file) %>%
    filter(model == "Original") %>%
    transmute(
      path,
      trimmed_std = std_coefficient,
      trimmed_p = p_value
    )

  untrimmed_paths_for_comparison <- structural_paths %>%
    transmute(
      path,
      untrimmed_std = std.all,
      untrimmed_p = pvalue
    )

  path_comparison <- full_join(
    trimmed_paths,
    untrimmed_paths_for_comparison,
    by = "path"
  ) %>%
    mutate(
      in_trimmed_model = !is.na(trimmed_std),
      in_untrimmed_model = !is.na(untrimmed_std),
      trimmed_significant = trimmed_p < 0.05,
      untrimmed_significant = untrimmed_p < 0.05,
      standardized_change = untrimmed_std - trimmed_std,
      significance_changed = ifelse(
        in_trimmed_model & in_untrimmed_model,
        trimmed_significant != untrimmed_significant,
        NA
      )
    ) %>%
    arrange(desc(significance_changed), path)
} else {
  path_comparison <- data.frame()
}

write.csv(fit_results, "IJHCS_untrimmed_fit_details.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(reviewer_fit_table, "IJHCS_untrimmed_reviewer_fit_table.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(untrimmed_trimmed_comparison, "IJHCS_untrimmed_vs_trimmed_fit_comparison.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(structural_paths, "IJHCS_untrimmed_structural_paths.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(measurement_loadings, "IJHCS_untrimmed_measurement_loadings.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(parameter_results, "IJHCS_untrimmed_all_parameters.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(r_squared, "IJHCS_untrimmed_r_squared.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cr_ave_table, "IJHCS_untrimmed_cr_ave.csv", row.names = FALSE, fileEncoding = "UTF-8")
if (nrow(path_comparison) > 0) {
  write.csv(
    path_comparison,
    "IJHCS_untrimmed_vs_trimmed_path_comparison.csv",
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

# Scaled chi-square difference test for the explicitly nested structural models.
# Both fitted objects use the yuan.bentler.mplus test, so the same test name must
# be supplied to lavTestLRT().
trimmed_fit_object_file <- "IJHCS_all_fitted_models.rds"
if (file.exists(trimmed_fit_object_file)) {
  trimmed_fit_object <- readRDS(trimmed_fit_object_file)$Original
  same_n <- identical(lavInspect(fit, "nobs"), lavInspect(trimmed_fit_object, "nobs"))
  same_observed_variables <- identical(
    sort(lavNames(fit, "ov")),
    sort(lavNames(trimmed_fit_object, "ov"))
  )
  same_estimator <- identical(
    lavInspect(fit, "options")$estimator,
    lavInspect(trimmed_fit_object, "options")$estimator
  )
  same_test <- identical(
    lavInspect(fit, "options")$test,
    lavInspect(trimmed_fit_object, "options")$test
  )

  scaled_lrt <- lavTestLRT(
    fit,
    trimmed_fit_object,
    test = "yuan.bentler.mplus"
  )

  chisq_difference_result <- data.frame(
    comparison = "Untrimmed vs trimmed structural model",
    method = "Satorra-Bentler 2001 scaled difference test",
    fitted_test = paste(lavInspect(fit, "options")$test, collapse = "; "),
    same_N = same_n,
    same_observed_variables = same_observed_variables,
    same_estimator = same_estimator,
    same_test_setting = same_test,
    untrimmed_df = scaled_lrt[1, "Df"],
    trimmed_df = scaled_lrt[2, "Df"],
    untrimmed_standard_chisq = scaled_lrt[1, "Chisq"],
    trimmed_standard_chisq = scaled_lrt[2, "Chisq"],
    scaled_chisq_difference = scaled_lrt[2, "Chisq diff"],
    df_difference = scaled_lrt[2, "Df diff"],
    p_value = scaled_lrt[2, "Pr(>Chisq)"],
    stringsAsFactors = FALSE
  )

  write.csv(
    chisq_difference_result,
    "IJHCS_untrimmed_vs_trimmed_scaled_chisq_difference.csv",
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  writeLines(
    capture.output(print(scaled_lrt)),
    "IJHCS_untrimmed_vs_trimmed_scaled_chisq_difference.txt"
  )
}

writeLines(
  capture.output(summary(
    fit,
    fit.measures = TRUE,
    standardized = TRUE,
    rsquare = TRUE,
    ci = TRUE
  )),
  "IJHCS_untrimmed_complete_summary.txt"
)
writeLines(model, "IJHCS_untrimmed_model_syntax.txt")
saveRDS(fit, "IJHCS_untrimmed_fitted_model.rds")

cat("\nUntrimmed model convergence\n")
print(c(converged = lavInspect(fit, "converged"), post_check = lavInspect(fit, "post.check")))
cat("\nReviewer fit table\n")
print(reviewer_fit_table, digits = 4)
cat("\nUntrimmed vs trimmed fit comparison\n")
print(untrimmed_trimmed_comparison, digits = 4)
if (exists("scaled_lrt")) {
  cat("\nScaled chi-square difference test\n")
  print(scaled_lrt)
}

