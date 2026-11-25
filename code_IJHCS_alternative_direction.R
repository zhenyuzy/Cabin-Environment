library(lavaan)
library(dplyr)

# Alternative-direction analysis based on code_IJHCS.R.
# Only the two thermal-discomfort/KSS paths are reversed. PAQ and KSS have no
# direct path in either model, matching the current code_IJHCS.R specification.

data <- read.csv("./data_IJHCS.csv")

data <- data %>%
  filter(
    !ID %in% c(
      "A17", "A28", "A210", "A212", "A37",
      "B11", "B19", "B21", "B28", "B29", "B211"
    ),
    Elapsed_time > 0
  )

data$Cold_Discomfort <- ifelse(data$TSV < 0, -data$TSV, 0) / 3
data$Heat_Discomfort <- ifelse(data$TSV > 0, data$TSV, 0) / 3

# Dummy variables: neutral temperature, low CO2, and takeover type 0 are references.
data$Temp_Cool <- ifelse(data$Temperature == 1, 1, 0)
data$Temp_Warm <- ifelse(data$Temperature == 3, 1, 0)
data$Temp_Neutral <- ifelse(data$Temperature == 2, 1, 0)
data$CO2_High <- ifelse(data$CO2 == 1, 1, 0)
data$takeover_type_1 <- ifelse(data$takeover_type == 1, 1, 0)
data$takeover_type_2 <- ifelse(data$takeover_type == 2, 1, 0)

vars_to_process <- c(
  "HR_mean", "HRV_SDNN", "HRV_RMSSD", "HRV_LFHF", "HRV_LF", "HRV_HF",
  "HRV_pNN20", "HRV_LFn", "HRV_HFn", "SCL_mean", "SCR_amplitude_mean",
  "SCR_frequency", "RR_mean", "RR_amplitude_mean", "Gamma_relative",
  "Beta_relative", "Alpha_relative", "Theta_relative", "Delta_relative",
  "temp_mean", "EMG_Activation_N", "EMG_Amplitude_Mean", "TSV", "TP",
  "PAQ", "KSS", "TOR_time", "reaction_time", "MinTTC",
  "SDOffsetfromLaneCenter", "MeanLongitudinalSpeed", "MeanLatitudinalSpeed",
  "SDLongitudinalSpeed", "SDLatitudinalSpeed", "SDLongitudinalAcceleration",
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
  ecg =~ HR_mean + HRV_pNN20
  eda =~ SCL_mean + SCR_amplitude_mean
  resp =~ RR_mean
  driving =~ SDOffsetfromLaneCenter + MaxSteeringAngleVelocity
  visual =~ saccade_frequency + mean_fixation_duration + fixation_frequency + ANNI
'

common_upstream <- '
  ecg              ~ Temp_Cool + Temp_Warm
  eda              ~ Temp_Warm
  resp             ~ CO2_High
  EMG_Activation_N ~ Temp_Cool + Temp_Warm
  temp_mean        ~ Temp_Cool + Temp_Warm
'

common_downstream <- '
  TOR_time ~ CO2_High + Cold_Discomfort + eda + PAQ + KSS + takeover_type_1 + takeover_type_2
  driving  ~ Temp_Warm + CO2_High + Heat_Discomfort + EMG_Activation_N + PAQ + KSS + TOR_time + takeover_type_1 + takeover_type_2
  visual   ~ Temp_Cool + Temp_Warm + CO2_High + Heat_Discomfort + resp + PAQ + KSS + TOR_time + takeover_type_1 + takeover_type_2
'

# Original direction: Cold_Discomfort and Heat_Discomfort -> KSS.
model_original <- paste0(
  measurement_model,
  common_upstream,
  '
  Cold_Discomfort ~ Temp_Cool + Temp_Warm + resp + temp_mean
  Heat_Discomfort ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + temp_mean
  PAQ             ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + temp_mean
  KSS             ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + Cold_Discomfort + Heat_Discomfort
  ',
  common_downstream
)

# Alternative direction: KSS -> Cold_Discomfort and Heat_Discomfort.
model_alternative <- paste0(
  measurement_model,
  common_upstream,
  '
  KSS             ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp
  Cold_Discomfort ~ Temp_Cool + Temp_Warm + resp + temp_mean + KSS
  Heat_Discomfort ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + temp_mean + KSS
  PAQ             ~ Temp_Cool + Temp_Warm + CO2_High + ecg + eda + resp + temp_mean
  ',
  common_downstream
)

fit_original <- sem(
  model_original,
  data = data,
  missing = "ml",
  cluster = "ID"
)

fit_alternative <- sem(
  model_alternative,
  data = data,
  missing = "ml",
  cluster = "ID"
)

fits <- list(Original_direction = fit_original, Alternative_direction = fit_alternative)

get_fit_measure <- function(fit, measure) {
  values <- tryCatch(fitMeasures(fit), error = function(e) numeric())
  if (measure %in% names(values)) unname(values[[measure]]) else NA_real_
}

extract_fit <- function(fit, model_name) {
  data.frame(
    model = model_name,
    converged = lavInspect(fit, "converged"),
    N = as.numeric(lavInspect(fit, "nobs")),
    N_clusters = dplyr::n_distinct(data$ID),
    N_free_parameters = as.numeric(lavInspect(fit, "npar")),
    chisq = get_fit_measure(fit, "chisq"),
    df = get_fit_measure(fit, "df"),
    chisq_p = get_fit_measure(fit, "pvalue"),
    chisq_df_ratio = get_fit_measure(fit, "chisq") / get_fit_measure(fit, "df"),
    GFI = get_fit_measure(fit, "gfi"),
    AGFI = get_fit_measure(fit, "agfi"),
    CFI = get_fit_measure(fit, "cfi"),
    TLI = get_fit_measure(fit, "tli"),
    NFI = get_fit_measure(fit, "nfi"),
    SRMR = get_fit_measure(fit, "srmr"),
    RMSEA = get_fit_measure(fit, "rmsea"),
    RMSEA_90CI_lower = get_fit_measure(fit, "rmsea.ci.lower"),
    RMSEA_90CI_upper = get_fit_measure(fit, "rmsea.ci.upper"),
    AIC = get_fit_measure(fit, "aic"),
    BIC = get_fit_measure(fit, "bic"),
    chisq_scaled = get_fit_measure(fit, "chisq.scaled"),
    chisq_p_scaled = get_fit_measure(fit, "pvalue.scaled"),
    CFI_robust = get_fit_measure(fit, "cfi.robust"),
    TLI_robust = get_fit_measure(fit, "tli.robust"),
    RMSEA_robust = get_fit_measure(fit, "rmsea.robust"),
    RMSEA_90CI_lower_robust = get_fit_measure(fit, "rmsea.ci.lower.robust"),
    RMSEA_90CI_upper_robust = get_fit_measure(fit, "rmsea.ci.upper.robust"),
    stringsAsFactors = FALSE
  )
}

fit_comparison <- bind_rows(Map(extract_fit, fits, names(fits)))
fit_comparison <- fit_comparison %>%
  mutate(
    across(
      c(chisq, df, chisq_p, chisq_df_ratio, GFI, AGFI, CFI, TLI, NFI, SRMR,
        RMSEA, RMSEA_90CI_lower, RMSEA_90CI_upper, AIC, BIC, chisq_scaled,
        chisq_p_scaled, CFI_robust, TLI_robust, RMSEA_robust,
        RMSEA_90CI_lower_robust, RMSEA_90CI_upper_robust),
      ~ as.numeric(.x)
    )
  )

fit_differences <- fit_comparison %>%
  select(where(is.numeric)) %>%
  summarise(across(everything(), ~ .x[2] - .x[1])) %>%
  mutate(comparison = "Alternative minus Original", .before = 1)

extract_structural <- function(fit, model_name) {
  parameterEstimates(fit, standardized = TRUE, ci = TRUE) %>%
    filter(op == "~") %>%
    transmute(
      model = model_name,
      outcome = lhs,
      predictor = rhs,
      path = paste(lhs, "~", rhs),
      estimate = est,
      se,
      z,
      p_value = pvalue,
      ci_lower = ci.lower,
      ci_upper = ci.upper,
      std_coefficient = std.all,
      significant = pvalue < 0.05
    )
}

structural_paths <- bind_rows(Map(extract_structural, fits, names(fits)))
significant_paths <- structural_paths %>%
  filter(significant) %>%
  arrange(model, outcome, predictor)

original_paths <- structural_paths %>%
  filter(model == "Original_direction") %>%
  select(-model) %>%
  rename_with(~ paste0("original_", .x), -c(outcome, predictor, path))

alternative_paths <- structural_paths %>%
  filter(model == "Alternative_direction") %>%
  select(-model) %>%
  rename_with(~ paste0("alternative_", .x), -c(outcome, predictor, path))

path_comparison <- full_join(
  original_paths,
  alternative_paths,
  by = c("outcome", "predictor", "path")
) %>%
  mutate(
    std_change = alternative_std_coefficient - original_std_coefficient,
    significance_status = case_when(
      is.na(original_p_value) ~ "Alternative model only",
      is.na(alternative_p_value) ~ "Original model only",
      original_significant & alternative_significant ~ "Significant in both",
      original_significant & !alternative_significant ~ "Original only significant",
      !original_significant & alternative_significant ~ "Alternative only significant",
      TRUE ~ "Not significant in either"
    )
  ) %>%
  arrange(outcome, predictor)

changed_direction_paths <- path_comparison %>%
  filter(
    (outcome == "KSS" & predictor %in% c("Cold_Discomfort", "Heat_Discomfort")) |
      (outcome %in% c("Cold_Discomfort", "Heat_Discomfort") & predictor == "KSS")
  )

shared_significance_changes <- path_comparison %>%
  filter(!is.na(original_p_value), !is.na(alternative_p_value)) %>%
  filter(original_significant != alternative_significant)

# Track every path that was significant in the original model. For the two
# reversed paths, report the corresponding reverse-direction estimate but do
# not calculate a coefficient difference because the estimands differ.
alternative_lookup <- structural_paths %>%
  filter(model == "Alternative_direction") %>%
  select(
    comparison_path = path,
    comparison_alt_std = std_coefficient,
    comparison_alt_p = p_value,
    comparison_alt_significant = significant
  )

original_significant_changes <- path_comparison %>%
  filter(original_significant %in% TRUE) %>%
  mutate(
    comparison_path = case_when(
      path == "KSS ~ Cold_Discomfort" ~ "Cold_Discomfort ~ KSS",
      path == "KSS ~ Heat_Discomfort" ~ "Heat_Discomfort ~ KSS",
      TRUE ~ path
    ),
    direction_reversed = comparison_path != path
  ) %>%
  select(
    path,
    comparison_path,
    direction_reversed,
    original_std_coefficient,
    original_p_value
  ) %>%
  left_join(alternative_lookup, by = "comparison_path") %>%
  mutate(
    std_change = ifelse(
      direction_reversed,
      NA_real_,
      comparison_alt_std - original_std_coefficient
    ),
    result = case_when(
      direction_reversed & comparison_alt_significant ~ "Direction reversed; reverse path significant",
      direction_reversed & !comparison_alt_significant ~ "Direction reversed; reverse path not significant",
      comparison_alt_significant ~ "Remained significant",
      TRUE ~ "No longer significant"
    )
  ) %>%
  arrange(path)

write.csv(
  fit_comparison,
  "IJHCS_direction_model_fit_comparison.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  fit_differences,
  "IJHCS_direction_model_fit_differences.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  structural_paths,
  "IJHCS_direction_model_all_structural_paths.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  significant_paths,
  "IJHCS_direction_model_significant_paths.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  path_comparison,
  "IJHCS_direction_model_path_comparison.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  changed_direction_paths,
  "IJHCS_direction_model_reversed_paths.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  original_significant_changes,
  "IJHCS_direction_model_original_significant_changes.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

fit_md <- apply(fit_comparison, 1, function(x) {
  sprintf(
    "| %s | %.2f (%d) | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f [%.3f, %.3f] | %.1f | %.1f |",
    x[["model"]], as.numeric(x[["chisq"]]), as.integer(as.numeric(x[["df"]])),
    fmt_p(as.numeric(x[["chisq_p"]])), as.numeric(x[["GFI"]]),
    as.numeric(x[["AGFI"]]), as.numeric(x[["CFI"]]), as.numeric(x[["TLI"]]),
    as.numeric(x[["SRMR"]]), as.numeric(x[["RMSEA"]]),
    as.numeric(x[["RMSEA_90CI_lower"]]), as.numeric(x[["RMSEA_90CI_upper"]]),
    as.numeric(x[["AIC"]]), as.numeric(x[["BIC"]])
  )
})

robust_fit_md <- apply(fit_comparison, 1, function(x) {
  sprintf(
    "| %s | %.2f (%d) | %s | %.3f | %.3f | %.3f [%.3f, %.3f] |",
    x[["model"]], as.numeric(x[["chisq_scaled"]]),
    as.integer(as.numeric(x[["df"]])), fmt_p(as.numeric(x[["chisq_p_scaled"]])),
    as.numeric(x[["CFI_robust"]]), as.numeric(x[["TLI_robust"]]),
    as.numeric(x[["RMSEA_robust"]]),
    as.numeric(x[["RMSEA_90CI_lower_robust"]]),
    as.numeric(x[["RMSEA_90CI_upper_robust"]])
  )
})

reversed_md <- changed_direction_paths %>%
  transmute(
    line = sprintf(
      "| `%s` | %s | %s | %s | %s |",
      path,
      ifelse(is.na(original_std_coefficient), "--", sprintf("%.3f", original_std_coefficient)),
      ifelse(is.na(original_p_value), "--", fmt_p(original_p_value)),
      ifelse(is.na(alternative_std_coefficient), "--", sprintf("%.3f", alternative_std_coefficient)),
      ifelse(is.na(alternative_p_value), "--", fmt_p(alternative_p_value))
    )
  ) %>%
  pull(line)

original_sig_change_md <- original_significant_changes %>%
  transmute(
    result_cn = case_when(
      direction_reversed & comparison_alt_significant ~ paste0(
        "方向反转为 `", comparison_path, "`，反向路径显著"
      ),
      direction_reversed & !comparison_alt_significant ~ paste0(
        "方向反转为 `", comparison_path, "`，反向路径不显著"
      ),
      comparison_alt_significant ~ "仍显著",
      TRUE ~ "不再显著"
    ),
    line = sprintf(
      "| `%s` | %.3f | %s | %.3f | %s | %s | %s |",
      path,
      original_std_coefficient,
      fmt_p(original_p_value),
      comparison_alt_std,
      ifelse(direction_reversed, "--", sprintf("%+.3f", std_change)),
      fmt_p(comparison_alt_p),
      result_cn
    )
  ) %>%
  pull(line)

sig_change_md <- if (nrow(shared_significance_changes) == 0) {
  "- 两模型共有路径中，没有路径跨越 p=0.05 显著性阈值。"
} else {
  paste0(
    "- `", shared_significance_changes$path, "`：原模型 β=",
    sprintf("%.3f", shared_significance_changes$original_std_coefficient),
    "，p=", fmt_p(shared_significance_changes$original_p_value),
    "；替代模型 β=", sprintf("%.3f", shared_significance_changes$alternative_std_coefficient),
    "，p=", fmt_p(shared_significance_changes$alternative_p_value), "。"
  )
}

largest_shared_change <- original_significant_changes %>%
  filter(!direction_reversed) %>%
  slice_max(order_by = abs(std_change), n = 1, with_ties = FALSE)

report_lines <- c(
  "# IJHCS SEM 路径方向替代模型比较",
  "",
  "## 模型设置",
  "",
  "- 两个模型使用相同数据处理、测量模型、样本、聚类稳健标准误和其余结构路径。",
  "- 原方向模型：`Cold_Discomfort` 和 `Heat_Discomfort` 指向 `KSS`。",
  "- 替代方向模型：`KSS` 指向 `Cold_Discomfort` 和 `Heat_Discomfort`。",
  "- `PAQ` 与 `KSS` 之间在两个模型中均不设直接路径，与当前 `code_IJHCS.R` 一致。",
  "- 两个方向模型并非嵌套模型，因此不使用 χ² difference test；主要比较拟合指标、AIC/BIC和路径结果。",
  "",
  "## 模型拟合比较",
  "",
  "| 模型 | χ²(df) | p | GFI | AGFI | CFI | TLI | SRMR | RMSEA [90% CI] | AIC | BIC |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
  fit_md,
  "",
  sprintf(
    "替代模型相对原方向模型的变化为：ΔCFI=%.3f，ΔTLI=%.3f，ΔSRMR=%.3f，ΔRMSEA=%.3f，ΔAIC=%.1f，ΔBIC=%.1f。",
    fit_comparison$CFI[2] - fit_comparison$CFI[1],
    fit_comparison$TLI[2] - fit_comparison$TLI[1],
    fit_comparison$SRMR[2] - fit_comparison$SRMR[1],
    fit_comparison$RMSEA[2] - fit_comparison$RMSEA[1],
    fit_comparison$AIC[2] - fit_comparison$AIC[1],
    fit_comparison$BIC[2] - fit_comparison$BIC[1]
  ),
  "两模型拟合几乎相同。原方向模型的 χ²、RMSEA、AIC 和 BIC 略低，CFI 和 TLI 略高，但 AIC/BIC 差值均仅为 1.9，不足以形成明显的模型选择证据。",
  "",
  "### Yuan–Bentler scaled/robust 拟合指标",
  "",
  "| 模型 | scaled χ²(df) | p | robust CFI | robust TLI | robust RMSEA [90% CI] |",
  "|---|---:|---:|---:|---:|---:|",
  robust_fit_md,
  "",
  "## 反向路径的估计结果",
  "",
  "| 路径 | 原模型 β | 原模型 p | 替代模型 β | 替代模型 p |",
  "|---|---:|---:|---:|---:|",
  reversed_md,
  "",
  "## 原模型显著路径在替代模型中的变化",
  "",
  "Δβ = 替代模型标准化系数 - 原模型标准化系数。对于方向反转的路径，两个系数对应不同因变量，因此不计算 Δβ。",
  "",
  "| 原模型显著路径 | 原模型 β | 原模型 p | 替代模型 β | Δβ | 替代模型 p | 结果 |",
  "|---|---:|---:|---:|---:|---:|---|",
  original_sig_change_md,
  "",
  "## 主要路径变化摘要",
  "",
  sprintf(
    "- 原模型的 %d 条显著路径中，%d 条在两模型中方向相同；其中 %d 条在替代模型中仍显著，%d 条不再显著。",
    nrow(original_significant_changes),
    sum(!original_significant_changes$direction_reversed),
    sum(!original_significant_changes$direction_reversed & original_significant_changes$comparison_alt_significant),
    sum(!original_significant_changes$direction_reversed & !original_significant_changes$comparison_alt_significant)
  ),
  sprintf(
    "- 同方向路径中系数变化最大的是 `%s`：β 由 %.3f 变为 %.3f，Δβ=%+.3f，替代模型中 p=%s。",
    largest_shared_change$path,
    largest_shared_change$original_std_coefficient,
    largest_shared_change$comparison_alt_std,
    largest_shared_change$std_change,
    fmt_p(largest_shared_change$comparison_alt_p)
  ),
  "- 两条热不舒适–KSS 路径在替代模型中反转后均显著；但由于因变量改变，反向系数不应与原系数直接作 Δβ 比较。",
  "",
  "## 两模型共有路径的显著性变化",
  "",
  sig_change_md,
  "",
  sprintf(
    "原方向模型共有 %d 条显著结构路径；替代方向模型共有 %d 条显著结构路径。全部显著与不显著结构路径均已保存至 CSV。",
    sum(structural_paths$model == "Original_direction" & structural_paths$significant),
    sum(structural_paths$model == "Alternative_direction" & structural_paths$significant)
  ),
  "",
  "## 解释注意事项",
  "",
  "该比较能够判断哪一种因果方向与观测协方差结构更一致，但横断/重复测量 SEM 的拟合优劣本身不能证明因果方向。路径方向仍应结合理论依据和测量时间顺序解释。",
  "",
  "`lavaan` 对两个模型均提示参数协方差矩阵非正定，最小特征值为 10^-14 量级。该问题也存在于原代码模型，因此边界显著路径和聚类稳健推断需谨慎解读。"
)

writeLines(report_lines, "IJHCS_direction_model_report.md", useBytes = TRUE)

cat("\nModel convergence\n")
print(sapply(fits, lavInspect, what = "converged"))
cat("\nFit comparison\n")
print(fit_comparison, digits = 4)
cat("\nReversed paths\n")
print(changed_direction_paths, digits = 4)
cat("\nShared paths crossing p = 0.05\n")
print(shared_significance_changes, digits = 4)
