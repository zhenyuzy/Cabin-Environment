library(lavaan)
library(dplyr)    
library(mice)     
library(naniar)   
library(semPlot)
library(semTools)  

library(survey)
library(lavaan.survey)

data <- read.csv("./data_IJHCS.csv")


data <- data %>% 
  filter(
         !ID %in% c('A17', 'A28', 'A210', 'A212', 'A37',
                    'B11', 'B19', 'B21', 'B28', 'B29', 'B211'),
         Elapsed_time > 0,
         )

data$Cold_Discomfort <- ifelse(data$TSV < 0, -data$TSV, 0)
data$Heat_Discomfort <- ifelse(data$TSV > 0,  data$TSV, 0)
data$Cold_Discomfort <- data$Cold_Discomfort / 3  
data$Heat_Discomfort <- data$Heat_Discomfort / 3

# Temperature dummy：reference = 2（Neutral）
data$Temp_Cool <- ifelse(data$Temperature == 1, 1, 0)
data$Temp_Warm    <- ifelse(data$Temperature == 3, 1, 0)
data$Temp_Neutral <- ifelse(data$Temperature == 2, 1, 0)
# CO2 dummy：reference = 0（Low）
data$CO2_High <- ifelse(data$CO2 == 1, 1, 0)
# takeover_type dummy：reference = 0
data$takeover_type_1 <- ifelse(data$takeover_type == 1, 1, 0)
data$takeover_type_2 <- ifelse(data$takeover_type == 2, 1, 0)



vars_to_process <- c(
  # "Temperature", "CO2",
  "HR_mean", "HRV_SDNN", "HRV_RMSSD", "HRV_LFHF", "HRV_LF", "HRV_HF", "HRV_pNN20", "HRV_LFn", "HRV_HFn",
  "SCL_mean", "SCR_amplitude_mean", "SCR_frequency",
  "RR_mean", "RR_amplitude_mean",
  "Gamma_relative", "Beta_relative", "Alpha_relative", "Theta_relative", "Delta_relative",
  "temp_mean",
  "EMG_Activation_N", "EMG_Amplitude_Mean",
  "TSV", "TP",
  # "Cold_Discomfort", "Heat_Discomfort",
  "PAQ", "KSS",
  "TOR_time", "reaction_time",
  "MinTTC", "SDOffsetfromLaneCenter", 
  "MeanLongitudinalSpeed", "MeanLatitudinalSpeed", "SDLongitudinalSpeed", "SDLatitudinalSpeed",
  "SDLongitudinalAcceleration", "MinLongitudinalAcceleration", "MeanLongitudinalAcceleration","MaxLongitudinalAcceleration",
  "SDLatitudinalAcceleration", "MeanOffsetfromLaneCenter", 
  "SDSteeringAngleVelocity", "RateofSteeringAngleReversal",
  "MaxBrakeDepth",
  "saccade_frequency", "saccade_peak_velocity", "saccade_mean_velocity", "mean_saccade_duration",
  "mean_fixation_duration", "fixation_frequency", "mean_pupil_diameter", "x_std", "y_std", "ANNI", "mean_blink_duration", "blink_frequency"
)

replace_outliers <- function(x) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- qnt[2] - qnt[1]
  if (iqr == 0) return(x)  # no replace outlier
  lower <- qnt[1] - 3 * iqr
  upper <- qnt[2] + 3 * iqr
  x[x < lower] <- lower
  x[x > upper] <- upper
  return(x)
}

data[vars_to_process] <- lapply(data[vars_to_process], replace_outliers)
data[vars_to_process] <- scale(data[vars_to_process])


# SEM
model <- '
  # measurement model
  ecg =~ HR_mean 
        + HRV_pNN20
  eda =~ SCL_mean
        + SCR_amplitude_mean
  resp =~ RR_mean
  driving =~ 
            SDOffsetfromLaneCenter
            + MaxSteeringAngleVelocity
  visual =~ saccade_frequency
            + mean_fixation_duration
            + fixation_frequency
            + ANNI
            
  # structural model
  ecg              ~ Temp_Cool  + Temp_Warm
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


fit <- sem(model, data = data, 
           # estimator = "MLR",
           missing = "ml",
           cluster = "ID")
summary(fit, fit.measures = TRUE, standardized = TRUE)
fitMeasures(fit, c("chisq", "df", "pvalue", "gfi", "cfi", "nfi", "agfi", "rmr", "srmr", "rmsea"))
# semPaths(fit, what = "est", layout =  "tree", edge.label.cex = 0.8)

# 第一步：提取测量模型部分
lambda_table <- parameterEstimates(fit, standardized = TRUE) %>%
  filter(op == "=~") 

# 第二步：找出有两个及以上观测指标的潜变量
valid_latents <- lambda_table %>%
  group_by(lhs) %>%
  summarise(n_indicators = n()) %>%
  filter(n_indicators >= 2) %>%
  pull(lhs)

# 第三步：提取标准化loading
loadings_std <- parameterEstimates(fit, standardized = TRUE) %>%
  filter(op == "=~", lhs %in% valid_latents) %>%
  select(lhs, rhs, std.all)

# 第四步：计算CR和AVE
cr_ave_table <- loadings_std %>%
  group_by(lhs) %>%
  summarise(
    CR = (sum(std.all))^2 / ((sum(std.all))^2 + sum(1 - std.all^2)),
    AVE = mean(std.all^2)
  )

# 第五步：查看结果
print(cr_ave_table)