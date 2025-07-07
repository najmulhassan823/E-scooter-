# Install packages
install.packages("readxl", "seminr")

# Loading the library
library(seminr)
library(readxl)

# Load and prepare your data
data <- read_excel("final survey data.xlsx")
# To inspect data
head(data)

# Define the Measurement Model
measurement_model <- constructs(
  composite("AT", c("AT1", "AT2", "AT3", "AT4", "AT5")),
  composite("SN", c("SN1", "SN2", "SN3")),
  composite("PBC", c("PBC1", "PBC2", "PBC3", "PBC4")),
  composite("EC", c("EC1", "EC2", "EC3", "EC4")),
  composite("PE", c("PE1", "PE2", "PE3", "PE4")),
  composite("EE", c("EE1", "EE2", "EE3", "EE4")),
  composite("FC", c("FC1", "FC2", "FC3", "FC4")),
  composite("HM", c("HM1", "HM2", "HM3")),
  composite("PR", c("PR1", "PR2", "PR3")),
  composite("H", c("H1", "H2", "H3", "H4")),
  composite("PV", c("PV1", "PV2", "PV3", "PV4")),
  composite("PS", c("PS1", "PS2", "PS3", "PS4")), 
  composite("RE", c("RE1", "RE2", "RE3")),
  composite("SE", c("SE1", "SE2", "SE3", "SE4")),
  composite("BI", c ("BI1", "BI2", "BI3", "BI4", "BI5"))
)

# Define structural model
structural_model <- relationships(
  paths (from = "PV", to = "AT"),
  paths (from = "PS", to = "AT"),
  paths (from = "RE", to = "PBC"),
  paths (from = "SE", to = "PBC"),
  paths (from = "AT", to = 'BI'),
  paths (from = "PBC", to = "BI"),
  paths (from = "SN", to = "BI"),
  paths (from = "EC", to = "AT"),
  paths (from = "PE", to = "AT"),
  paths (from = "EE", to = "AT"),
  paths (from = "FC", to = "PBC"),
  paths (from = "HM", to = "AT"),
  paths (from = "PR", to = "AT"),
  paths (from = "H", to = "AT")
)

# Estimate the model
simple_model <- estimate_pls(
  data = data,
  measurement_model = measurement_model,
  structural_model = structural_model
)

# For printing all rows
options(max.print = 999999)

# Print the summary to get measurement model factor loading, reliability and validity measures.
summary_model <-summary(simple_model)
summary_model
summary_model$loadings


# For discriminant validity
summary_model$validity$fl_criteria
summary_model$validity$htmt
summary_model$validity$cross_loadings


#Bootstrap the model to get the which latent variable is significant, R-square values and path co-efficient
boot_model <-bootstrap_model(
  seminr_model= simple_model,
  nboot= 1000,
  cores= NULL,
  seed=123
)
# Print the summary of boot strap model
summary_boot <-summary(boot_model) 
summary_boot
summary_boot$bootstrapped_paths
summary_boot$bootstrapped_loadings

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Measurement invariance test (MIT) using MICOM
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Split by gender
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gender groups

male_data <- subset(data, Gender == "Male")
female_data <- subset(data, Gender == "Female")

# Check the first few rows of the separated dataset
head(male_data)
head(female_data)


# Estimate the model for Gender
model_male <- estimate_pls(
  data = male_data,
  measurement_model = measurement_model,
  structural_model = structural_model
)

model_female <- estimate_pls(
  data = female_data,
  measurement_model = measurement_model,
  structural_model = structural_model
)

# Convert matrices to data frames
male_scores <- as.data.frame(model_male$construct_scores)
female_scores <- as.data.frame(model_female$construct_scores)

# Identify common constructs
constructs <- intersect(colnames(male_scores), colnames(female_scores))

# Initialize results table
results <- data.frame(
  Construct = constructs,
  Configural_Invariance = rep("Yes", length(constructs)),
  Original_Correlation = NA,
  Quantile_5pct = NA,
  Corr_p_value = NA,
  Partial_MI = NA,
  Mean_Diff = NA,
  Mean_p = NA,
  Mean_Equal = NA,
  Var_Diff = NA,
  Var_p = NA,
  Var_Equal = NA,
  Full_MI = NA
)

# MICOM loop
set.seed(123)  # for reproducibility
n_boot <- 500  # number of bootstrap samples

for (i in seq_along(constructs)) {
  con <- constructs[i]
  
  male_vec <- male_scores[[con]]
  female_vec <- female_scores[[con]]
  
  n_common <- min(length(male_vec), length(female_vec))
  male_sample <- male_vec[1:n_common]
  female_sample <- female_vec[1:n_common]
  
  # Compositional Invariance
  cor_val <- abs(cor(male_sample, female_sample, method = "pearson", use = "complete.obs"))
  boot_corrs <- replicate(n_boot, {
    abs(cor(
      sample(male_sample, replace = TRUE),
      sample(female_sample, replace = TRUE),
      method = "pearson"
    ))
  })
  quant_5pct <- quantile(boot_corrs, 0.05)
  p_val_corr <- mean(boot_corrs >= cor_val)
  partial_MI <- ifelse(cor_val >= quant_5pct, "Yes", "No")
  
  # Equality of Means using matched-length samples
  ttest <- t.test(male_sample, female_sample, var.equal = FALSE, paired = FALSE)
  mean_diff <- round(mean(male_sample, na.rm = TRUE) - mean(female_sample, na.rm = TRUE), 3)
  mean_equal <- ifelse(ttest$p.value > 0.05, "Yes", "No")
  
  # Equality of Variances using matched-length samples
  var_diff <- round(var(male_sample, na.rm = TRUE) - var(female_sample, na.rm = TRUE), 3)
  vartest <- var.test(male_sample, female_sample)
  var_equal <- ifelse(vartest$p.value > 0.05, "Yes", "No")
  
  
  # Now Full MI after defining mean_equal and var_equal
  full_MI <- ifelse(partial_MI == "Yes" & mean_equal == "Yes" & var_equal == "Yes", "Yes", "No")
  
  # Store results
  results$Original_Correlation[i] <- round(cor_val, 3)
  results$Quantile_5pct[i] <- round(quant_5pct, 3)
  results$Corr_p_value[i] <- round(p_val_corr, 4)
  results$Partial_MI[i] <- partial_MI
  results$Mean_Diff[i] <- mean_diff
  results$Mean_p[i] <- round(ttest$p.value, 4)
  results$Mean_Equal[i] <- mean_equal
  results$Var_Diff[i] <- var_diff
  results$Var_p[i] <- round(vartest$p.value, 4)
  results$Var_Equal[i] <- var_equal
  results$Full_MI[i] <- full_MI
}

# Display final MICOM results
print(results)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Split by Age
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Age groups
age_below_and_equal_30 <- subset(data, Ageyears <= 30)
age_above_30 <- subset(data, Ageyears > 30)

# Check the first few rows of each group
head(age_below_and_equal_30)
head(age_above_30)

# Estimate the model for Age groups
model_age_young <- estimate_pls(
  data = age_below_and_equal_30,
  measurement_model = measurement_model,
  structural_model = structural_model
)

model_age_old <- estimate_pls(
  data = age_above_30,
  measurement_model = measurement_model,
  structural_model = structural_model
)

# Convert construct score matrices to data frames
young_scores <- as.data.frame(model_age_young$construct_scores)
old_scores   <- as.data.frame(model_age_old$construct_scores)

# Identify common constructs
constructs <- intersect(colnames(young_scores), colnames(old_scores))

# Initialize results table
results <- data.frame(
  Construct = constructs,
  Configural_Invariance = rep("Yes", length(constructs)),
  Original_Correlation = NA,
  Quantile_5pct = NA,
  Corr_p_value = NA,
  Partial_MI = NA,
  Mean_Diff = NA,
  Mean_p = NA,
  Mean_Equal = NA,
  Var_Diff = NA,
  Var_p = NA,
  Var_Equal = NA,
  Full_MI = NA
)

# MICOM loop
set.seed(123)  # for reproducibility
n_boot <- 500  # number of bootstrap samples

for (i in seq_along(constructs)) {
  con <- constructs[i]
  
  young_vec <- young_scores[[con]]
  old_vec   <- old_scores[[con]]
  
  n_common <- min(length(young_vec), length(old_vec))
  young_sample <- young_vec[1:n_common]
  old_sample   <- old_vec[1:n_common]
  
  # Compositional Invariance
  cor_val <- abs(cor(young_sample, old_sample, method = "pearson", use = "complete.obs"))
  boot_corrs <- replicate(n_boot, {
    abs(cor(
      sample(young_sample, replace = TRUE),
      sample(old_sample, replace = TRUE),
      method = "pearson"
    ))
  })
  quant_5pct <- quantile(boot_corrs, 0.05)
  p_val_corr <- mean(boot_corrs >= cor_val)
  partial_MI <- ifelse(cor_val >= quant_5pct, "Yes", "No")
  
  # Equality of Means using matched-length samples
  ttest <- t.test(young_sample, old_sample, var.equal = FALSE, paired = FALSE)
  mean_diff <- round(mean(young_sample, na.rm = TRUE) - mean(old_sample, na.rm = TRUE), 3)
  mean_equal <- ifelse(ttest$p.value > 0.05, "Yes", "No")
  
  # Equality of Variances using matched-length samples
  var_diff <- round(var(young_sample, na.rm = TRUE) - var(old_sample, na.rm = TRUE), 3)
  vartest <- var.test(young_sample, old_sample)
  var_equal <- ifelse(vartest$p.value > 0.05, "Yes", "No")
  
  # Full Measurement Invariance
  full_MI <- ifelse(partial_MI == "Yes" & mean_equal == "Yes" & var_equal == "Yes", "Yes", "No")
  
  # Store results
  results$Original_Correlation[i] <- round(cor_val, 3)
  results$Quantile_5pct[i] <- round(quant_5pct, 3)
  results$Corr_p_value[i] <- round(p_val_corr, 4)
  results$Partial_MI[i] <- partial_MI
  results$Mean_Diff[i] <- mean_diff
  results$Mean_p[i] <- round(ttest$p.value, 4)
  results$Mean_Equal[i] <- mean_equal
  results$Var_Diff[i] <- var_diff
  results$Var_p[i] <- round(vartest$p.value, 4)
  results$Var_Equal[i] <- var_equal
  results$Full_MI[i] <- full_MI
}

# Display final MICOM results
print(results)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Split by Monthly Income
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Monthly income group
low_income <- subset(data, MonthlyincomePKR == "Less than 37,000")
middle_high_income <- subset(data, MonthlyincomePKR %in% c("37,000 to 49,999", "50,000 to 64,999", "65,000 to 79,999", "80,000 to 94,999", "95,000 to 109,999", "110,000 or above")) 

# Check the first few rows of each group
head(low_income)
head(middle_high_income)

model_low_income <- estimate_pls(
  data = low_income,
  measurement_model = measurement_model,
  structural_model = structural_model
)

model_mid_high_income <- estimate_pls(
  data = middle_high_income,
  measurement_model = measurement_model,
  structural_model = structural_model
)

# Convert matrices to data frames
low_scores <- as.data.frame(model_low_income$construct_scores)
high_scores <- as.data.frame(model_mid_high_income$construct_scores)

# Identify common constructs
constructs <- intersect(colnames(low_scores), colnames(high_scores))

# Initialize results table
results <- data.frame(
  Construct = constructs,
  Configural_Invariance = rep("Yes", length(constructs)),
  Original_Correlation = NA,
  Quantile_5pct = NA,
  Corr_p_value = NA,
  Partial_MI = NA,
  Mean_Diff = NA,
  Mean_p = NA,
  Mean_Equal = NA,
  Var_Diff = NA,
  Var_p = NA,
  Var_Equal = NA,
  Full_MI = NA
)

# MICOM loop
set.seed(123)
n_boot <- 500

for (i in seq_along(constructs)) {
  con <- constructs[i]
  
  low_vec <- low_scores[[con]]
  high_vec <- high_scores[[con]]
  
  n_common <- min(length(low_vec), length(high_vec))
  low_sample <- low_vec[1:n_common]
  high_sample <- high_vec[1:n_common]
  
  # Step 2: Compositional Invariance
  cor_val <- abs(cor(low_sample, high_sample, method = "pearson", use = "complete.obs"))
  boot_corrs <- replicate(n_boot, {
    abs(cor(
      sample(low_sample, replace = TRUE),
      sample(high_sample, replace = TRUE),
      method = "pearson"
    ))
  })
  quant_5pct <- quantile(boot_corrs, 0.05)
  p_val_corr <- mean(boot_corrs >= cor_val)
  partial_MI <- ifelse(cor_val >= quant_5pct, "Yes", "No")
  
  # Means
  ttest <- t.test(low_sample, high_sample, var.equal = FALSE, paired = FALSE)
  mean_diff <- round(mean(low_sample, na.rm = TRUE) - mean(high_sample, na.rm = TRUE), 3)
  mean_equal <- ifelse(ttest$p.value > 0.05, "Yes", "No")
  
  # Variances
  var_diff <- round(var(low_sample, na.rm = TRUE) - var(high_sample, na.rm = TRUE), 3)
  vartest <- var.test(low_sample, high_sample)
  var_equal <- ifelse(vartest$p.value > 0.05, "Yes", "No")
  
  # Full MI
  full_MI <- ifelse(partial_MI == "Yes" & mean_equal == "Yes" & var_equal == "Yes", "Yes", "No")
  
  # Store results
  results$Original_Correlation[i] <- round(cor_val, 3)
  results$Quantile_5pct[i] <- round(quant_5pct, 3)
  results$Corr_p_value[i] <- round(p_val_corr, 4)
  results$Partial_MI[i] <- partial_MI
  results$Mean_Diff[i] <- mean_diff
  results$Mean_p[i] <- round(ttest$p.value, 4)
  results$Mean_Equal[i] <- mean_equal
  results$Var_Diff[i] <- var_diff
  results$Var_p[i] <- round(vartest$p.value, 4)
  results$Var_Equal[i] <- var_equal
  results$Full_MI[i] <- full_MI
}

# View results
print(results)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Split by Transit use frequency
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Public transit use frequency
Occasional_user <- subset(data, Frequencyofpublictransituse %in% c ("Less than once per month", "1 or 2 days per month","1 day per week", "2 or 3 days per week"))
Frequent_user <- subset(data, Frequencyofpublictransituse %in% c ("4 or 5 days per week", "6 or 7 days per week"))

#check the first few rows of each group
head(Occasional_user)
head(Frequent_user)

model_occasional <- estimate_pls(
  data = Occasional_user,
  measurement_model = measurement_model,
  structural_model = structural_model
)

model_frequent <- estimate_pls(
  data = Frequent_user,
  measurement_model = measurement_model,
  structural_model = structural_model
)

# Convert matrices to data frames
occasional_scores <- as.data.frame(model_occasional$construct_scores)
frequent_scores <- as.data.frame(model_frequent$construct_scores)

# Identify common constructs
constructs <- intersect(colnames(occasional_scores), colnames(frequent_scores))

# Initialize results table
results <- data.frame(
  Construct = constructs,
  Configural_Invariance = rep("Yes", length(constructs)),
  Original_Correlation = NA,
  Quantile_5pct = NA,
  Corr_p_value = NA,
  Partial_MI = NA,
  Mean_Diff = NA,
  Mean_p = NA,
  Mean_Equal = NA,
  Var_Diff = NA,
  Var_p = NA,
  Var_Equal = NA,
  Full_MI = NA
)

# MICOM loop
set.seed(123)
n_boot <- 500

for (i in seq_along(constructs)) {
  con <- constructs[i]
  
  occasional_vec <- occasional_scores[[con]]
  frequent_vec <- frequent_scores[[con]]
  
  n_common <- min(length(occasional_vec), length(frequent_vec))
  occasional_sample <- occasional_vec[1:n_common]
  frequent_sample <- frequent_vec[1:n_common]
  
  # Compositional Invariance
  cor_val <- abs(cor(occasional_sample, frequent_sample, method = "pearson", use = "complete.obs"))
  boot_corrs <- replicate(n_boot, {
    abs(cor(
      sample(occasional_sample, replace = TRUE),
      sample(frequent_sample, replace = TRUE),
      method = "pearson"
    ))
  })
  quant_5pct <- quantile(boot_corrs, 0.05)
  p_val_corr <- mean(boot_corrs >= cor_val)
  partial_MI <- ifelse(cor_val >= quant_5pct, "Yes", "No")
  
  # Means
  ttest <- t.test(occasional_sample, frequent_sample, var.equal = FALSE, paired = FALSE)
  mean_diff <- round(mean(occasional_sample, na.rm = TRUE) - mean(frequent_sample, na.rm = TRUE), 3)
  mean_equal <- ifelse(ttest$p.value > 0.05, "Yes", "No")
  
  # Variances
  var_diff <- round(var(occasional_sample, na.rm = TRUE) - var(frequent_sample, na.rm = TRUE), 3)
  vartest <- var.test(occasional_sample, frequent_sample)
  var_equal <- ifelse(vartest$p.value > 0.05, "Yes", "No")
  
  # Full MI
  full_MI <- ifelse(partial_MI == "Yes" & mean_equal == "Yes" & var_equal == "Yes", "Yes", "No")
  
  # Store results
  results$Original_Correlation[i] <- round(cor_val, 3)
  results$Quantile_5pct[i] <- round(quant_5pct, 3)
  results$Corr_p_value[i] <- round(p_val_corr, 4)
  results$Partial_MI[i] <- partial_MI
  results$Mean_Diff[i] <- mean_diff
  results$Mean_p[i] <- round(ttest$p.value, 4)
  results$Mean_Equal[i] <- mean_equal
  results$Var_Diff[i] <- var_diff
  results$Var_p[i] <- round(vartest$p.value, 4)
  results$Var_Equal[i] <- var_equal
  results$Full_MI[i] <- full_MI
}

# Show Results
print(results)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Multi-Group Analysis (MGA) for groups
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# MGA for Age

# Split your data
young_data <- subset(data, Ageyears <= 30)
old_data <- subset(data, Ageyears > 30)

# Estimate model for young group
model_young <- estimate_pls(
  data = young_data,
  measurement_model = measurement_model,
  structural_model = structural_model
)

boot_young <- bootstrap_model(
  seminr_model = model_young,
  nboot = 1000,
  cores = NULL,
  seed = 123
)

# Get the bootstrapped path coefficients for young group
young_paths <- summary(boot_young)$bootstrapped_paths


# Estimate model for old group
model_old <- estimate_pls(
  data = old_data,
  measurement_model = measurement_model,
  structural_model = structural_model
)

boot_old <- bootstrap_model(
  seminr_model = model_old,
  nboot = 1000,
  cores = NULL,
  seed = 123
)

# Get the bootstrapped path coefficients for old group
old_paths <- summary(boot_old)$bootstrapped_paths


# Get common paths
common_paths <- intersect(rownames(young_paths), rownames(old_paths))


# Initialize results
results <- data.frame(
  Path = character(),
  Young_Mean = numeric(),
  Old_Mean = numeric(),
  Difference = numeric(),
  Z_value = numeric(),
  P_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop to compare each path using Z-test
for (path in common_paths) {
  b1 <- young_paths[path, "Bootstrap Mean"]
  b2 <- old_paths[path, "Bootstrap Mean"]
  se1 <- young_paths[path, "Bootstrap SD"]
  se2 <- old_paths[path, "Bootstrap SD"]
  
  diff <- b1- b2
  z <- (b1 - b2) / sqrt(se1^2 + se2^2)
  p <- 2 * (1 - pnorm(abs(z)))
  sig <- ifelse(p < 0.05, "Yes", "No")
  
  results <- rbind(results, data.frame(
    Path = path,
    Young_Mean = round(b1, 3),
    Old_Mean = round(b2, 3),
    Difference = round(diff, 3),
    Z_value = round(z, 3),
    P_value = round(p, 3),
    Significant = sig
  ))
}

# Show results
print(results)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MGA for Income

# Split your data
low_income <- subset(data, MonthlyincomePKR == "Less than 37,000")
mid_high_income <- subset(data, MonthlyincomePKR %in% c(
  "37,000 to 49,999", "50,000 to 64,999", "65,000 to 79,999",
  "80,000 to 94,999", "95,000 to 109,999", "110,000 or above"))


# Estimate model for low income group
model_low <- estimate_pls(
  data = low_income,
  measurement_model = measurement_model,
  structural_model = structural_model
)

boot_low <- bootstrap_model(
  seminr_model = model_low,
  nboot = 1000,
  cores = NULL,
  seed = 123
)

low_paths <- summary(boot_low)$bootstrapped_paths

# Estimate model for mid-high income group
model_mid_high <- estimate_pls(
  data = mid_high_income,
  measurement_model = measurement_model,
  structural_model = structural_model
)

boot_mid_high <- bootstrap_model(
  seminr_model = model_mid_high,
  nboot = 1000,
  cores = NULL,
  seed = 123
)

mid_high_paths <- summary(boot_mid_high)$bootstrapped_paths

#  Get common paths
common_paths_income <- intersect(rownames(low_paths), rownames(mid_high_paths))

# Initialize results table
income_results <- data.frame(
  Path = character(),
  Low_Income_Mean = numeric(),
  Mid_High_Income_Mean = numeric(),
  Difference = numeric(),
  Z_value = numeric(),
  P_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop through and calculate Z-test for each path
for (path in common_paths_income) {
  b1 <- low_paths[path, "Bootstrap Mean"]
  b2 <- mid_high_paths[path, "Bootstrap Mean"]
  se1 <- low_paths[path, "Bootstrap SD"]
  se2 <- mid_high_paths[path, "Bootstrap SD"]
  
  diff <- b1 - b2
  z <- diff / sqrt(se1^2 + se2^2)
  p <- 2 * (1 - pnorm(abs(z)))
  sig <- ifelse(p < 0.05, "Yes", "No")
  
  income_results <- rbind(income_results, data.frame(
    Path = path,
    Low_Income_Mean = round(b1, 3),
    Mid_High_Income_Mean = round(b2, 3),
    Difference = round(diff, 3),
    Z_value = round(z, 3),
    P_value = round(p, 3),
    Significant = sig
  ))
}

# View the result
print(income_results)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MGA for Transit use frequency

# Split your data
occasional_user <- subset(data, Frequencyofpublictransituse %in% c(
  "Less than once per month", "1 or 2 days per month", 
  "1 day per week", "2 or 3 days per week"))

frequent_user <- subset(data, Frequencyofpublictransituse %in% c(
  "4 or 5 days per week", "6 or 7 days per week"))


# Estimate model for occasional users
model_occasional <- estimate_pls(
  data = occasional_user,
  measurement_model = measurement_model,
  structural_model = structural_model
)

boot_occasional <- bootstrap_model(
  seminr_model = model_occasional,
  nboot = 1000,
  cores = NULL,
  seed = 123
)

occasional_paths <- summary(boot_occasional)$bootstrapped_paths

# Estimate model for frequent users
model_frequent <- estimate_pls(
  data = frequent_user,
  measurement_model = measurement_model,
  structural_model = structural_model
)

boot_frequent <- bootstrap_model(
  seminr_model = model_frequent,
  nboot = 1000,
  cores = NULL,
  seed = 123
)

frequent_paths <- summary(boot_frequent)$bootstrapped_paths

# Find common paths between both groups
common_paths_freq <- intersect(rownames(occasional_paths), rownames(frequent_paths))

# Initialize results table
frequency_results <- data.frame(
  Path = character(),
  Occasional_User_Mean = numeric(),
  Frequent_User_Mean = numeric(),
  Difference = numeric(),
  Z_value = numeric(),
  P_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop through each path and perform Z-test
for (path in common_paths_freq) {
  b1 <- occasional_paths[path, "Bootstrap Mean"]
  b2 <- frequent_paths[path, "Bootstrap Mean"]
  se1 <- occasional_paths[path, "Bootstrap SD"]
  se2 <- frequent_paths[path, "Bootstrap SD"]
  
  diff <- b1 - b2
  z <- diff / sqrt(se1^2 + se2^2)
  p <- 2 * (1 - pnorm(abs(z)))
  sig <- ifelse(p < 0.05, "Yes", "No")
  
  frequency_results <- rbind(frequency_results, data.frame(
    Path = path,
    Occasional_User_Mean = round(b1, 3),
    Frequent_User_Mean = round(b2, 3),
    Difference = round(diff, 3),
    Z_value = round(z, 3),
    P_value = round(p, 3),
    Significant = sig
  ))
}

# Print the results
print(frequency_results)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Importance Performance Map Analysis (IPMA)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract importance values (Total path) form your boot mode, saved in summary

summary_boot$bootstrapped_total_paths

# separate the original_Est column,to get the paths ending on BI to get total effect of each construct on BI
importance <- summary_boot$bootstrapped_total_paths[grepl("->  BI$", rownames(summary_boot$bootstrapped_total_paths)), "Original Est."]

# Convert it into data frame
importance <- data.frame(
  Construct = sub(" ->  BI", "", rownames(summary_boot$bootstrapped_total_paths)[grepl("->  BI$", rownames(summary_boot$bootstrapped_total_paths))]),
  Importance = importance
)
print(importance)

# Extract performance values

# Get the latent variable scores from the model
lv_scores <- simple_model$construct_scores

# Remove the target construct (BI) from performance calculation
lv_scores <- lv_scores[, colnames(lv_scores) != "BI"]

# Min-max scale each construct's scores individually to 0-100
scale_0_100 <- function(x) {
  if (min(x, na.rm = TRUE) == max(x, na.rm = TRUE)) {
    return(rep(50, length(x)))  # avoid divide by zero; assume neutral
  } else {
    return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * 100)
  }
}

lv_scores_scaled <- as.data.frame(apply(lv_scores, 2, scale_0_100))

# Take the mean for each construct's scaled values
performance <- data.frame(
  Construct = colnames(lv_scores_scaled),
  Performance = round(colMeans(lv_scores_scaled, na.rm = TRUE), 8)
)

# Print final performance table
print(performance)

importance$Construct <- trimws(as.character(importance$Construct))
performance$Construct <- trimws(as.character(performance$Construct))

# Merge again
ipma_table <- merge(importance, performance, by = "Construct", all.x = TRUE)
print(ipma_table)

