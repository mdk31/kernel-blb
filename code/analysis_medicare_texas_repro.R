rm(list=ls())
library(pbapply)
library(reticulate)
library(data.table)
library(devtools)
library(fst)
library(tidyverse)
library(MatchIt)
library(cobalt)

# Set Working Directory 
source('code/helper_functions.R')
source('scale_function.R')

# Kernel specific parts
use_virtualenv("r-reticulate", required = TRUE)
#py_install("numpy<2.0")
#py_install("scikit-learn")
np <- import("numpy")
source_python("code/gp_simu_gate.py")

# Load Data
aggregate_data <- as.data.frame(read_fst("/n/dominici_nsaph_l3/Lab/projects/fbargaglistoffi_casual_rule_ensamble/aggregate_medicare_data_2010to2016.fst"))


# Filter for region of interest: Texas
TX_data <- aggregate_data %>%
  filter(statecode == "TX") %>%
  filter(popdensity != 0)

# Texas Data
dim(TX_data)

# Select confounders and treatment
TX_data <- TX_data %>%
  select(pm25_12, dead_in_5, sex, race, age, dual,
         mean_bmi, smoke_rate, hispanic, pct_blk,
         medhouseholdincome, medianhousevalue, poverty, education, popdensity,
         pct_owner_occ, summer_tmmx, winter_tmmx, summer_rmax, winter_rmax) %>%
  rename(pct_hispanic = hispanic)

# Perform log transformation on popdensity
TX_data$popdensity <- log(TX_data$popdensity)

# Extract data
X <- TX_data %>% select(-dead_in_5, -pm25_12) 
covariates <- names(X)
TX_data$outcome <- as.numeric(TX_data$dead_in_5)
TX_data$treat <- TX_data$pm25_12 

# Apply standardization
result <- standardize_mixed_data(X, max_unique_for_discrete = 10)

# Access standardized data
TX_data_standardized <- TX_data
TX_data_standardized[names(X)] <- result$data

# Data Balancing (1 to 1)

# Check current balance
table(TX_data_standardized$treat)
prop.table(table(TX_data_standardized$treat))

# Downsample to balance the treatment variable
# set.seed(123)  # for reproducibility
# 
# # Split data by treatment group
# treat_0 <- TX_data_standardized[TX_data_standardized$treat == 0, ]
# treat_1 <- TX_data_standardized[TX_data_standardized$treat == 1, ]
# 
# # Find the minority class size
# min_size <- min(nrow(treat_0), nrow(treat_1))
# 
# # Randomly sample from the majority class to match minority class size
# if (nrow(treat_0) > nrow(treat_1)) {
#   # Downsample treat_0
#   treat_0_balanced <- treat_0[sample(nrow(treat_0), min_size), ] 
#   TX_data_standardized_balanced <- rbind(treat_0_balanced, treat_1)
# } else {
#   # Downsample treat_1
#   treat_1_balanced <- treat_1[sample(nrow(treat_1), min_size), ]
#   TX_data_standardized_balanced <- rbind(treat_0, treat_1_balanced)
# }

# Shuffle the balanced dataset
TX_data_standardized_balanced <- TX_data_standardized_balanced[sample(nrow(TX_data_standardized_balanced)), ]

# Check new balance
table(TX_data_standardized_balanced$treat)
prop.table(table(TX_data_standardized_balanced$treat))

### Causal BLB
# Set Up parameters
degree1 <- 1
degree2 <- 1
k1 <- "poly"
k2 <- "poly"
operator <- "single"
penal <- log(2)
subsets <- 10
gamma <- calculate_gamma(nrow(TX_data_standardized_balanced), subsets)
b <- floor(nrow(TX_data_standardized_balanced)^gamma)
B <- 100

# Run BLB AIPW
system.time(
  aipw_blb_results_1 <-  causal_blb_aipw(data = TX_data_standardized_balanced, y = 'outcome', Tr = 'treat', confounders = covariates,
                                         b = b, degree1, degree2, k1, k2, operator, penal, subsets = subsets)
)

aipw_blb_results_1

# AIPW Coefficient to Probability Interpretation with Confidence Intervals

# Extract AIPW results (replace with your actual AIPW output)
aipw_coefficient <- aipw_blb_results_1$estim
aipw_se <- aipw_blb_results_1$se

# Calculate baseline mortality rate (control group)
baseline_mortality <- mean(TX_data$outcome[TX_data$treat == 0], na.rm = TRUE)

# Calculate treated group mortality (baseline + AIPW effect)
treated_mortality <- baseline_mortality + aipw_coefficient

# 95% confidence intervals for AIPW coefficient
aipw_ci_lower <- aipw_blb_results_1$lower_ci
aipw_ci_upper <- aipw_blb_results_1$upper_ci

# Confidence intervals for treated mortality probability
treated_ci_lower <- baseline_mortality + aipw_ci_lower
treated_ci_upper <- baseline_mortality + aipw_ci_upper


cat("AIPW coefficient (marginal effect):", round(aipw_coefficient, 6), "\n")
cat("Standard error:", round(aipw_se, 6), "\n")
cat("95% CI for marginal effect: [", round(aipw_ci_lower, 6), ",", round(aipw_ci_upper, 6), "]\n")
cat("In percentage points: [", round(aipw_ci_lower * 100, 3), ",", round(aipw_ci_upper * 100, 3), "]\n\n")

cat("Baseline mortality (PM2.5 ≤ 12 μg/m³):", round(baseline_mortality * 100, 2), "%\n")
cat("Treated mortality (PM2.5 > 12 μg/m³):", round(treated_mortality * 100, 2), "%\n")
cat("95% CI for treated mortality: [", round(treated_ci_lower * 100, 2), "%,", round(treated_ci_upper * 100, 2), "%]\n\n")

cat("Being above NAAQS threshold changes 5-year mortality by:", round(aipw_coefficient * 100, 3), "percentage points\n")
cat("95% CI for causal effect: [", round(aipw_ci_lower * 100, 3), ",", round(aipw_ci_upper * 100, 3), "] percentage points\n\n")


cat("PM2.5 above NAAQS INCREASES 5-year mortality risk\n")
cat("For every 1000 people exposed above threshold:\n")
cat("  - Expected additional deaths:", round(aipw_coefficient * 1000, 1), "\n")
cat("  - 95% CI: [", round(aipw_ci_lower * 1000, 1), ",", round(aipw_ci_upper * 1000, 1), "] additional deaths\n")

# Calculate Texas population estimates
tx_total_population <- nrow(TX_data)  # Total sample size
tx_above_naaqs <- sum(TX_data$treat == 1, na.rm = TRUE)  # People above NAAQS
tx_population_estimate <- 31290831  # Texas population as of July 2024 (latest Census estimate)# AIPW Coefficient to Probability Interpretation with Confidence Intervals

cat("For Texas population (~", round(tx_population_estimate/1000000, 1), " million people):\n")
# Assume proportion above NAAQS in sample applies to population
prop_above_naaqs <- tx_above_naaqs / tx_total_population
tx_exposed_population <- tx_population_estimate * prop_above_naaqs

cat("  - People exposed above NAAQS:", round(tx_exposed_population/1000000, 1), "million (", round(prop_above_naaqs * 100, 1), "% of population)\n")
cat("  - Expected additional deaths:", round(aipw_coefficient * tx_exposed_population, 0), "\n")
cat("  - 95% CI: [", round(aipw_ci_lower * tx_exposed_population, 0), ",", round(aipw_ci_upper * tx_exposed_population, 0), "] additional deaths\n")


# Create a summary table of key results
results_summary <- data.frame(
  Metric = c("AIPW Coefficient", "Standard Error", "95% CI Lower", "95% CI Upper",
             "Baseline Mortality (%)", "Treated Mortality (%)",
             "Additional Deaths per 1000", "TX Exposed Population (millions)",
             "TX Expected Additional Deaths"),
  Value = c(round(aipw_coefficient, 6), round(aipw_se, 6), 
            round(aipw_ci_lower, 6), round(aipw_ci_upper, 6),
            round(baseline_mortality * 100, 2), round(treated_mortality * 100, 2),
            round(aipw_coefficient * 1000, 1), round(tx_exposed_population/1000000, 1),
            round(aipw_coefficient * tx_exposed_population, 0))
)
results_summary$Value <- round(results_summary$Value, 5)
results_summary

# Save results tValue# Save results to CSV
write.csv(results_summary, file="/n/dominici_nsaph_l3/Lab/projects/fbargaglistoffi_ccit/bootstrap/aipw_results_summary_treated_controls_same_size.csv", row.names = FALSE)

# Data Balancing (1 to 2)

# Downsample to balance the treatment variable
set.seed(123)  # for reproducibility

# Randomly sample from the majority class to match minority class size
if (nrow(treat_0) > nrow(treat_1)) {
  # Downsample treat_0
  treat_0_balanced <- treat_0[sample(nrow(treat_0), min_size * 2), ] 
  TX_data_standardized_balanced <- rbind(treat_0_balanced, treat_1)
} else {
  # Downsample treat_1
  treat_1_balanced <- treat_1[sample(nrow(treat_1), min_size), ]
  TX_data_standardized_balanced <- rbind(treat_0, treat_1_balanced)
}

# Shuffle the balanced dataset
TX_data_standardized_balanced <- TX_data_standardized_balanced[sample(nrow(TX_data_standardized_balanced)), ]

# Check new balance
table(TX_data_standardized_balanced$treat)
prop.table(table(TX_data_standardized_balanced$treat))

### Causal BLB
# Set Up parameters
degree1 <- 1
degree2 <- 1
k1 <- "poly"
k2 <- "poly"
operator <- "single"
penal <- log(2)
subsets <- 10
gamma <- calculate_gamma(nrow(TX_data_standardized_balanced), subsets)
b <- floor(nrow(TX_data_standardized_balanced)^gamma)
B <- 100

# Run BLB AIPW
system.time(
  aipw_blb_results_2 <-  causal_blb_aipw(data = TX_data_standardized_balanced, y = 'outcome', Tr = 'treat', confounders = covariates,
                                         b = b, degree1, degree2, k1, k2, operator, penal, subsets = subsets)
)

aipw_blb_results_2

# AIPW Coefficient to Probability Interpretation with Confidence Intervals

# Extract AIPW results (replace with your actual AIPW output)
aipw_coefficient <- aipw_blb_results_2$estim
aipw_se <- aipw_blb_results_2$se

# Calculate baseline mortality rate (control group)
baseline_mortality <- mean(TX_data$outcome[TX_data$treat == 0], na.rm = TRUE)

# Calculate treated group mortality (baseline + AIPW effect)
treated_mortality <- baseline_mortality + aipw_coefficient

# 95% confidence intervals for AIPW coefficient
aipw_ci_lower <- aipw_blb_results_2$lower_ci
aipw_ci_upper <- aipw_blb_results_2$upper_ci

# Confidence intervals for treated mortality probability
treated_ci_lower <- baseline_mortality + aipw_ci_lower
treated_ci_upper <- baseline_mortality + aipw_ci_upper


cat("AIPW coefficient (marginal effect):", round(aipw_coefficient, 6), "\n")
cat("Standard error:", round(aipw_se, 6), "\n")
cat("95% CI for marginal effect: [", round(aipw_ci_lower, 6), ",", round(aipw_ci_upper, 6), "]\n")
cat("In percentage points: [", round(aipw_ci_lower * 100, 3), ",", round(aipw_ci_upper * 100, 3), "]\n\n")

cat("Baseline mortality (PM2.5 ≤ 12 μg/m³):", round(baseline_mortality * 100, 2), "%\n")
cat("Treated mortality (PM2.5 > 12 μg/m³):", round(treated_mortality * 100, 2), "%\n")
cat("95% CI for treated mortality: [", round(treated_ci_lower * 100, 2), "%,", round(treated_ci_upper * 100, 2), "%]\n\n")

cat("Being above NAAQS threshold changes 5-year mortality by:", round(aipw_coefficient * 100, 3), "percentage points\n")
cat("95% CI for causal effect: [", round(aipw_ci_lower * 100, 3), ",", round(aipw_ci_upper * 100, 3), "] percentage points\n\n")


cat("PM2.5 above NAAQS INCREASES 5-year mortality risk\n")
cat("For every 1000 people exposed above threshold:\n")
cat("  - Expected additional deaths:", round(aipw_coefficient * 1000, 1), "\n")
cat("  - 95% CI: [", round(aipw_ci_lower * 1000, 1), ",", round(aipw_ci_upper * 1000, 1), "] additional deaths\n")

# Calculate Texas population estimates
tx_total_population <- nrow(TX_data)  # Total sample size
tx_above_naaqs <- sum(TX_data$treat == 1, na.rm = TRUE)  # People above NAAQS
tx_population_estimate <- 31290831  # Texas population as of July 2024 (latest Census estimate)# AIPW Coefficient to Probability Interpretation with Confidence Intervals

cat("For Texas population (~", round(tx_population_estimate/1000000, 1), " million people):\n")
# Assume proportion above NAAQS in sample applies to population
prop_above_naaqs <- tx_above_naaqs / tx_total_population
tx_exposed_population <- tx_population_estimate * prop_above_naaqs

cat("  - People exposed above NAAQS:", round(tx_exposed_population/1000000, 1), "million (", round(prop_above_naaqs * 100, 1), "% of population)\n")
cat("  - Expected additional deaths:", round(aipw_coefficient * tx_exposed_population, 0), "\n")
cat("  - 95% CI: [", round(aipw_ci_lower * tx_exposed_population, 0), ",", round(aipw_ci_upper * tx_exposed_population, 0), "] additional deaths\n")


# Create a summary table of key results
results_summary <- data.frame(
  Metric = c("AIPW Coefficient", "Standard Error", "95% CI Lower", "95% CI Upper",
             "Baseline Mortality (%)", "Treated Mortality (%)",
             "Additional Deaths per 1000", "TX Exposed Population (millions)",
             "TX Expected Additional Deaths"),
  Value = c(round(aipw_coefficient, 6), round(aipw_se, 6), 
            round(aipw_ci_lower, 6), round(aipw_ci_upper, 6),
            round(baseline_mortality * 100, 2), round(treated_mortality * 100, 2),
            round(aipw_coefficient * 1000, 1), round(tx_exposed_population/1000000, 1),
            round(aipw_coefficient * tx_exposed_population, 0))
)
results_summary$Value <- round(results_summary$Value, 5)
results_summary

# Save results tValue# Save results to CSV
write.csv(results_summary, file="/n/dominici_nsaph_l3/Lab/projects/fbargaglistoffi_ccit/bootstrap/aipw_results_summary_treated_controls_1_to_2_size.csv", row.names = FALSE)

