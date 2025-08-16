rm(list=ls())
library(pbapply)
library(reticulate)
library(data.table)
library(devtools)
library(fst)
library(dplyr)

# DGP
synth_tx_dataset <- function(n, seed = 123) {
  stopifnot(n > 0)
  set.seed(seed)
  
  # helper ---------------------------------------------------------------
  clamp01 <- function(x) pmin(pmax(x, 0), 1)
  rlnorm_tuned <- function(n, meanlog, sdlog, min_val = NULL, max_val = NULL) {
    x <- rlnorm(n, meanlog = meanlog, sdlog = sdlog)
    if (!is.null(min_val)) x[x < min_val] <- min_val
    if (!is.null(max_val)) x[x > max_val] <- max_val
    x
  }
  rtruncnorm <- function(n, mean, sd, lo, hi) {
    x <- rnorm(n, mean, sd)
    x <- pmin(pmax(x, lo), hi)
    x
  }
  
  # latent drivers -------------------------------------------------------
  # Urbanicity / density, socioeconomics, climate heterogeneity
  z_pop    <- rnorm(n, 0, 1)   # higher -> more urban
  z_pov    <- rnorm(n, 0, 1)   # higher -> poorer
  z_clim   <- rnorm(n, 0, 1)   # east/west & gulf gradients
  
  # demographics ---------------------------------------------------------
  sex  <- rbinom(n, 1, 0.45)  # 1 = male, 0 = female
  # Race: 1 White, 2 Black, 3 Other (TX has separate hispanic share below)
  race <- sample.int(3, n, replace = TRUE, prob = c(0.70, 0.12, 0.18))
  
  age  <- rtruncnorm(n, mean = 76, sd = 7, lo = 65, hi = 99)
  
  # population density (persons/km^2), strictly > 0
  popdensity <- rlnorm_tuned(n,
                             meanlog = log(400) - 0.5^2/2 + 0.6 * z_pop,
                             sdlog   = 0.9,
                             min_val = 1)  # never 0 (your code logs this)
  
  # climate --------------------------------------------------------------
  # Max temperature in °C; rmax as % relative humidity
  summer_tmmx <- rtruncnorm(n, mean = 35 + 1.5*z_clim - 0.5*z_pop, sd = 2.5, lo = 28, hi = 43)
  winter_tmmx <- rtruncnorm(n, mean = 16 + 1.0*z_clim - 0.3*z_pop, sd = 3.0, lo = 5,  hi = 25)
  summer_rmax <- rtruncnorm(n, mean = 62 + 6*z_clim - 2*z_pop, sd = 7, lo = 35, hi = 90)
  winter_rmax <- rtruncnorm(n, mean = 68 + 5*z_clim - 1*z_pop, sd = 6, lo = 35, hi = 95)
  
  # socioeconomics -------------------------------------------------------
  # pct_blk increases with urbanicity; hispanic higher overall in TX, mild urban tilt
  pct_blk      <- clamp01(plogis(-1.2 + 0.8*z_pop + rnorm(n, 0, 0.5)))
  hispanic     <- clamp01(plogis( 0.3 + 0.3*z_pop + rnorm(n, 0, 0.6))) # will be renamed pct_hispanic later
  poverty      <- clamp01(plogis(-0.2 + 0.8*z_pov - 0.3*z_pop + rnorm(n, 0, 0.5)))
  education    <- clamp01(plogis( 0.5 - 1.0*z_pov + 0.2*z_pop + rnorm(n, 0, 0.5))) # higher is more educated
  pct_owner_occ<- clamp01(plogis( 0.8 - 0.9*z_pop - 0.5*z_pov + rnorm(n, 0, 0.5)))
  
  # income/house value (lognormal; increase with education & urbanicity)
  medhouseholdincome <- rlnorm_tuned(
    n,
    meanlog = log(60000) - 0.35^2/2 + 0.25*education - 0.25*poverty + 0.2*z_pop + rnorm(n,0,0.05),
    sdlog   = 0.35,
    min_val = 25000, max_val = 160000
  )
  medianhousevalue <- rlnorm_tuned(
    n,
    meanlog = log(220000) - 0.45^2/2 + 0.25*education + 0.35*z_pop - 0.15*poverty + rnorm(n,0,0.06),
    sdlog   = 0.45,
    min_val = 50000, max_val = 800000
  )
  
  dual <- rbinom(
    n, 1,
    clamp01(plogis(-0.5 + 1.1*poverty - 0.6*education - 0.2*z_pop + rnorm(n,0,0.4))))
  
  # health behaviors / clinical ----------------------------------------
  smoke_rate <- clamp01(plogis(-0.6 + 0.9*poverty - 0.3*education - 0.1*z_pop + rnorm(n,0,0.5)))
  mean_bmi   <- rtruncnorm(n,
                           mean = 28.5 + 1.2*poverty - 0.6*education + 0.2*smoke_rate + rnorm(n,0,0.8),
                           sd   = 3.5, lo = 18, hi = 48)
  
  # treatment: PM2.5 (µg/m^3) ------------------------------------------
  pm25_12 <- 8.0 +
    1.0*scale(popdensity)[,1] +
    0.5*scale(winter_rmax)[,1] -
    0.3*scale(summer_tmmx)[,1] +
    0.2*scale(poverty)[,1] +
    rnorm(n, 0, 0.9)
  pm25_12 <- pmin(pmax(pm25_12, 3.0), 18.0)
  
  # outcome: death within 5 years (binary) ------------------------------
  male <- sex
  linpred <- -11.0 +
    0.095*(age - 75) +
    0.35*male +
    0.60*dual +
    1.10*smoke_rate +
    0.55*poverty +
    0.03*(mean_bmi - 28) +
    0.04*(pm25_12 - 8) +
    0.15*pct_blk -
    0.10*education +
    rnorm(n, 0, 0.35)
  
  p_dead <- clamp01(plogis(linpred))
  dead_in_5 <- rbinom(n, 1, p_dead)
  
  # assemble in your expected order -------------------------------------
  out <- data.frame(
    pm25_12 = as.numeric(pm25_12),
    dead_in_5 = as.integer(dead_in_5),
    sex = as.integer(sex),               # 0 = female, 1 = male
    race = as.integer(race),             # 1 White, 2 Black, 3 Other
    age = as.numeric(age),
    dual = as.integer(dual),
    mean_bmi = as.numeric(mean_bmi),
    smoke_rate = as.numeric(smoke_rate),
    hispanic = as.numeric(hispanic),     # your pipeline renames -> pct_hispanic
    pct_blk = as.numeric(pct_blk),
    medhouseholdincome = as.numeric(medhouseholdincome),
    medianhousevalue = as.numeric(medianhousevalue),
    poverty = as.numeric(poverty),
    education = as.numeric(education),
    popdensity = as.numeric(popdensity), # strictly > 0 (you log this next)
    pct_owner_occ = as.numeric(pct_owner_occ),
    summer_tmmx = as.numeric(summer_tmmx),
    winter_tmmx = as.numeric(winter_tmmx),
    summer_rmax = as.numeric(summer_rmax),
    winter_rmax = as.numeric(winter_rmax)
  )
  
  out
}


# Set Working Directory 
source('code/helper_functions.R')
source('code/scale_function.R')

# Kernel specific parts
use_virtualenv("r-reticulate", required = TRUE)
#py_install("numpy<2.0")
#py_install("scikit-learn")
np <- import("numpy")
source_python("code/gp_simu_gate.py")

# Load Data
TX_data <- synth_tx_dataset(n = 1.1e6)

# Texas Data
dim(TX_data)

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

