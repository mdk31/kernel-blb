library(data.table)
library(ggplot2)
library(vroom)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(knitr)
library(kableExtra)
library(reticulate)


source('code/helper_functions.R')

base_nm <- 'case_study'
image_path <- 'images'
dat_path <- 'data'

# Create the temporary directory
temp_dir <- file.path(dat_path, paste0(base_nm, '_tmp'))
if(!file.exists(temp_dir)){
  dir.create(temp_dir, recursive = TRUE)
}

escape_latex <- function(x){
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

build_latex_table <- function(df, caption, label, align = NULL){
  if(is.null(align)){
    align <- paste(rep("l", ncol(df)), collapse = "")
  }

  df_chr <- as.data.frame(lapply(df, escape_latex), stringsAsFactors = FALSE)
  header <- paste(names(df_chr), collapse = " & ")
  rows <- apply(df_chr, 1, function(r) paste(r, collapse = " & "))

  lines <- c(
    "\\begin{table}[!h]",
    "\\centering",
    sprintf("\\caption{\\label{%s}%s}", label, escape_latex(caption)),
    "\\centering",
    sprintf("\\begin{tabular}[t]{%s}", align),
    "\\toprule",
    paste0(header, "\\\\"),
    "\\midrule",
    paste0(rows, "\\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )

  paste(lines, collapse = "\n")
}

# Code-----
PATH_DIR <- "~/Documents/Projects/kernel-blb/data"  # change if needed
DCT_URL  <- "https://data.nber.org/nvss/natality/programs/dct/natality2023.dct"

# locate the TXT
txt_candidates <- list.files(PATH_DIR, pattern="\\.txt$", full.names=TRUE, ignore.case=TRUE)
stopifnot(length(txt_candidates) >= 1)
TXT <- txt_candidates[1]

# load and parse the 2023 .dct
DCT_FILE <- file.path(tempdir(), "natality2023.dct")
download.file(DCT_URL, DCT_FILE, mode="wb", quiet=TRUE)
dct <- read_lines(DCT_FILE)

rx <- "_column\\((\\d+)\\)\\s+([a-z0-9]+)\\s+([A-Za-z0-9_]+)\\s+%([0-9]+)"
m  <- str_match(dct, rx)
dct_df <- data.frame(
  start = as.integer(m[,2]),
  type  = m[,3],
  name  = m[,4],
  width = as.integer(m[,5]),
  stringsAsFactors = FALSE
)
dct_df <- dct_df[!is.na(dct_df$start), ]
dct_df <- dct_df[order(dct_df$start), ]
dct_df$end <- dct_df$start + dct_df$width - 1
stopifnot(!is.na(TXT), file.exists(TXT))


# Byte positions (inclusive) from your 2023 .dct snippet
# dob_yy 9–12; mager 75–76; mrace6 107; mhisp_r 115; dmar 120; meduc 124;
# lbo_rec 179; precare 224–225; cig_1 255–256; cig_2 257–258; cig_3 259–260;
# dplural 454; dbwt 504–507; pay 435
pos <- fwf_positions(
  start = c(   9,  75, 107, 115, 120, 124, 179, 224, 255, 257, 259, 454, 504, 435),
  end   = c(  12,  76, 107, 115, 120, 124, 179, 225, 256, 258, 260, 454, 507, 435),
  col_names = c("dob_yy","mager","mrace6","mhisp_r","dmar","meduc","lbo_rec",
                "precare","cig_1","cig_2","cig_3","dplural","dbwt","pay")
)

# Read only those columns (fast, robust to .dct quirks)
raw <- vroom_fwf(TXT, col_positions = pos, progress = TRUE, altrep = TRUE)

# Minimal hygiene and variable construction
nv <- raw |>
  mutate(across(everything(), \(x) na_if(trimws(as.character(x)), ""))) |>
  transmute(
    # outcomes
    dbwt    = as.integer(ifelse(dbwt == "9999", NA, dbwt)),
    lbw     = as.integer(dbwt < 2500),

    # exposure
    cig1    = as.integer(ifelse(cig_1 %in% c("", "99"), NA, cig_1)),
    cig2    = as.integer(ifelse(cig_2 %in% c("", "99"), NA, cig_2)),
    cig3    = as.integer(ifelse(cig_3 %in% c("", "99"), NA, cig_3)),
    any_smoke = as.integer((coalesce(cig1, 0) > 0) | (coalesce(cig2, 0) > 0) | (coalesce(cig3, 0) > 0)),

    # confounders
    age     = as.integer(ifelse(mager %in% c("", "99"), NA, mager)),
    educ    = as.integer(ifelse(meduc %in% c("", "9"), NA, meduc)),
    race6   = as.integer(ifelse(mrace6 %in% c("", "9"), NA, mrace6)),
    hisp    = as.integer(ifelse(mhisp_r %in% c("", "9"), NA, mhisp_r)),
    married = as.integer(ifelse(dmar %in% c("", "9"), NA, dmar)),
    parity  = as.integer(ifelse(lbo_rec %in% c("", "9"), NA, lbo_rec)),
    precare = as.integer(ifelse(precare %in% c("", "99"), NA, precare)),
    pay     = as.integer(ifelse(pay %in% c("", "9"), NA, pay)),
    plural1 = as.integer(dplural) == 1,
    year    = as.integer(dob_yy)
  ) |>
  filter(plural1, dbwt >= 350, dbwt <= 6000) |>
  drop_na(dbwt, age, educ, race6, hisp, married, parity, precare, pay, any_smoke, year) %>%
  select(-c(cig1, cig2, cig3))

# Baseline characteristics by smoking status -----
nv_baseline <- nv %>%
  mutate(
    smoking_status = if_else(any_smoke == 1, "Smoker", "Non-smoker"),
    age_group = cut(
      age,
      breaks = c(-Inf, 19, 24, 29, 34, 39, Inf),
      labels = c("<20", "20-24", "25-29", "30-34", "35-39", "40+"),
      right = TRUE
    )
  ) %>%
  mutate(
    # Mother’s race recode (6 categories)
    race6 = factor(
      race6,
      levels = c(1, 2, 3, 4, 5, 6),
      labels = c(
        "White (only)",
        "Black (only)",
        "AIAN (only)",
        "Asian (only)",
        "NHOPI (only)",
        "More than one race"
      )
    ),
    
    # Mother’s Hispanic origin recode
    hisp = factor(
      hisp,
      levels = c(0, 1, 2, 3, 4, 5),
      labels = c(
        "Non-Hispanic",
        "Mexican",
        "Puerto Rican",
        "Cuban",
        "Central and South American",
        "Other and Unknown Hispanic origin"
      )
    ),
    
    # Mother’s education
    educ = factor(
      educ,
      levels = c(1, 2, 3, 4, 5, 6, 7, 8),
      labels = c(
        "8th grade or less",
        "9th through 12th grade with no diploma",
        "High school graduate or GED completed",
        "Some college credit, but not a degree",
        "Associate degree (AA, AS)",
        "Bachelor’s degree (BA, AB, BS)",
        "Master’s degree (MA, MS, MEng, MEd, MSW, MBA)",
        "Doctorate (PhD, EdD) or Professional Degree (MD, DDS, DVM, LLB, JD)"
      )
    ),
    
    # Live birth order recode
    parity = factor(
      parity,
      levels = c(1, 2, 3, 4, 5, 6, 7, 8),
      labels = c(
        "1st live birth",
        "2nd live birth",
        "3rd live birth",
        "4th live birth",
        "5th live birth",
        "6th live birth",
        "7th live birth",
        "8 or more live births"
      )
    ),
    
    # Marital status (covers both the U.S. and Puerto Rico encodings that can appear in the file)
    married = factor(
      married,
      levels = c(1, 2, 3),
      labels = c("Married", "Unmarried", "Unmarried")
    ),
    
    # Month prenatal care began
    precare = case_when(
      precare == 0 ~ "No prenatal care",
      precare >= 1 & precare <= 10 ~ paste0("Month prenatal care began: ", precare),
      TRUE ~ as.character(precare)
    ),
    
    # Payment source for delivery
    pay = factor(
      pay,
      levels = c(1, 2, 3, 4, 5, 6, 8),
      labels = c(
        "Medicaid",
        "Private Insurance",
        "Self-Pay",
        "Indian Health Service",
        "CHAMPUS/TRICARE",
        "Other Government (Federal, State, Local)",
        "Other"
      )
    )
  )

# Sample size + birthweight summary by smoking status
baseline_by_smoking <- nv_baseline %>%
  group_by(smoking_status) %>%
  summarise(
    sample_size = n(),
    pct_of_total = 100 * n() / nrow(nv_baseline),
    mean_birthweight = mean(dbwt),
    sd_birthweight = sd(dbwt),
    mean_age = mean(age),
    sd_age = sd(age),
    .groups = "drop"
  )

# Covariate distributions by smoking status (n and % within smoking group)
covariates_for_distribution <- c(
  "age_group", "race6", "educ", "parity", "hisp", "married", "precare", "pay"
)

baseline_covariate_distribution <- nv_baseline %>%
  select(smoking_status, all_of(covariates_for_distribution)) %>%
  mutate(across(all_of(covariates_for_distribution), as.character)) %>%
  pivot_longer(
    cols = -smoking_status,
    names_to = "covariate",
    values_to = "level"
  ) %>%
  group_by(smoking_status, covariate, level) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

# Wide-format table for reporting
baseline_covariate_distribution_wide <- baseline_covariate_distribution %>%
  mutate(value = sprintf("%s (%.1f%%)", format(n, big.mark = ","), pct)) %>%
  select(covariate, level, smoking_status, value) %>%
  pivot_wider(names_from = smoking_status, values_from = value)

baseline_by_smoking_tbl <- baseline_by_smoking %>%
  transmute(
    `Smoking status` = smoking_status,
    `N` = format(sample_size, big.mark = ","),
    `% of sample` = sprintf("%.1f", pct_of_total),
    `Mean birthweight` = sprintf("%.3f", mean_birthweight),
    `SD birthweight` = sprintf("%.3f", sd_birthweight),
    `Mean age` = sprintf("%.3f", mean_age),
    `SD age` = sprintf("%.3f", sd_age)
  )

baseline_covariate_tbl <- baseline_covariate_distribution_wide %>%
  transmute(
    Covariate = covariate,
    Level = level,
    `Non-smoker` = `Non-smoker`,
    Smoker = Smoker
  )

baseline_smoking_tex <- build_latex_table(
  df = baseline_by_smoking_tbl,
  caption = "Baseline characteristics by smoking status (NVSS)",
  label = "tab:baseline_smoking",
  align = "lllllll"
)

baseline_covariate_tex <- build_latex_table(
  df = baseline_covariate_tbl,
  caption = "Covariate distributions by smoking status (NVSS)",
  label = "tab:baseline_covariates",
  align = "llll"
)

writeLines(
  baseline_smoking_tex,
  file.path(temp_dir, "baseline_by_smoking_table.tex")
)
writeLines(
  baseline_covariate_tex,
  file.path(temp_dir, "baseline_covariate_distribution_table.tex")
)

# Run-----
n <- nrow(nv)
subsets <- 2500
gamma <- calculate_gamma(n, subsets)
b <- floor(n^gamma)
B <- 100

# Scale all confounders to be on the same scale
confounders <- c('educ', 'race6', 'hisp', 'married', 'pay', 'age', 'parity', 'precare')

if(file.exists(file.path(temp_dir, 'policy_bootstrap.rds'))){
  out <- readRDS(file.path(temp_dir, 'policy_bootstrap.rds'))
} else{
  time2 <- system.time({
    out <- causal_blb_policy(data = nv,
                             y = 'dbwt',
                             A = 'any_smoke',
                             r_tilde_form = dbwt ~ educ + race6 + hisp + married + pay + age + parity + precare + any_smoke,
                             covariates = confounders,
                             initial_params =  c(rep(0, b), 0),
                             lambda = 0.01,
                             b = b, subsets = subsets,
                             balance_treated = TRUE, randomized = FALSE)
  })
  out[, elapsed := unname(time2[["elapsed"]])]
  saveRDS(out, file.path(temp_dir, 'policy_bootstrap.rds'))
}

if(file.exists(file.path(temp_dir, 'aipw_bootstrap.rds'))){
  out <- readRDS(file.path(temp_dir, 'aipw_bootstrap.rds'))
  
} else{
  # Kernel specific parts
  use_virtualenv("r-reticulate", required = TRUE)
  np <- import("numpy")
  source_python("code/gp_simu_gate.py")
  degree1 <- 1
  degree2 <- 1
  k1 <- "poly"
  k2 <- "poly"
  operator <- "single"
  penal <- log(2)
  time2 <- system.time({
    out <- causal_blb_aipw(data = nv, 
                           y = 'dbwt', 
                           Tr = 'any_smoke', confounders = confounders,
                           b = b, subsets = subsets, disjoint = TRUE, weight_function = aipw_kernel_weights,
                           balance_treated = TRUE, return_ey0 = FALSE,
                           degree1, degree2, k1, k2, operator, penal)    

  })
  out[, elapsed := unname(time2[["elapsed"]])]
  saveRDS(out, file.path(temp_dir, 'aipw_bootstrap.rds'))
}

if(file.exists(file.path(temp_dir, 'aipw_bootstrap_ey0.rds'))){
  out <- readRDS(file.path(temp_dir, 'aipw_bootstrap_ey0.rds'))
  
} else{
  # Kernel specific parts
  use_virtualenv("r-reticulate", required = TRUE)
  np <- import("numpy")
  source_python("code/gp_simu_gate.py")
  degree1 <- 1
  degree2 <- 1
  k1 <- "poly"
  k2 <- "poly"
  operator <- "single"
  penal <- log(2)
  out <- causal_blb_aipw(data = nv, 
                         y = 'dbwt', 
                         Tr = 'any_smoke', confounders = confounders,
                         b = b, subsets = subsets, disjoint = TRUE, weight_function = aipw_kernel_weights,
                         balance_treated = TRUE, return_ey0 = TRUE,
                         degree1, degree2, k1, k2, operator, penal)    
  saveRDS(out, file.path(temp_dir, 'aipw_bootstrap_ey0.rds'))
}

# # FULL BOOTSTRAP
# if(file.exists(file.path(temp_dir, 'aipw_bootstrap_full.rds'))){
#   out <- readRDS(file.path(temp_dir, 'aipw_bootstrap_full.rds'))
# } else{
#   subsets_full <- 1
#   gamma_full <- calculate_gamma(n, subsets_full)
#   b_full <- floor(n^gamma_full)
#   # Kernel specific parts
#   use_virtualenv("r-reticulate", required = TRUE)
#   np <- import("numpy")
#   source_python("code/gp_simu_gate.py")
#   degree1 <- 1
#   degree2 <- 1
#   k1 <- "poly"
#   k2 <- "poly"
#   operator <- "single"
#   penal <- log(2)
#   time2 <- system.time({
#     out <- causal_blb_aipw(data = nv, 
#                            y = 'dbwt', 
#                            Tr = 'any_smoke', confounders = confounders,
#                            b = b_full, subsets = subsets_full, disjoint = FALSE, weight_function = aipw_kernel_weights,
#                            balance_treated = FALSE, degree1, degree2, k1, k2, operator, penal)    
#     
#   })
#   out[, elapsed := unname(time2[["elapsed"]])]
#   saveRDS(out, file.path(temp_dir, 'aipw_bootstrap_full.rds'))
# }

if(file.exists(file.path(temp_dir, 'dml_bootstrap.rds'))){
  out <- readRDS(file.path(temp_dir, 'dml_bootstrap.rds'))
} else{
  subsets <- 300
  gamma <- calculate_gamma(n, subsets)
  b <- floor(n^gamma)
  time2 <- system.time({
    out <- causal_blb(data = nv, 
                      y = 'dbwt', 
                      Tr = 'any_smoke', confounders = confounders,,
                      b = b, subsets = subsets, balance_treated = TRUE,
                      return_ey0 = FALSE)
  })
  out[, elapsed := unname(time2[["elapsed"]])]
  saveRDS(out, file.path(temp_dir, 'dml_bootstrap.rds'))
}

if(file.exists(file.path(temp_dir, 'dml_bootstrap_ey0.rds'))){
  out <- readRDS(file.path(temp_dir, 'dml_bootstrap_ey0.rds'))
} else{
  subsets <- 300
  gamma <- calculate_gamma(n, subsets)
  b <- floor(n^gamma)
  time2 <- system.time({
    out <- causal_blb(data = nv, 
                      y = 'dbwt', 
                      Tr = 'any_smoke', confounders = confounders,,
                      b = b, subsets = subsets, balance_treated = TRUE,
                      return_ey0 = TRUE)
  })
  out[, elapsed := unname(time2[["elapsed"]])]
  saveRDS(out, file.path(temp_dir, 'dml_bootstrap_ey0.rds'))
}


# if(file.exists(file.path(temp_dir, 'dml_bootstrap_full.rds'))){
#   out <- readRDS(file.path(temp_dir, 'dml_bootstrap_full.rds'))
# } else{
#   subsets_full <- 1
#   gamma_full <- calculate_gamma(n, subsets_full)
#   b_full <- floor(n^gamma_full)
#   time2 <- system.time({
#     out <- causal_blb(data = nv,
#                       y = 'dbwt',
#                       Tr = 'any_smoke', confounders = confounders,
#                       b = b_full, subsets = subsets_full, balance_treated = FALSE,
#                       disjoint = FALSE)
#   })
#   out[, elapsed := unname(time2[["elapsed"]])]
#   saveRDS(out, file.path(temp_dir, 'dml_bootstrap_full.rds'))
# }

# NO BLB
# if(file.exists(file.path(temp_dir, 'full_bootstrap.rds'))){
#   out <- readRDS(file.path(temp_dir, 'full_bootstrap.rds'))
# } else{
#   nv2 <- copy(nv) |>
#     slice_sample(prop = 0.01)
#   time2 <- system.time({
#     y <- 'dbwt'
#     A <- 'any_smoke'
#     n <- nrow(nv2)
#     ps_dat <- copy(nv2)
#     lambda <- 0.01
#     ps_dat$A01 <- ifelse(ps_dat[[A]] == 1, 1, 0)
#     
#     ps_formula <- as.formula(
#       paste("A01 ~", paste(confounders, collapse = " + "))
#     )
#     ps_fit <- glm(ps_formula, data = ps_dat, family = binomial)
#     pi_1 <- predict(ps_fit, newdata = ps_dat, type = "response")
#     pi_1 <- pmin(pmax(pi_1, 1e-6), 1 - 1e-6)
#     pi_A <- ifelse(nv[[A]] == 1, pi_1, 1 - pi_1)
#     
#     
#     estim_opt_regime <- estimate_optimal_regime(data = nv2,
#                                                 r_tilde_form = dbwt ~ educ + race6 + hisp + married + pay + age + parity + precare + any_smoke,
#                                                 covariates = confounders,
#                                                 A = A,
#                                                 y = y,
#                                                 initial_params = c(rep(0, b), 0), 
#                                                 lambda = lambda,
#                                                 boundary = 0,
#                                                 pi_1 = pi_1)
#     
#     M <- rmultinom(n = B, size = n, prob = rep(1, n))
#     
#     blb_reps <- sapply(seq_len(B), function(bt){
#       sum(M[, bt]*nv2[[y]]/pi_A*(nv2[[A]] == estim_opt_regime))/n
#     })
#     
#     perc_ci <- boot:::perc.ci(blb_reps)
#   })
#   out[, elapsed := unname(time2[["elapsed"]])]
#   saveRDS(out, file.path(temp_dir, 'cblb_bootstrap.rds'))
# }

policy_out <- readRDS(file.path(temp_dir, "policy_bootstrap.rds"))
aipw_out   <- readRDS(file.path(temp_dir, "aipw_bootstrap.rds"))
dml_out    <- readRDS(file.path(temp_dir, "dml_bootstrap.rds"))

# Read ey0 from dedicated runs when available and append to existing data
aipw_ey0_file <- file.path(temp_dir, "aipw_bootstrap_ey0.rds")
dml_ey0_file  <- file.path(temp_dir, "dml_bootstrap_ey0.rds")
if (file.exists(aipw_ey0_file)) {
  aipw_ey0 <- readRDS(aipw_ey0_file)
  if ("ey0" %in% names(aipw_ey0)) aipw_out[, ey0 := aipw_ey0$ey0]
}
if (file.exists(dml_ey0_file)) {
  dml_ey0 <- readRDS(dml_ey0_file)
  if ("ey0" %in% names(dml_ey0)) dml_out[, ey0 := dml_ey0$ey0]
}

if (!("ey0" %in% names(policy_out))) policy_out[, ey0 := NA_real_]
if (!("ey0" %in% names(aipw_out))) aipw_out[, ey0 := NA_real_]
if (!("ey0" %in% names(dml_out))) dml_out[, ey0 := NA_real_]

policy_out[, output_type := "Kernelized AOL"]
aipw_out[,   output_type := "Kernel minimax weights"]
dml_out[,    output_type := "DML with SVM"]

# Combine into one table (3 rows)
tab <- data.table::rbindlist(
  list(aipw_out, dml_out, policy_out),
  use.names = TRUE,
  fill = TRUE
)

results_tbl <- as.data.frame(tab)
col_map <- c(
  output_type = "Model",
  ey0 = "E(Y(0))",
  lower_ci = "Lower CI",
  upper_ci = "Upper CI",
  estim = "Estim.",
  se = "S.E.",
  elapsed = "Time (seconds)"
)

present_cols <- names(col_map)[names(col_map) %in% names(results_tbl)]
results_tbl <- results_tbl[, present_cols, drop = FALSE]
names(results_tbl) <- unname(col_map[present_cols])

num_cols <- setdiff(names(results_tbl), "Model")
results_tbl[num_cols] <- lapply(
  results_tbl[num_cols],
  function(x) sprintf("%.3f", as.numeric(x))
)

results_tex <- build_latex_table(
  df = results_tbl,
  caption = "Treatment Effects and Optimal Outcome for NVSS dataset",
  label = "tab:casestudy",
  align = paste(rep("l", ncol(results_tbl)), collapse = "")
)

writeLines(results_tex, file.path(temp_dir, "cblb_bootstrap_table.tex"))
