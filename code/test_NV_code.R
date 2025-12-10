library(reticulate)
library(data.table)
library(ggplot2)
library(vroom)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

source('code/helper_functions.R')

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
    educ    = ifelse(meduc %in% c("", "9"), NA, meduc),
    race6   = (ifelse(mrace6 %in% c("", "9"), NA, mrace6)),
    hisp    = (ifelse(mhisp_r %in% c("", "9"), NA, mhisp_r)),
    married = (ifelse(dmar %in% c("", "9"), NA, dmar)),
    parity  = as.integer(ifelse(lbo_rec %in% c("", "9"), NA, lbo_rec)),
    precare = as.integer(ifelse(precare %in% c("", "99"), NA, precare)),
    pay     = (ifelse(pay %in% c("", "9"), NA, pay)),
    plural1 = as.integer(dplural) == 1,
    year    = as.integer(dob_yy)
  ) |>
  filter(plural1, dbwt >= 350, dbwt <= 6000) |>
  drop_na(dbwt, age, educ, race6, hisp, married, parity, precare, pay, any_smoke, year)

# Run-----
n <- nrow(nv)
subsets <- 2500
gamma <- calculate_gamma(n, subsets)
b <- floor(n^gamma)
B <- 100

# Scale all confounders to be on the same scale
confounders <- c('educ', 'race6', 'hisp', 'married', 'pay', 'age', 'parity', 'precare')


# time1 <- system.time({
#   out <- causal_blb_aipw(data = nv,
#                          weight_function = aipw_kernel_weights,
#                          y = 'dbwt',
#                          Tr = 'any_smoke',
#                          confounders = confounders,
#                          b = b,
#                          subsets = subsets,
#                          disjoint = TRUE,
#                          degree1, degree2, k1, k2, operator, penal)
# })

time2 <- system.time({
  out <- causal_blb_policy(data = nv,
                           y = 'dbwt',
                           A = 'any_smoke',
                           r_tilde_form = dbwt ~ educ + race6 + hisp + married + pay + age + parity + precare + any_smoke,
                           covariates = confounders,
                           initial_params =  c(rep(0, b), 0),
                           lambda = 0.01,
                           b = b, subsets = subsets)
})

# time2 <- system.time({
#   out2 <- causal_blb_aipw(data = nv,
#                          weight_function = parametric_weights,
#                          y = 'dbwt',
#                          Tr = 'any_smoke',
#                          confounders = confounders,
#                          b = b,
#                          subsets = subsets)
# })



# ggplot(out, aes(x = estim_satt)) +
#   geom_histogram(position = 'identity', bins = 50) +
#   geom_vline(xintercept = mean(out$true_satt), color = 'yellow') +
#   theme_minimal() + 
#   ggtitle('Kernel Weighted Estimator of SATT over 1000 replications (n = 5000)')
