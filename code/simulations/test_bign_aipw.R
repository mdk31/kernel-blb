library(pbapply)
library(reticulate)
library(data.table)

source('code/helper_functions.R')

kangschafer3_test <- function(n, te, sigma, beta_overlap = 0.5) {
  ## ── Group 1: 3 independent standard normals ───────────────────────────────
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  
  ## ── Group 2: 2 independent uniforms on (-1,1) ─────────────────────────────
  X4 <- runif(n, -1, 1)
  X5 <- runif(n, -1, 1)
  
  ## ── Group 3: 3-variate normal, ρ = 0.6 ────────────────────────────────────
  Sigma3 <- matrix(0.6, 3, 3); diag(Sigma3) <- 1
  X678  <- MASS::mvrnorm(n, mu = rep(0, 3), Sigma = Sigma3)
  X6 <- X678[, 1];  X7 <- X678[, 2];  X8 <- X678[, 3]
  
  ## ── Group 4: 2 correlated Bernoulli via latent normal, ρ = 0.5 ───────────
  Sigma2 <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  latent <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma2)
  X9  <- as.numeric(latent[, 1] > 0)
  X10 <- as.numeric(latent[, 2] > 0)
  
  ## ── Propensity score & treatment assignment ───────────────────────────────
  prt <- 1 / (1 + exp(-beta_overlap * (X1 + X2 + X4)))  # mix of groups
  Tr  <- rbinom(n, 1, prt)
  
  ## ── Potential outcomes & observed outcome ─────────────────────────────────
  Y0 <- X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + rnorm(n, 0, sigma)
  Y1 <- Y0 + te
  y  <- ifelse(Tr == 1, Y1, Y0)
  
  ## ── Assemble and return data.table ────────────────────────────────────────
  out <- data.table::as.data.table(
    cbind(y, Tr,
          X1, X2, X3, X4, X5, X6, X7, X8, X9, X10,
          Y1, Y0)
  )
  return(out)
}
te <- 0.8
sigma <- 1
replications <- 5
B <- 100 #bootstrap resamples

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

# Values for simulations
n_values <- c(10*1e6)
subset_values <- c(20000)

# cBLB SIMULATIONS----
grid_vals <- as.data.table(expand.grid(n = n_values,
                         subsets = subset_values))
grid_vals[, `:=`(gamma = calculate_gamma(n, subsets))]
seq_row <- seq_len(nrow(grid_vals))

cblb <- lapply(seq_row, function(i){
  grid_val <- grid_vals[i]
  n <- grid_val$n
  subsets <- grid_val$subsets
  gamma <- grid_val$gamma
  b <- floor(n^gamma)
  
  out <- pblapply(seq_len(replications), function(rp){
    set.seed(rp)
    dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
    return(causal_blb_aipw(data = dat, y = 'y', Tr = 'Tr', confounders = c('X1', 'X2'),
                           b = b, subsets = subsets, degree1, degree2, k1, k2, operator, penal))
  }, cl = 1)
  
  out <- rbindlist(out)
  out[, `:=`(n = n, subsets = subsets, gamma = gamma)]
  out
})
cblb <- rbindlist(cblb)

