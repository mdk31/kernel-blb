library(pbapply)
library(reticulate)
library(data.table)
library(ggplot2)


source('code/helper_functions.R')

kangschafer3_test <- function(n, te, sigma, beta_overlap = 0.5){
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
  
  Tr <- rbinom(n, 1, prt)
  
  Y0  <- X1 + X2 + rnorm(n, 0, sigma)
  Y1  <- Y0 + te + 0.2*te*X2
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X1, X2, Y1, Y0)
  out <- data.table::as.data.table(out)
  return(out)
}

te <- 0.8
sigma <- 1
replications <- 1000

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
n_values <- c(2500)

# SIMULATIONS----

out <- pblapply(seq_len(replications), function(rp){
  set.seed(rp)
  dat <- kangschafer3_test(n = n_values, te = te, sigma = sigma, beta_overlap = 0.5)
  satt <- mean(dat$Y1[dat$Tr == 1]) - mean(dat$Y0[dat$Tr == 1])
  
  att_estim <- att_kernel_weights(dat, 
                     Tr = 'Tr', 
                     Y = 'y', 
                     confounder_names = c('X1', 'X2'), 
                     degree1, 
                     degree2, 
                     k1, 
                     k2, 
                     operator, 
                     penal)
  
  return(data.table(true_satt = satt, estim_satt = att_estim))

}, cl = 1)

out <- rbindlist(out)

ggplot(out, aes(x = estim_satt)) +
  geom_histogram(position = 'identity', bins = 50) +
  geom_vline(xintercept = mean(out$true_satt), color = 'yellow') +
  theme_minimal() + 
  ggtitle('Kernel Weighted Estimator of SATT over 1000 replications (n = 5000)')
