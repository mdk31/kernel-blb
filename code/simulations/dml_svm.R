library(pbapply)
library(data.table)

source('code/helper_functions.R')


te <- 0.8
sigma <- 1
replications <- 1000
K <- 10
B <- 100 #bootstrap resamples

base_nm <- 'svm_dml'
image_path <- 'images'
dat_path <- 'data'

# Create the temporary directory
temp_dir <- file.path(dat_path, paste0(base_nm, '_tmp'))
if(!file.exists(temp_dir)){
  dir.create(temp_dir, recursive = TRUE)
}

img_tmp_dir <- file.path(image_path, paste0(base_nm, '_tmp'))
if(!file.exists(img_tmp_dir)){
  dir.create(img_tmp_dir, recursive = TRUE)
}

# Values for simulations
n_values <- c(10000, 50000)
subset_values <- c(5, 10, 15)

# FULL SIMULATIONS----
grid_vals <- as.data.table(expand.grid(n = n_values))
seq_row <- seq_len(nrow(grid_vals))

if(file.exists(file.path(temp_dir, 'full_bootstrap.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'full_bootstrap.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- grid_vals[i]
    n <- grid_val$n

    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      crossfit <- crossfit_estimator(dat, K = K)
      M <- rmultinom(n = B, size = n, prob = rep(1, n))
      
      boot_reps <- sapply(seq_len(B), function(bt){
        phi1 <- M[, bt]*((crossfit$Tr/crossfit$prop_score)*(crossfit$y - crossfit$m1) + crossfit$m1)
        phi0 <- M[, bt]*((1 - crossfit$Tr)/(1 - crossfit$prop_score)*(crossfit$y - crossfit$m0) + crossfit$m0)
        sum(phi1)/n - sum(phi0)/n
      })
      
      boot_ci <- boot:::perc.ci((boot_reps))
      blb_out <- data.table(lower_ci = boot_ci[4],
                            upper_ci = boot_ci[5],
                            estim = mean(boot_reps),
                            se = mean(boot_reps))
      blb_out
    }, cl = 1)
    out <- rbindlist(out)
    out[, `:=`(n = n)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'full_bootstrap.rds'))
}


# cBLB SIMULATIONS----
grid_vals <- as.data.table(expand.grid(n = n_values,
                         subsets = subset_values))
grid_vals[, `:=`(gamma = calculate_gamma(n, subsets))]
seq_row <- seq_len(nrow(grid_vals))

if(file.exists(file.path(temp_dir, 'cblb_bootstrap.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'cblb_bootstrap.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- grid_vals[i]
    n <- grid_val$n
    subsets <- grid_val$subsets
    gamma <- grid_val$gamma
    b <- round(n^gamma)
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      crossfit <- crossfit_estimator(dat)
      M <- rmultinom(n = B, size = n, prob = rep(1, n))
      return(causal_blb(data = dat, b = b, subsets = subsets))
      
    }, cl = 1)
    
    out <- rbindlist(out)
    out[, `:=`(n = n)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'cblb_bootstrap.rds'))
}
