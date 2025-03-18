library(pbapply)
library(data.table)

# TODO: install kernel packages
source('code/helper_functions.R')

te <- 0.8
sigma <- 1
replications <- 1000
B <- 100 #bootstrap resamples
optimal_val <- 1.001 #extracted from PAPER


base_nm <- 'optimal_policy'
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
n_values <- c(10000)
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
      dat <- aol_dgp(n = n)
      lambda <- 0.01
      initial_params <- c(rep(0, n), 0)  # Initial v and b
      estim_opt_regime <- estimate_optimal_regime(data, initial_params, lambda) 
      M <- rmultinom(n = B, size = n, prob = rep(1, n))
      
      boot_reps <- sapply(seq_len(B), function(bt){
        sum(M[, bt]*dat$y/0.5*(dat$A == estim_opt_regime))/n
      })
      
      perc_ci <- boot:::perc.ci(boot_reps)
      return(data.table(lower_ci = perc_ci[4],
                        upper_ci = perc_ci[5],
                        estim = mean(boot_reps),
                        se = sd(boot_reps)))
    }, cl = 4)
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
      return(causal_blb_policy(data = dat,
                               y = 'y',
                               A = 'A',
                               r_tilde_form = y ~ A,
                               covariates = c('x1', 'x2'),
                               initial_params =  c(rep(0, b), 0),
                               lambda = 0.01,
                               b = b, subsets = subsets))
    }, cl = 1)
    
    out <- rbindlist(out)
    out[, `:=`(n = n)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'cblb_bootstrap.rds'))
}
