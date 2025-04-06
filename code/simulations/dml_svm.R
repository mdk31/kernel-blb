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
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      crossfit <- crossfit_estimator(dat, K = K)
      M <- rmultinom(n = B, size = n, prob = rep(1, n))
      
      boot_reps <- sapply(seq_len(B), function(bt){
        phi1 <- M[, bt]*((crossfit$Tr/crossfit$prop_score)*(crossfit$y - crossfit$m1) + crossfit$m1)
        phi0 <- M[, bt]*((1 - crossfit$Tr)/(1 - crossfit$prop_score)*(crossfit$y - crossfit$m0) + crossfit$m0)
        sum(phi1)/n - sum(phi0)/n
      })
      
      boot_ci <- boot:::perc.ci((boot_reps))
      return(data.table(lower_ci = boot_ci[4],
                            upper_ci = boot_ci[5],
                            estim = mean(boot_reps),
                            se = mean(boot_reps)))
    }, cl = 4)
    out <- rbindlist(out)
    out[, `:=`(n = n, subsets = subsets, gamma = gamma)]
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
    b <- floor(n^gamma)
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      return(causal_blb(data = dat, b = b, subsets = subsets))
      
    }, cl = 4)
    
    out <- rbindlist(out)
    out[, `:=`(n = n)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'cblb_bootstrap.rds'))
}


zip_plot_obj <- zip_plots_helper(data = cblb, type = 'blb', te = te)
for(idx in seq_len(nrow(grid_vals))){
  grid_row <- grid_vals[idx, ]
  n <- grid_row$n
  subsets <- grid_row$subsets
  gamma <- grid_row$gamma
  zip_plots(data = zip_plot_obj$zip, zip_labels = zip_plot_obj$zip_labels, n = n, 
            type = 'cblb', use_case = 'dml', plot_title = 'DML SVM', 
            te = te, image_path = img_tmp_dir, text_x = 0.89)
}
