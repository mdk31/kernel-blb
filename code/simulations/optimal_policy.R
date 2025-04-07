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
n_values <- c(5000)
subset_values <- c(5, 10, 15)

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
      dat <- aol_dgp(n = n)
      return(causal_blb_policy(data = dat,
                               y = 'y',
                               A = 'A',
                               r_tilde_form = y ~ x1 + x2 + x3 + x4 + x5 + A + A:x1 + A:x2,
                               covariates = c('x1', 'x2', 'x3', 'x4', 'x5'),
                               initial_params =  c(rep(0, b), 0),
                               lambda = 0.01,
                               b = b, subsets = subsets))
    }, cl = 4)
    
    out <- rbindlist(out)
    out[, `:=`(n = n, subsets = subsets, gamma = gamma)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'cblb_bootstrap.rds'))
}

zip_plot_obj <- zip_plots_helper(data = cblb, type = 'blb', te = optimal_val)
for(idx in seq_len(nrow(grid_vals))){
  grid_row <- grid_vals[idx, ]
  n <- grid_row$n
  subsets <- grid_row$subsets
  gamma <- grid_row$gamma
  zip_plots(data = zip_plot_obj$zip, zip_labels = zip_plot_obj$zip_labels, n = n, 
            type = 'cblb', use_case = 'policy', plot_title = 'Optimal value', te = optimal_val, image_path = img_tmp_dir, text_x = 0.89)
}


