library(pbapply)
library(reticulate)
library(data.table)

source('code/helper_functions.R')

te <- 0.8
sigma <- 1
replications <- 1000
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

base_nm <- 'aipw_kernel_balancing'
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
  saveRDS(cblb, file.path(temp_dir, 'cblb_bootstrap.rds'))
}

zip_plot_obj <- zip_plots_helper(data = cblb, type = 'blb', te = te)
for(idx in seq_len(nrow(grid_vals))){
  grid_row <- grid_vals[idx, ]
  n <- grid_row$n
  subsets <- grid_row$subsets
  gamma <- grid_row$gamma
  zip_plots(data = zip_plot_obj$zip, zip_labels = zip_plot_obj$zip_labels, n = n, 
            type = 'cblb', use_case = 'aipw', plot_title = 'AIPW', 
            te = te, image_path = img_tmp_dir, text_x = 0.89)
}
