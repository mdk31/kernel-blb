library(pbapply)
library(data.table)

source('code/helper_functions.R')
Rcpp::sourceCpp("code/RBF_kernel_C_parallel.cpp")

te <- 0.8
sigma <- 1
replications <- 25
B <- 100 # bootstrap resamples

base_nm <- 'stable_balancing_noaugmentation_timing'
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
n_values <- c(100000)
subset_values <- c(5, 10)

# FULL SIMULATIONS----
grid_vals <- as.data.table(expand.grid(n = n_values, 
                                       kernel_approx = c(TRUE),
                                       kernel_type = c('rbfdot', 'vanilladot')))
grid_vals <- grid_vals[!(kernel_type == 'vanilladot' & kernel_approx)]
seq_row <- seq_len(nrow(grid_vals))

if(file.exists(file.path(temp_dir, 'full_bootstrap.rds'))){
  full <- readRDS(file.path(temp_dir, 'full_bootstrap.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- grid_vals[i]
    n <- grid_val$n
    kernel_approx <- grid_val$kernel_approx
    kernel_type <- grid_val$kernel_type
    if(kernel_type == 'vanilladot'){
      eig_clip <- 1e-10
    } else{
      eig_clip <- NULL
    }

    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      time <- system.time({
        output <- osqp_kernel_sbw_twofit(X = as.matrix(dat[, c('X1', 'X2')]),
                                         A = dat$Tr,
                                         Y = dat$y,
                                         delta.v=1e-4,
                                         kernel.approximation = kernel_approx,
                                         kernel_type = kernel_type, eig_clip = eig_clip)
        
        weights0 <- output[[1]]$res0$x
        weights1 <- output[[1]]$res1$x
        dat$full_weights <- NA
        dat$full_weights[dat$Tr == 1] <- weights1
        dat$full_weights[dat$Tr == 0] <- weights0
        dat$full_weights <- dat$full_weights*n
        # dat$full_weights <- output[[1]]$res$x*n
        
        phi1 <- (dat$Tr)*dat$full_weights*(dat$y)
        phi0 <- (1-dat$Tr)*dat$full_weights*(dat$y)
        M <- rmultinom(n = B, size = n, prob = rep(1, n))
        
        blb_reps <- sapply(seq_len(B), function(bt){
          boot_phi1 <- M[, bt]*phi1
          boot_phi0 <- M[, bt]*phi0
          mean(boot_phi1) - mean(boot_phi0)
        })
        perc_ci <- boot:::perc.ci(blb_reps)
      })
      
      return(data.table(time = time['elapsed']))
    }, cl = 1)
    out <- rbindlist(out)
    out[, `:=`(n = n, kernel_approx = kernel_approx, kernel_type = kernel_type)]
    out
  })
  full <- rbindlist(cblb)
  saveRDS(full, file.path(temp_dir, 'full_bootstrap.rds'))
}


# cBLB SIMULATIONS----
grid_vals <- as.data.table(expand.grid(n = n_values,
                         subsets = subset_values,
                         kernel_approx = c(TRUE),
                         kernel_type = c('rbfdot', 'vanilladot')))
grid_vals <- grid_vals[!(kernel_type == 'vanilladot' & kernel_approx)]
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
    kernel_approx <- grid_val$kernel_approx
    kernel_type <- grid_val$kernel_type
    if(kernel_type == 'vanilladot'){
      eig_clip <- 1e-10
    } else{
      eig_clip <- NULL
    }
    
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      time <- system.time({
        causal_blb_stable(data = dat, b = b, subsets = subsets, kernel_approx = kernel_approx,
                          augment = FALSE, eig_clip = eig_clip, kernel_type = kernel_type)
      })
      return(data.table(time = time['elapsed']))
    }, cl = 1)
    
    out <- rbindlist(out)
    out[, `:=`(n = n, subsets = subsets, gamma = gamma, kernel_approx = kernel_approx,
               kernel_type = kernel_type)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'cblb_bootstrap.rds'))
}

summary(full$time)
summary(cblb$time)
