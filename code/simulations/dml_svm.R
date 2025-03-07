

source('code/helper_funcs.R')


te <- 0.8
sigma <- 1
replications <- 1000
K <- 10
r <- 100

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
B_vals <- 100 #bootstrap resamples


# FULL SIMULATIONS----
if(file.exists(file.path(temp_dir, 'full_bootstrap.rds'))){
  cblb <- readRDS(file.path(temp_dir, 'full_bootstrap.rds'))
} else{
  cblb <- lapply(seq_row, function(i){
    grid_val <- hyper_grid[i]
    n <- grid_val$n
    B <- grid_val$B
    subsets <- grid_val$subsets
    
    out <- pblapply(seq_len(replications), function(rp){
      set.seed(rp)
      dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)
      idx <- seq_len(nrow(dat))
      folds <- split(idx, sample(rep(1:K, length.out = length(idx))))
      crossfit <- lapply(folds, function(test_idx){
        train_idx <- setdiff(idx, test_idx)
        train_dat <- dat[train_idx]
        test_dat <- dat[-train_idx]
        
        m <- svm(y ~ Tr + X1 + X2, data = train_dat, kernel = 'linear', type = 'eps-regression')
        g <- svm(Tr ~ X1 + X2, data = train_dat, kernel = 'linear', type = 'C-classification', probability = TRUE)
        
        prop_score <- attr(predict(g, newdata = test_dat, probability = TRUE), 'probabilities')[, 1]
        
        newdata <- data.frame(Tr = 1, X1 = test_dat$X1, X2 = test_dat$X2, Z1 = test_dat$Z1, Z2 = test_dat$Z2)
        m1 <- predict(m, newdata = newdata)
        
        newdata <- data.frame(Tr = 0, X1 = test_dat$X1, X2 = test_dat$X2, Z1 = test_dat$Z1, Z2 = test_dat$Z2)
        m0 <- predict(m, newdata = newdata)
        
        test_dat$prop_score <- prop_score
        test_dat$m1 <- m1
        test_dat$m0 <- m0
        
        test_dat
      })
      crossfit <- rbindlist(crossfit)
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
    }, cl = 4)
    
    out <- rbindlist(out)
    out[, `:=`(n = n,
               B = B)]
    out
  })
  cblb <- rbindlist(cblb)
  saveRDS(cblb, file.path(temp_dir, 'full_bootstrap.rds'))
}
