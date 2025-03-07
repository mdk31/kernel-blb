
calculate_gamma <- function(n, subsets){
  soln <- 1 - log(subsets)/log(n)
  return(truncate_to_n(soln, 5))
}

causal_blb <- function(data, b, subsets, disjoint = TRUE, K = 10){
  assertthat::assert_that(is.integer(subsets))
  if(data.table::is.data.table(data) == FALSE){
    data <- data.table::as.data.table(data)
  }
  n <- nrow(data)
  partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = disjoint)
  idx <- seq_len(b)
  
  blb_out <- lapply(partitions, function(i){
    tmp_dat <- dat[i]
    folds <- split(idx, sample(rep(1:K, length.out = length(idx))))
    crossfit <- crossfit_estimator(data = tmp_dat, K = K)

    M <- rmultinom(n = B, size = n, prob = rep(1, b))
    blb_reps <- sapply(seq_len(B), function(bt){
      phi1 <- M[, bt]*((crossfit$Tr/crossfit$prop_score)*(crossfit$y - crossfit$m1) + crossfit$m1)
      phi0 <- M[, bt]*((1 - crossfit$Tr)/(1 - crossfit$prop_score)*(crossfit$y - crossfit$m0) + crossfit$m0)
      sum(phi1)/n - sum(phi0)/n
    })
    
    perc_ci <- boot:::perc.ci(blb_reps)
    return(data.table(lower_ci = perc_ci[4],
                      upper_ci = perc_ci[5],
                      estim = mean(blb_reps),
                      se = sd(blb_reps)))
  })
  
  blb_out <- rbindlist(blb_out)
  blb_out <- blb_out[, .(lower_ci = mean(lower_ci),
                         upper_ci = mean(upper_ci),
                         estim = mean(estim),
                         se = mean(se))]
  return(blb_out)
}

crossfit_estimator <- function(data, K = 10){
  if(data.table::is.data.table(data) == FALSE){
    data <- data.table::as.data.table(data)
  }
  
  fold_idx <- seq_len(nrow(data))
  folds <- split(fold_idx, sample(rep(1:K, length.out = length(fold_idx))))
  
  crossfit_dt <- lapply(folds, function(test_idx){
    train_idx <- setdiff(fold_idx, test_idx)
    train_dat <- data[train_idx]
    test_dat <- data[-train_idx]
    
    m <- e1071::svm(y ~ Tr + X1 + X2, data = train_dat, kernel = 'linear', type = 'eps-regression')
    g <- e1071::svm(Tr ~ X1 + X2, data = train_dat, kernel = 'linear', type = 'C-classification', probability = TRUE)
    
    prop_score <- attr(predict(g, newdata = test_dat, probability = TRUE), 'probabilities')[, 1]
    
    newdata <- data.frame(Tr = 1, X1 = test_dat$X1, X2 = test_dat$X2)
    m1 <- predict(m, newdata = newdata)
    
    newdata <- data.frame(Tr = 0, X1 = test_dat$X1, X2 = test_dat$X2)
    m0 <- predict(m, newdata = newdata)
    
    test_dat$prop_score <- prop_score
    test_dat$m1 <- m1
    test_dat$m0 <- m0
    
    test_dat
  })
  return(rbindlist(crossfit_dt))
}


kangschafer3 <- function(n, te, sigma, beta_overlap = 0.5){
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1)
  prt <- 1/(1 + exp(-beta_overlap*(X1 + X2)))
  
  Tr <- rbinom(n, 1, prt)
  
  Y0  <- X1 + X2 + rnorm(n, 0, sigma)
  Y1  <- Y0 + te
  y   <- Y0*(1-Tr) + Y1*(Tr)
  out <- cbind(y, Tr, X1, X2, Y1, Y0)
  out <- data.table::as.data.table(out)
  return(out)
}

make_partition <- function(n, subsets, b, disjoint = TRUE){
  part_idx <- seq(1, n, by = 1)
  if(disjoint){
    # Generate disjoint sets
    # Permute indices
    partition <- sample(part_idx, n)
    totality <- b*subsets
    stopifnot(totality <= n)
    # Exclude rest of sample
    partition <- partition[1:totality]
    partition <- split(partition, f = rep(1:subsets, each = b))
  } else{
    partition <- replicate(subsets, {
      sample(part_idx, size = b, replace = FALSE)
    }, simplify = FALSE)
  }
  partition
}

zip_plots <- function(data, n, use_case, subsets = NULL, gamma = NULL){
  assertthat::assert_that(!((is.null((subsets) & !is.null(gamma)) | (!is.null(subsets) & is.null(gamma)))))
  if(data.table::is.data.table(data) == FALSE){
    data <- data.table::as.data.table(data)
  }
  nm_prefix <- paste0(use_case, '_zip_plot_n', n)
  if(is.null(subsets)){
    nm <- paste0(nm_prefix, '_full.pdf')
    title <- bquote(n == .(n))
    ggsub <- data[n == n]
    label_sub <- zip_labels[n == n_val & subsets == subsets_val & gamma == gamma_val & B == 25]
    
  } else{
    nm <- paste0(nm_prefix, '_subset_', subsets, '_gamma_', gamma, '_cblb.pdf')
    title <- bquote(paste(s == .(subsets_val), ' and ', gamma == .(gamma_val), ' and ', n == .(n_val)))
    ggsub <- data[n == n & subsets == subsets & gamma == gamma]
    label_sub <- zip_labels[n == n_val & subsets == subsets_val & gamma == gamma_val & B == 25]
    
  }
  for(row_idx in seq_len(nrow(value_grid))){
    ggsub <- ggdat[n == n_val & subsets == subsets_val & gamma == gamma_val & B == 25]
    p <- ggplot(ggsub, aes(y = rank)) +
      geom_segment(aes(x = lower_ci, y = rank, xend = upper_ci, yend = rank, color = covered)) +
      facet_grid(`Propensity Model` ~ `Outcome Model`, labeller = label_both) +
      geom_vline(aes(xintercept = te), color = 'yellow', size = 0.5, linetype = 'dashed') +
      ylab('Fractional Centile of |z|') +
      xlab('95% Confidence Intervals') +
      theme_bw() +
      scale_y_continuous(breaks = c(5, 50, 95)) +
      scale_color_discrete(name = "Coverage") +
      geom_text(x = 0.75, y = 50, aes(label = perc_cover), data = label_sub, size = 4) +
      ggtitle(title)
    
    print(p)
    ggsave(file.path(image_path, nm), height = 9, width = 7)
  }
}

