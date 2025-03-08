aipw_kernel_weights <- function(data, degree1, degree2, k1, k2, operator, penal, bootstrap_size=length(data)){
  
  intervention <- data$Tr
  outcome <- data$y
  confounders <- data.frame(data$X1,data$X2)
  n <- nrow(data)
  
  t1 <- as.integer(intervention)
  t0 <- as.integer((1-intervention))
  y <- outcome
  X <- data.frame(confounders, y, intervention)
  
  X1t <- X[X$intervention == 1, ]
  X0t <- X[X$intervention == 0, ]
  y1 <- y[X$intervention == 1]
  y0 <- y[X$intervention == 0]
  n1 <- sum(X$intervention == 1)
  n0 <- sum(X$intervention == 0)
  
  mX0 <- as.matrix(X0t[, 1:(dim(X0t)[2]-2)])
  mX1 <- as.matrix(X1t[, 1:(dim(X1t)[2]-2)])
  mY0 <- as.matrix(y0)
  mY1 <- as.matrix(y1)
  
  pyX0   <- np_array(np$array(mX0), dtype = "float")
  pyX1   <- np_array(np$array(mX1), dtype = "float")
  pyY0   <- np$array(mY0, dtype = "float")
  pyY1   <- np$array(mY1, dtype = "float")
  
  res.optim2_0 <- tryCatch(gp(pyX0, pyY0,
                                 degree1 = degree1,
                                 degree2 = degree2,
                                 k1 = k1,
                                 k2 = k2,
                                 operator = operator),
                              error=function(e) NULL)
  
  res.optim2_1 <- tryCatch(gp(pyX1,pyY1,
                                 degree1 = degree1,
                                 degree2 = degree2,
                                 k1 = k1,
                                 k2 = k2,
                                 operator=operator),
                              error=function(e) NULL)
  
  #compute K
  matrix_eva <- as.matrix( confounders )
  # evaluation matrix
  res.optim2_1$par <- exp(res.optim2_1$gpr$kernel_$theta)
  res.optim2_0$par <- exp(res.optim2_0$gpr$kernel_$theta)
  # Gram matrices
  K1 <- res.optim2_1$gpr$kernel_(matrix_eva)
  K0 <- res.optim2_0$gpr$kernel_(matrix_eva)
  
  p1 <- as.numeric(res.optim2_1$gpr$predict(matrix_eva))
  p0 <- as.numeric(res.optim2_0$gpr$predict(matrix_eva))
  
  V <- rep(1/n,n)
  
  #Quadratic part
  I1 <- diag( t1 )
  I0 <- diag( t0 )
  I1KI1 <- I1%*%K1%*%I1
  I0KI0 <- I0%*%K0%*%I0
  
  KI1 <- diag(t1)%*%K1
  KI0 <- diag(t0)%*%K0
  
  VKI1 <- t(V)%*%KI1
  VKI0 <- t(V)%*%KI0
  
  tol <- 1e-08
  
  sigma1 <- res.optim2_1$par[3]^2
  sigma0 <- res.optim2_0$par[3]^2
  
  Sigma <- sigma1*diag(t1) + sigma0*diag(t0)
  #Update Q
  Q <- (1/n^2)*( I1KI1 + I0KI0 + penal*Sigma )
  
  #Update c
  c <- -2*(1/n^2)*(VKI1 + VKI0)
  
  rm(list = c("VKI1","VKI0"))
  
  model <- list()
  model$A          <- matrix(c( t1/n ,t0/n), nrow=2, byrow=T)
  model$rhs        <- c(1,1)
  model$modelsense <- "min"
  model$Q          <- Q
  model$obj        <- c
  model$sense      <- c("=")
  model$lb <- rep(tol,n)
  model$vtypes <- "C"
  Dmat <- Q  # Symmetric positive-definite matrix for the quadratic term
  dvec <- c  # Linear coefficients
  Amat <- t(matrix(c(t1/n, t0/n), nrow = 2, byrow = TRUE))
  bvec <- c(1, 1)  # Right-hand side values for the equality constraints
  meq <- 2  # N
  
  params <- list(Presolve=2,OutputFlag=0,QCPDual=0)
  
  res <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq)
  
  phi0 <- (1-data$Tr)*res$solution*(data$y - p0) + p0
  phi1 <- (data$Tr)*res$solution*(data$y - p1) + p1
  
  return(list(phi1 = phi1, phi0 = phi0))
  
  
  
}

aol_dgp <- function(n){
  # Copied from paper
  X <- lapply(1:5, function(i){
    dt <- data.frame(runif(n, -1, 1))
    names(dt) <- paste0('x', i)
    dt
  })
  X <- do.call(cbind, X)
  A <- ifelse(rbinom(n, size = 1, prob = 0.5) == 1, 1, -1)
  mu_y <- 0.5 + as.matrix(X) %*% c(0.5, 0.8, 0.3, -0.5, 0.7) + A*(0.2 - 0.6*X$x1 - 0.8*X$x2)
  y <- rnorm(n, mean = mu_y[, 1])
  X$y <- y
  X$A <- A
  return(X)
}

calculate_gamma <- function(n, subsets){
  soln <- 1 - log(subsets)/log(n)
  return(truncate_to_n(soln, 5))
}

causal_blb <- function(data, b, subsets, disjoint = TRUE, K = 10){
  if(data.table::is.data.table(data) == FALSE){
    data <- data.table::as.data.table(data)
  }
  n <- nrow(data)
  partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = disjoint)
  idx <- seq_len(b)
  
  blb_out <- lapply(partitions, function(i){
    tmp_dat <- data[i]
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

causal_blb_aipw <- function(data, b, subsets,  degree1, degree2, k1, k2, operator, penal, disjoint = TRUE){
  if(data.table::is.data.table(data) == FALSE){
    data <- data.table::as.data.table(data)
  }
  n <- nrow(data)
  partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = disjoint)
  idx <- seq_len(b)
  
  blb_out <- lapply(partitions, function(i){
    tmp_dat <- data[i]
    output <- kernel_weights(tmp_dat, degree1, degree2, k1, k2, operator, penal)
    phi1 <- output$phi1
    phi0 <- output$phi0
    
    M <- rmultinom(n = B, size = n, prob = rep(1, b))
    blb_reps <- sapply(seq_len(B), function(bt){
      boot_phi1 <- M[, bt]*phi1
      boot_phi0 <- M[, bt]*phi0
      sum(boot_phi1)/n - sum(boot_phi0)/n
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

truncate_to_n <- function(number, n) {
  factor <- 10^n
  truncated <- trunc(number * factor) / factor
  return(truncated)
}

zip_plots <- function(data, zip_labels, n, use_case, image_path, subsets = NULL, gamma = NULL){
  assertthat::assert_that(!((is.null((subsets) & !is.null(gamma)) | (!is.null(subsets) & is.null(gamma)))))
  if(data.table::is.data.table(data) == FALSE){
    data <- data.table::as.data.table(data)
  }

  nm_prefix <- paste0(use_case, '_zip_plot_n', n)
  if(is.null(subsets)){
    nm <- paste0(nm_prefix, '_full.pdf')
    title <- bquote(n == .(n))
    ggsub <- data[n == n]
    label_sub <- zip_labels[n == n]
    
  } else{
    nm <- paste0(nm_prefix, '_subset_', subsets, '_gamma_', gamma, '_cblb.pdf')
    title <- bquote(paste(s == .(subsets), ' and ', gamma == .(gamma), ' and ', n == .(n)))
    ggsub <- data[n == n & subsets == subsets & gamma == gamma]
    label_sub <- zip_labels[n == n & subsets == subsets & gamma == gamma]
  }
  
  p <- ggplot2::ggplot(ggsub, ggplot2::aes(y = rank)) +
    ggplot2::geom_segment(ggplot2::aes(x = lower_ci, y = rank, xend = upper_ci, yend = rank, color = covered)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = te), color = 'yellow', size = 0.5, linetype = 'dashed') +
    ggplot2::ylab('Fractional Centile of |z|') +
    ggplot2::xlab('95% Confidence Intervals') +
    ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(breaks = c(5, 50, 95)) +
    ggplot2::scale_color_discrete(name = "Coverage") +
    ggplot2::geom_text(x = 0.75, y = 50, ggplot2::aes(label = perc_cover), data = label_sub, size = 4) +
    ggplot2::ggtitle(title)
  
  if(is.null(subsets)){
    p <- p + ggplot2::facet_grid(n ~ ., labeller = ggplot2::label_both)
  } else{
    p <- p + ggplot2::facet_grid(n ~ subsets, labeller = ggplot2::label_both)
  }
  
  ggplot2::ggsave(file.path(image_path, nm), plot = p, height = 9, width = 7)
}

zip_plots_helper <- function(data, type){
  data[, `:=`(cent = abs(estim - te)/se)]
  if(type == 'full'){
    group_cols <- c('n')
  } else{
    group_cols <- c('n', 'gamma', 'subsets')
  }
  data[, `:=`(rank = as.integer(cut(cent, quantile(cent, probs = seq(0, 1, by = 0.01)), include.lowest = TRUE))), 
      by = group_cols]
  
  data[, `:=`(n = format(n, scientific = FALSE, big.mark = ','),
             covered = fifelse(lower_ci <= te & upper_ci >= te, 'Coverer', 'Non-coverer'))]
  zip_labels <- data[, .(perc_cover = round(mean(covered == 'Coverer'), 3)), by = group_cols]
  
  return(list(zip = data, zip_labels = zip_labels))
}

