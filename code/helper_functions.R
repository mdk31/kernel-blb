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

aol_loss <- function(params, X, A, r_tilde, K_matrix, lambda) {
  n <- length(A)
  v <- params[1:n]  # Extract v coefficients
  b <- params[n+1]   # Extract intercept
  
  f_x <- K_matrix %*% v + b
  
  hinge_loss <- huber_hinge(A * sign(r_tilde) * f_x)
  
  loss_term <- mean(abs(r_tilde) / 0.5 * hinge_loss)
  reg_term <- (lambda / 2) * sum(v %*% K_matrix %*% v)
  
  return(loss_term + reg_term)
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
    output <- aipw_kernel_weights(tmp_dat, degree1, degree2, k1, k2, operator, penal)
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

causal_blb_stable <- function(data, b, subsets, kernel_approx = TRUE, disjoint = TRUE, augment = FALSE,
                              delta.v = 1e-4){
  if(data.table::is.data.table(data) == FALSE){
    data <- data.table::as.data.table(data)
  }
  n <- nrow(data)
  partitions <- make_partition(n = n, subsets = subsets, b = b, disjoint = disjoint)
  idx <- seq_len(b)
  
  blb_out <- lapply(partitions, function(i){
    tmp_dat <- data[i]
    output <- osqp_kernel_sbw_twofit(X = as.matrix(tmp_dat[, c('X1', 'X2')]),
                                                   A = tmp_dat$Tr,
                                                   Y = tmp_dat$y,
                                                   delta.v=delta.v,
                                                   kernel.approximation = kernel_approx,
                                                   c = 100)
    
    weights0 <- output[[1]]$res0$x
    weights1 <- output[[1]]$res1$x
    tmp_dat$full_weights <- NA
    tmp_dat$full_weights[tmp_dat$Tr == 1] <- weights1
    tmp_dat$full_weights[tmp_dat$Tr == 0] <- weights0
    tmp_dat$full_weights <- tmp_dat$full_weights*b
    if(augment){
      m0 <- output[[1]]$m0
      m1 <- output[[1]]$m1
      
      phi1 <- (tmp_dat$Tr)*tmp_dat$full_weights*(tmp_dat$y - m1) + m1
      phi0 <- (1-tmp_dat$Tr)*tmp_dat$full_weights*(tmp_dat$y - m0) + m0
    } else{
      phi1 <- (tmp_dat$Tr)*tmp_dat$full_weights*(tmp_dat$y)
      phi0 <- (1-tmp_dat$Tr)*tmp_dat$full_weights*(tmp_dat$y)
    }

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

estimate_optimal_regime <- function(data, initial_params, lambda){
  
  num_params <- nrow(data)
  X <- data[, c('x1', 'x2', 'x3', 'x4', 'x5')]
  A <- data$A
  K_matrix <- kernelMatrix(vanilladot(), as.matrix(X))
  r_tilde <- train_aol(dat)
  
  opt_result <- optim(
    par = initial_params,
    fn = aol_loss,
    X = X,
    A = A,
    r_tilde = r_tilde,
    K_matrix = K_matrix,
    lambda = lambda,
    method = "L-BFGS-B"
  )
  
  # Extract optimized parameters
  v_opt <- opt_result$par[1:num_params]
  b_opt <- opt_result$par[num_params + 1]
  decision_boundary <- K_matrix %*% v_opt + b_opt
  estim_opt_regime <- ifelse(decision_boundary > 0, 1, -1)
  return(estim_opt_regime)
}

huber_hinge <- function(u, delta = 1) {
  ifelse(u >= 1, 0, ifelse(u >= -1, (1 - u)^2 / 4, -u))
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


kernel.basis <- function(X,A,Y, 
                         kernel.approximation=TRUE,
                         dim.reduction=FALSE,
                         c = NULL, l=NULL, s=NULL, gamma=NULL, U.mat = NULL) {
  n <- nrow(X)
  if (kernel.approximation) {
    if (is.null(c)) {
      # use some heuristics for choice of c 
      c = ifelse(n<1e+5, 100, 250) 
    }
    if (is.null(l)) {
      # set arbitrarily similar to Wang, 2019  
      l=round(c/2)
      s=round(l/2)
    }
    if (is.null(gamma)) {
      # heuristics by Hazlett, 2020  
      gamma = 1/(2*ncol(X)) 
    }
    
    set_c = sample(x = 1:n, size = c, replace = FALSE)
    set_c = sort(set_c)
    X <- scale(X)
    # --------------------------------------------------------
    # # slow version:  
    # f <- function(x,y) sum(x*y)
    # C <- outer(
    #   1:n, set_c,
    #   Vectorize( function(i,j) f(X[i,], X[j,]) )
    # )
    # --------------------------------------------------------
    # C <- RBF_kernel_C(X, c, set_c)
    C <- RBF_kernel_C_parallel(X, c, set_c)
    # --------------------------------------------------------
    W <- C[set_c,]
    SVD_W <- RSpectra::svds(A=W, k=l)
    if (dim.reduction) {
      R <- C %*% SVD_W$u %*% diag(1/sqrt(SVD_W$d))
      SVD_R <- svds(A=R, k=s) 
      X_ <- R %*% SVD_R$v # B
    } else {
      X_ <- C %*% SVD_W$u # %*% diag(1/SVD_W$d)
      # lambda_ <- 1/SVD_W$d
    }
  } else{
    # O(n^2) spoce-, O(n^3) time-complexity
    gram.mat <- makeK_noparallel(X)
    res = RSpectra::eigs_sym(gram.mat, c, which = "LM")
    # CLIP NEGATIVE
    # res$values[res$values < 1e-10] <- 1e-10
    X_ <- if(is.null(U.mat)) {
      res$vectors %*% diag(1/sqrt(res$values))
    } else{
      U.mat %*% diag(1/sqrt(res$values))
    }
  }
  return(X_)
}

makeK_noparallel <- function(allx, useasbases=NULL, b=NULL, linkernel = FALSE, scale = TRUE){
  N=nrow(allx)
  # If no "useasbasis" given, assume all observations are to be used.
  if(is.null(useasbases)) {useasbases = rep(1, N)}
  
  #default b is set to 2ncol to match kbal for now
  if (is.null(b)){ b=2*ncol(allx) }
  
  if(scale) {
    allx = scale(allx)
  } 
  bases = allx[useasbases==1, ]
  
  if (linkernel == TRUE) {
    K <- allx  # Linear kernel case
  } else {
    if (sum(useasbases) == N) {
      K <- kernlab::kernelMatrix(kernlab::vanilladot(), allx)
    } 
  }
  return(K)
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

# From JOSE paper
osqp_kernel_sbw_twofit <- function(X,A,Y,
                                   delta.v=0.005, 
                                   X_=NULL,
                                   osqp.setting=NULL,
                                   basis="kernel", kernel.approximation=TRUE,
                                   c = NULL, l=NULL, gamma=NULL, U.mat = NULL,
                                   dim.reduction=FALSE, s=NULL,
                                   K=2, interactions=FALSE) {
  
  #------------------------------------------------------------------------------------------
  # @X: covariates (n x d matrix)
  # @A: treatment (n x 1 vector)
  # @Y: outcome (n x 1 vector)
  # @X_: the basis matrix where (approximate) mean balancing between the treated and control groups is to be achieved
  #      (ncol(X_) = n)
  # @osqp.setting: osqp settings
  # @delta.v: tolerance level (either a single number or a numeric vector);
  #           the numeric vector is used for multiple iterations, each with different tolerance level
  # @basis: support kernel- or power series- (moment-) based constructions
  # + kernel basis (basis=='kernel')
  #     @kernel.approximation: if TRUE, use the Nystrom kernel approximation method proposed in \insertRef{wang2019scalable}
  #                            if FALSE, use the full gram matrix to compute eigenvectors as in \insertRef{hazlett2018}
  #     @c: sketch size for the standard Nystrom approximation. 
  #     @l: regularization parameter l < c. 
  #     @s: target rank for s < l the rank-restricted Nyström approximation
  #     @gamma: RBK parameter; set default value as 1/(2*dim(X))
  #     @U.mat: externel input for nxc eigenvector matrix; only used when kernel.approximation==FALSE
  # + power series basis (basis=='power')
  #     @K: use up to the K-th moment
  #     @interactions: add up to the K-th order interaction terms if TRUE
  #------------------------------------------------------------------------------------------
  
  
  res.list <- list()
  if (is.null(X_)) {
    if (basis=="kernel"){
      X_ <- kernel.basis(X,A,Y, 
                         kernel.approximation=kernel.approximation, 
                         c=c, l=l, gamma=gamma, U.mat=U.mat,
                         s=s, dim.reduction=dim.reduction)
    } 
    
    else if (basis=="power") {
      X_ <- power.basis(X,A,Y, 
                        K=K, interactions=interactions)
    } 
    
    else {
      stop("not available yet")
    }
  }
  
  nX <- ncol(X_)
  n1 <- sum(A); n0 <- sum(1-A)
  Xt <- X_[A==1,]; Xc <- X_[A==0,]
  Yt <- Y[A==1]; Yc <- Y[A==0]
  
  # Xc_df <- as.data.frame(Xc)
  # Xt_df <- as.data.frame(Xt)
  # model_c <- lm(Yc ~ ., data = Xc_df)
  # model_t <- lm(Yt ~ ., data = Xt_df)
  # 
  # X_df <- as.data.frame(X_)
  # colnames(X_df) <- colnames(Xt_df)
  # 
  # m0 <- predict(model_c, newdata = X_df)
  # m1 <- predict(model_t, newdata = X_df)
  m0 <- 0
  m1 <- 1
  
  P.mat0 <- as(Matrix::.symDiagonal(n=n0, x=1.), "dgCMatrix")
  q.vec0 <- rep(-1./n0,n0)
  
  P.mat1 <- as(Matrix::.symDiagonal(n=n1, x=1.), "dgCMatrix")
  q.vec1 <- rep(-1./n1,n1)
  
  A.mat0 <- Matrix::Matrix(rbind(rep(1.,n0),
                                 P.mat0,
                                 t(Xc)), sparse = TRUE)
  A.mat1 <- Matrix::Matrix(rbind(rep(1.,n1),
                                 P.mat1,
                                 t(Xt)), sparse = TRUE)
  
  if (is.null(osqp.setting)) {
    # use default one
    settings <- osqp::osqpSettings(alpha = 1.5, verbose = FALSE)  
  } else {
    settings <- osqp.setting
  }
  
  for (j in 1:length(delta.v)) {
    l.vec0 <- c(1., rep(0.,n0),
                colMeans(Xt) - delta.v[j] * rep(1,nX))
    u.vec0 <- c(1., rep(1.,n0),
                colMeans(Xt) + delta.v[j] * rep(1,nX))
    
    l.vec1 <- c(1., rep(0.,n1),
                colMeans(Xc) - delta.v[j] * rep(1,nX))
    u.vec1 <- c(1., rep(1.,n1),
                colMeans(Xc) + delta.v[j] * rep(1,nX))
    if (j==1) {
      model0 <- osqp::osqp(P.mat0, q.vec0, A.mat0, l.vec0, u.vec0, settings)
      model1 <- osqp::osqp(P.mat1, q.vec1, A.mat1, l.vec1, u.vec1, settings)
    } else {
      model$Update(l = l.vec, u = u.vec)
    }
    res0 <- model0$Solve()
    res1 <- model1$Solve()
    # if (res$info$status != "solved") {
    #   warning(res$info$status)
    # }
    res.list[[j]] <- list(res0 = res0, res1 = res1, m0 = m0, m1 = m1)
  }
  return(res.list)
}

truncate_to_n <- function(number, n) {
  factor <- 10^n
  truncated <- trunc(number * factor) / factor
  return(truncated)
}

zip_plots <- function(data, zip_labels, n, use_case, image_path, subsets = NULL, gamma = NULL,
                      text_x = 0.75, text_y = 50){
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
    ggplot2::geom_text(x = text_x, y = text_y, ggplot2::aes(label = perc_cover), data = label_sub, size = 4) +
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

