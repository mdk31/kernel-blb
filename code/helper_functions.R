

crossfit_estimator <- function(data, K = 10){
  assertthat::assert_that(is.integer(K))
  if(data.table::is.data.table(data) == FALSE){
    data <- data.table::as.data.table(data)
  }
  
  fold_idx <- seq_len(nrow(data))
  folds <- split(fold_idx, sample(rep(1:K, length.out = length(fold_idx))))
  
  crossfit_dt <- lapply(folds, function(test_idx){
    train_idx <- setdiff(fold_idx, test_idx)
    train_dat <- data[train_idx]
    test_dat <- data[-train_idx]
    
    m <- svm(y ~ Tr + X1 + X2, data = train_dat, kernel = 'linear', type = 'eps-regression')
    g <- svm(Tr ~ X1 + X2, data = train_dat, kernel = 'linear', type = 'C-classification', probability = TRUE)
    
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
  return(rbindlist(crossfit))
}


kangschafer3 <- function(){
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

