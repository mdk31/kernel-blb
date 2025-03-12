library(pbapply)
library(data.table)
library(sbw)

source('code/helper_functions.R')
Rcpp::sourceCpp("code/RBF_kernel_C_parallel.cpp")

te <- 0.8
sigma <- 1
replications <- 1000
B <- 100 # bootstrap resamples

dat <- kangschafer3(n = n, te = te, sigma = sigma, beta_overlap = 0.5)


weights0 <- output[[1]]$res0$x
weights1 <- output[[1]]$res1$x
dat$full_weights <- NA
dat$full_weights[dat$Tr == 1] <- weights1
dat$full_weights[dat$Tr == 0] <- weights0
dat$full_weights <- dat$full_weights*n
# dat$full_weights <- output[[1]]$res$x*n

phi1 <- (dat$Tr)*dat$full_weights*(dat$y)
phi0 <- (1-dat$Tr)*dat$full_weights*(dat$y)
