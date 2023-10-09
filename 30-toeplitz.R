library(tidyverse)
theme_set(theme_bw())
library(iprior)
set.seed(123)

gen_cov_mat <- function(m, x = 1) {
  rho <- function(k, x) {
    return(exp(-k / x))
  }
  # Create the first column of the Toeplitz covariance matrix
  sigma_column <- sapply(0:(m-1), function(k) rho(k, x))
  
  # Generate the Toeplitz covariance matrix
  covariance_matrix <- toeplitz(sigma_column)
  
  return(covariance_matrix)
}

m <- 5
Sigma0 <- gen_cov_mat(m)
Theta0 <- solve(Sigma0)
n <- 10000
y <- mvtnorm::rmvnorm(n = n, sigma = Sigma0)
x <- 0:(m-1)
H <- kern_fbm(x)

# Posterior mean of I-prior estimate for Theta
S <- crossprod(y - mean(y)) / n
What <- solve(n * H %*% S %*% H + solve(Theta0))
nuhat <- nrow(What) + n
Theta_hat <- H %*% (What * nuhat) %*% H

round(Theta0, 3)
round(Theta_hat, 3)
