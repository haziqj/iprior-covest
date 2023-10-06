library(iprior)
# source("10-data.R")
set.seed(123)

# Generate data
# m <- 10  
n <- 100
# x <- rchisq(m, df = 100)
# Theta_true <- solve(x %*% t(x) + diag(rep(0.0001, m)))
# # R <- tcrossprod(rnorm(m, sd = 0.5))
# Theta <- Theta_true #+ R + diag(rep(0.1, m))
# y <- mvtnorm::rmvnorm(n, rep(0, m), solve(Theta))

m <- 2
Sigma <- matrix(c(5, 2, 2, 10), nrow = m)
Theta <- solve(Sigma)
y <- mvtnorm::rmvnorm(n = n, sigma = Sigma)


# Posterior distribution of Theta
S <- crossprod(y - mean(y)) / n
W0 <- 1000 * diag(m)
nu <- m
What <- solve(W0) + n * S
nuhat <- nu + n
Thetahat <- nuhat * solve(What) 

# AR(1) process ----------------------------------------------------------------
# Parameters for the AR(1) process
phi <- 0.8
sigma2 <- 1
m <- 10

# Initialize the covariance matrix
R <- matrix(0, nrow = m, ncol = m)

# Compute the covariance values
for (i in 1:m) {
  for (j in 1:m) {
    R[i, j] <- (phi ^ abs(i - j)) * sigma2
  }
}

Sigma <- R
Theta <- solve(Sigma)
y <- mvtnorm::rmvnorm(n = n, sigma = Sigma)

# Posterior distribution of Theta
S <- crossprod(y - mean(y)) / n
W0 <- 1000 * diag(m)
nu <- m
What <- solve(W0) + n * S
nuhat <- nu + n
Thetahat <- nuhat * solve(What) 



