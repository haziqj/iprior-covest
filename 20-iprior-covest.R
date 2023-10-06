library(iprior)
# source("10-data.R")
set.seed(123)

# Generate data
m <- 10  
n <- 250  
x <- rchisq(m, df = 100)
Theta_true <- solve(x %*% t(x) + diag(rep(0.0001, m)))
# R <- tcrossprod(rnorm(m, sd = 0.5))
Theta <- Theta_true #+ R + diag(rep(0.1, m))
y <- mvtnorm::rmvnorm(n, rep(0, m), solve(Theta))

# Posterior distribution of Theta
S <- crossprod(y - mean(y)) / n
W0 <- 1000 * diag(m)
nu <- m
What <- solve(W0) + n * S
nuhat <- nu + n
Thetahat <- nuhat * solve(What) 
