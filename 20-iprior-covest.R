library(tidyverse)
theme_set(theme_bw())
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

# I-prior process --------------------------------------------------------------
n <- 500
m <- 4
x <- runif(m, 1, 2)
H <- kern_fbm(x)
Theta_true <- x %*% t(x) + diag(x)
y <- mvtnorm::rmvnorm(n, rep(0, m), solve(Theta_true))

# write a function to get the diagonal and off-diagonals of W
vec_mat <- function(x) {
  n <- nrow(x)
  diag <- diag(x)
  offdiag <- x[upper.tri(x)]
  c(diag, offdiag)
}

rTheta <- function(Theta0, H, nu = nrow(Theta0), vectorise = FALSE) {
  W <- rWishart(1, nu, Theta0)[, , 1]
  Theta <- H %*% W %*% H
  if (isTRUE(vectorise)) {
    vec_mat(Theta)
  } else {
    Theta
  }
}

# Prior plots
B <- 250
map(1:B, \(x) tibble(i = x, Theta = list(rTheta(Theta_true, H, vec = TRUE)))) |>
  list_rbind() |>
  unnest(Theta) |>
  group_by(i) |>
  mutate(name = row_number()) |>
  ggplot(aes(name, Theta, group = i)) +
  geom_line(linewidth = 0.1, col = "gray") +
  scale_x_continuous(breaks = seq_len(m * (m + 1) / 2)) +
  theme(panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) 

# Posterior plots
B <- 250
S <- crossprod(y - mean(y)) / n
What <- solve(n * H %*% S %*% H + solve(S))
nuhat <- nrow(What) + n
map(1:B, \(x) tibble(i = x, 
                     Theta = list(rTheta(What, H, nuhat, TRUE)))) |>
  list_rbind() |>
  unnest(Theta) |>
  group_by(i) |>
  mutate(name = row_number()) |>
  ggplot(aes(name, Theta, group = i)) +
  geom_line(linewidth = 0.1, col = "gray") +
  geom_line(data = tibble(Theta = vec_mat(Theta_true), 
                          name = seq_len(m * (m + 1) / 2)), 
            col = "red", size = 1) +
  scale_x_continuous(breaks = seq_len(m * (m + 1) / 2)) +
  theme(panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) 








