library(tidyverse)
set.seed(42)

# Define parameters
m <- 10  # number of spatial units
n <- 250  # sample size

# Generate a covariate x that affects the covariance matrix
x <- runif(m, 10, 20)  # Values between 10 and 20

# Generate the covariance matrix using xx^T
Theta <- x %*% t(x)
R <- tcrossprod(rnorm(m, sd = 2))
Theta <- Theta + R
Theta <- Theta + diag(rep(0.1, m))

# Sample n observations from this covariance matrix
ys <- mvtnorm::rmvnorm(n, rep(0, m), Theta)

# Bayesian modeling with Inverse-Wishart prior
# Prior parameters
nu_prior <- m + 2  # degrees of freedom
Psi_prior <- diag(m)  # scale matrix, identity for simplicity

# Posterior parameters
nu_posterior <- nu_prior + n
sample_mean <- colMeans(ys)
S <- Reduce(`+`, lapply(1:n, function(i) (ys[i,] - sample_mean) %*% t(ys[i,] - sample_mean)))
Psi_posterior <- Psi_prior + S

# Compute the posterior mean for Theta
Theta_posterior_mean <- Psi_posterior / (nu_posterior + m + 1)

# Compute fit index: Frobenius norm between the true covariance matrix and the posterior mean
fit_index <- norm(as.matrix(Theta - Theta_posterior_mean), type="F")

# Plotting functions (similar to Python heatmaps)
plot_heatmap <- function(mat, title) {
  df <- as.data.frame(as.table(mat))
  ggplot(df, aes(Var1, Var2)) + 
    geom_tile(aes(fill = Freq), color = "white") +
    geom_text(aes(label = sprintf("%.2f", Freq)), vjust = 1) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = mean(df$Freq), limit = c(min(df$Freq), max(df$Freq)), 
                         space = "Lab", name="Value") + 
    theme_minimal() + 
    labs(title = title, x = "Neighborhood", y = "Neighborhood") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Plot the heatmaps
plot1 <- plot_heatmap(Theta, "True Covariance Matrix")
plot2 <- plot_heatmap(Theta_posterior_mean, "Posterior Mean")

# Display plots side by side
library(gridExtra)
grid.arrange(plot1, plot2, ncol=2)

fit_index




library(ggplot2)
library(reshape2)

# Set parameters
m <- 10  # Dimension of covariance matrix
B <- 100  # Number of samples

# Simulate B covariance matrices from a Wishart distribution
samples <- lapply(1:B, function(i) {
  stats::rWishart(1, m, diag(m))
})

# Compute average and variance across samples
average_matrix <- Reduce("+", samples) / B
variance_matrix <- Reduce("+", lapply(samples, function(mat) {
  (mat - average_matrix)^2
})) / B

# 1. Heatmap of Averages
ggplot(melt(average_matrix), aes(Var1, Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(name="Average Value") + 
  labs(title="Average Covariance Matrix", x="Row", y="Column") +
  theme_minimal()

# 2. Heatmap of Variances
ggplot(melt(variance_matrix), aes(Var1, Var2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(name="Variance") + 
  labs(title="Variance of Covariance Matrix Entries", x="Row", y="Column") +
  theme_minimal()

# 3. Matrix Variability Plot
upper_tri_indices <- which(upper.tri(diag(m)))
upper_tri_values <- do.call(rbind, lapply(samples, function(mat) {
  as.numeric(mat)[upper_tri_indices]
}))

df <- as.data.frame(upper_tri_values)
# colnames(df) <- sprintf("Entry %s,%s", row(average_matrix)[upper_tri_indices], col(average_matrix)[upper_tri_indices])

df_melted <- melt(mutate(df, id = row_number()), id.vars = "id")
ggplot(df_melted, aes(x=variable, y=value, group = id)) +
  geom_line(aes(color=id)) + 
  labs(title="Variability of Matrix Entries Across Samples", x="Sample Index", y="Matrix Entry Value") +
  theme_minimal() + 
  theme(legend.position="top")
