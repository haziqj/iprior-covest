import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wishart

# Set parameters
m = 10  # Dimension of covariance matrix
B = 100  # Number of samples

# Simulate B covariance matrices from a Wishart distribution
# Using identity matrix for the scale matrix and m degrees of freedom
samples = [wishart.rvs(df=m, scale=np.eye(m)) for _ in range(B)]

# Compute average and variance across samples
average_matrix = np.mean(samples, axis=0)
variance_matrix = np.var(samples, axis=0)

# Visualization

# 1. Heatmap of Averages
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
sns.heatmap(average_matrix, cmap="coolwarm", cbar=True, annot=True, fmt=".2f")
plt.title("Average Covariance Matrix")

# 2. Heatmap of Variances
plt.subplot(1, 2, 2)
sns.heatmap(variance_matrix, cmap="coolwarm", cbar=True, annot=True, fmt=".2f")
plt.title("Variance of Covariance Matrix Entries")
plt.tight_layout()
plt.show()

# 3. Matrix Variability Plot
# Extract upper triangle of each matrix and visualize variability
upper_tri_indices = np.triu_indices(m, 1)
upper_tri_values = np.array([sample[upper_tri_indices] for sample in samples])

plt.figure(figsize=(14, 7))
for i in range(upper_tri_values.shape[1]):
    sns.lineplot(x=range(B), y=upper_tri_values[:, i], label=f"Entry {upper_tri_indices[0][i]+1},{upper_tri_indices[1][i]+1}");
    
plt.title("Variability of Matrix Entries Across Samples")
plt.xlabel("Sample Index")
plt.ylabel("Matrix Entry Value")
plt.legend(title="Matrix Entry", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
