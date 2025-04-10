# Load necessary library
library(ggplot2)

# Function to generate a sample partition from the Gnedin Process
gnedin_process <- function(n, gamma = 0.5) {
  clusters <- c(1)  # First element starts in a new cluster
  num_clusters <- 1  # Track the number of clusters
  
  for (i in 2:n) {
    # Compute probabilities for joining existing clusters
    existing_probs <- sapply(1:num_clusters, function(k) {
      sum(clusters == k) * (i - num_clusters + gamma)
    })
    
    # Compute probability for forming a new cluster
    new_cluster_prob <- num_clusters^2 - num_clusters * gamma
    
    # Normalize probabilities
    probs <- c(existing_probs, new_cluster_prob)
    probs <- probs / sum(probs)
    
    # Sample a cluster assignment
    new_assignment <- sample(1:(num_clusters + 1), size = 1, prob = probs)
    
    # Assign the element to the selected cluster
    clusters <- c(clusters, new_assignment)
    
    # If a new cluster is created, update counter
    if (new_assignment > num_clusters) {
      num_clusters <- new_assignment
    }
  }
  
  return(clusters)
}

# Generate and visualize multiple draws from the Gnedin process
set.seed(6)
n_samples <- 100  # Number of elements
num_draws <- 5    # Number of different partitions

# Create a data frame to store draws
df <- data.frame(
  Element = rep(1:n_samples, num_draws),
  Partition = unlist(lapply(1:num_draws, function(x) gnedin_process(n_samples, gamma = 1.0))),
  Draw = rep(1:num_draws, each = n_samples)
)

# Plot the partitions
ggplot(df, aes(x = Element, y = Draw, color = factor(Partition))) +
  geom_point(size = 3) +
  labs(title = "",
       x = "Element Index",
       y = "Sampled Partition",
       color = "Cluster") +
  theme_minimal() +
  scale_color_brewer(palette = "Set3")


### DP 

# Function to generate a sample partition from the Dirichlet Process (DP)
dirichlet_process <- function(n, alpha = 1.0) {
  clusters <- c(1)  # First element starts in a new cluster
  num_clusters <- 1  # Track the number of clusters
  
  for (i in 2:n) {
    # Compute probabilities for joining existing clusters
    existing_probs <- sapply(1:num_clusters, function(k) {
      sum(clusters == k)
    })
    
    # Compute probability for forming a new cluster
    new_cluster_prob <- alpha
    
    # Normalize probabilities
    probs <- c(existing_probs, new_cluster_prob)
    probs <- probs / sum(probs)
    
    # Sample a cluster assignment
    new_assignment <- sample(1:(num_clusters + 1), size = 1, prob = probs)
    
    # Assign the element to the selected cluster
    clusters <- c(clusters, new_assignment)
    
    # If a new cluster is created, update counter
    if (new_assignment > num_clusters) {
      num_clusters <- new_assignment
    }
  }
  
  return(clusters)
}

# Generate and visualize multiple draws from the Dirichlet Process
set.seed(123)
n_samples <- 100  # Number of elements
num_draws <- 5    # Number of different partitions

# Create a data frame to store draws
df <- data.frame(
  Element = rep(1:n_samples, num_draws),
  Partition = unlist(lapply(1:num_draws, function(x) dirichlet_process(n_samples, alpha = 0.5))),
  Draw = rep(1:num_draws, each = n_samples)
)

# Plot the partitions
ggplot(df, aes(x = Element, y = Draw, color = factor(Partition))) +
  geom_point(size = 3) +
  labs(title = "",
       x = "Element Index",
       y = "Sampled Partition",
       color = "Cluster") +
  theme_minimal() +
  scale_color_brewer(palette = "Set3")

## bayesian hist 

# Load required packages
# Load required packages
library(ggplot2)
library(dplyr)
library(MCMCpack) # For Dirichlet distribution

# Step 1: Generate synthetic data (Mixture of Gaussians)
set.seed(42)
n <- 500
data <- c(rnorm(n/2, mean = -2, sd = 1), rnorm(n/2, mean = 2, sd = 1.2))
data_disc <- c(rbinom(n, 1:num_bins, 0.5))

# Step 2: Define histogram bins
num_bins <- 20
bin_edges <- seq(min(data), max(data), length.out = num_bins + 1)

# Step 3: Assign prior probabilities using a Dirichlet distribution
pi_prior <- rep(1, num_bins) # Non-informative prior (uniform)
dirichlet_sample <- rdirichlet(1, pi_prior) # Sample from Dirichlet prior

# Step 4: Compute histogram counts
bin_counts <- hist(data, breaks = bin_edges, plot = FALSE)$counts

# Step 5: Compute Bayesian posterior by updating the Dirichlet prior
posterior_pi <- bin_counts + pi_prior # Update posterior parameters
posterior_probs <- as.vector(rdirichlet(1, posterior_pi)) # Sample posterior probabilities

# Step 6: Create a data frame for visualization
hist_df <- data.frame(
  bin_midpoints = (bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2,
  posterior_density = posterior_probs * length(data) / diff(bin_edges) # Normalize
)

# Step 7: Visualize Bayesian Histogram Model
ggplot(hist_df, aes(x = bin_midpoints, y = posterior_density)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "",
       x = "Data Value",
       y = "Density") +
  theme_minimal()


### disc data bayesian histogram ex 
set.seed(123)

n_trials <- 10
data_disc <- c(rbinom(n, n_trials, 0.5))

# Step 3: Assign prior probabilities using a Dirichlet distribution
pi_prior <- rep(1, n_trials) # Non-informative prior (uniform)
dirichlet_sample <- rdirichlet(1, pi_prior) # Sample from Dirichlet prior

# Step 4: Compute histogram counts
bins = 1:n_trials
bin_counts <- table(factor(data_disc, levels = bins))


# Step 5: Compute Bayesian posterior by updating the Dirichlet prior
posterior_pi_disc <- bin_counts + pi_prior # Update posterior parameters
posterior_probs_disc <- as.vector(rdirichlet(1, posterior_pi_disc)) # Sample posterior probabilities

# Step 6: Create a data frame for visualization
hist_df_disc <- data.frame(
  bin = bins,
  posterior_density = posterior_probs_disc
)

# Step 7: Visualize Bayesian Histogram Model
ggplot(hist_df_disc, aes(x = bin, y = posterior_density)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "",
       x = "Data Value",
       y = "Density") +
  theme_minimal()

