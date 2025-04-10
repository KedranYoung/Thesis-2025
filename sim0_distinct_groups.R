
setwd("/Users/kedranyoung/Desktop/THESIS 2025/Code/Simulation Studies")
source("Source/code_MFMSBM.R") # from Geng's Github
source("Source/esbm.R") # from Durante's Github
source("Source/traditional_SBM.R")
Rcpp::sourceCpp('Source/stirling.cpp') # from Durante's Github

library(LaplacesDemon)
library(igraph)
library(reshape)
library(gdata)
#library(mcclust.ext)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(coda)
#install.packages("/Users/kedranyoung/Desktop/THESIS 2025/Code/Durante Code/dummies_1.5.6.tar.gz", repos = NULL, type = "source")
library(dummies)
library(randnet)
library(greed)
library(RColorBrewer)
library(scales)
library(clue)
library(dplyr)
library(extrafont)
library(ggplot2)
library(reshape2)

load("simulated_data0.RData")

# ------------------------------------

# ========== Visualization ===========

# ------------------------------------

graph_sim2 <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)

V(graph_sim2)$z_0 <- z_0
group_colors <- brewer.pal(length(unique(z_0)), "Dark2") 
V(graph_sim2)$color <- group_colors[as.factor(V(graph_sim2)$z_0)]
V(graph_sim2)$frame.color <- "black"  
E(graph_sim2)$color <- alpha("gray50", 0.6)  

plot(graph_sim2, 
     vertex.size = 10,  
     vertex.label = NA, 
     vertex.color = V(graph_sim2)$color,
     vertex.frame.color = V(graph_sim2)$frame.color,  
     edge.width = 1,  
     edge.arrow.size = 0.5,
     layout = layout_with_fr,  
     main = "Community Structure: Simulated Data (2)")

legend("topright", legend = levels(as.factor(V(graph_sim2)$z_0)), 
       col = group_colors, pch = 21, pt.bg = group_colors,  
       pt.cex = 2, cex = 1.2, bty = "n", text.col = "black") 

library(kernlab)  # for specc
k <- 4  # number of blocks

# Use spectral clustering to reorder the adjacency matrix
clusters <- specc(A, centers = k)
ordering <- order(clusters)
A_reordered <- A[ordering, ordering]

Y = A_reordered
V = N = dim(A)[1]
diag(Y) <- 0

z_0_reordered <- z_0[ordering]
row_plot_Y <- data.frame(z_0 = as.factor(z_0_reordered))
names(row_plot_Y) <-"z_0"
rownames(Y) <- rownames(row_plot_Y)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_Y$z_0)
mycolors <- list(z_0 = mycolors)

Network <- pheatmap(Y,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_Y,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)


# ------------------------------------

# ============= MFM SBM ==============

# ------------------------------------

niters <- 10000
burn <- 1000
N <- dim(A)[1]

fit1 = CDMFM_new(data = A, data1 = AAA, niterations = niters, beta.a = 1, beta.b = 1, 
                 GAMMA=1, LAMBDA = 1.5, initNClusters = 10)

# CDMFM_new() found in code_MFMSBM.R (Geng)

# ======== posterior, ratios ===========

#load("sim0_results.RData")

Z_MFM <- matrix(0, (niters - burn), N)  

for(i in (burn + 1):niters){  
  Z_MFM[i - burn, ] <- fit1$Iterates[[i]]$zout 
  Z_MFM[i - burn, ] <- sort(Z_MFM[i - burn, ])  
}

sorted_z0 <- sort(z_0, decreasing = FALSE)
result2 <- round(colMeans(Z_MFM))

# will change depending on output, to fix visual analysis due to label switching
result2_reassigned <- recode(result2, `1` = 4, `2` = 1, `3` = 3, `4` = 2)
result2 <- sort(result2_reassigned)
table(result2)

result1_reassigned <- recode(result2, `1` = 4, `2` = 1, `3` = 3, `4` = 2)
table(result1_reassigned)

table(sorted_z0)

# cluster confusion matrices to visualize z_hat output
conf_matrix <- as.data.frame(table(sorted_z0, result2))

true_labels <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
names(true_labels) <- c(1, 2, 3, 4)

confusion_plot_MFM <- ggplot(conf_matrix, aes(x = factor(sorted_z0, levels = c(1, 2, 3, 4)), 
                                              y = factor(result2, levels = rev(levels(factor(result2)))), 
                                              fill = Freq)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "lightblue") +
  scale_x_discrete(labels = true_labels) +
  labs(
    x = "True Clusters",
    y = "Inferred Clusters"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )


ggsave("confusion_matrix_mfm.pdf", 
       plot = confusion_plot_MFM, 
       width = 7, 
       height = 5, 
       units = "in")

# ============= Posterior checking and model comparison ==============

# WAIC
set.seed(1)
Z_MFM_WAIC <- matrix(0, (niters - burn), N)  
for(i in (burn + 1):niters){  
  Z_MFM_WAIC[i - burn, ] <- fit1$Iterates[[i]]$zout
}

Z_MFM_WAIC <- t(Z_MFM_WAIC)

V = N
a <- b <- 1
LL <- matrix(nrow=V*(V-1)/2,ncol=9000)

for (t in 1:dim(Z_MFM_WAIC)[2]){
  LL[,t]<-sampleLL(as.factor(Z_MFM_WAIC[,t]),A,a,b)
  if (t%%10000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC

quantile(apply(t(Z_MFM)[,1:9000],2,max))[c(2:4)]

# TRACE PLOT
# randomly select
set.seed(1)
index_traceplot <- sample(c(1:(V*(V-1)/2)),1)
plot(ts(LL[index_traceplot,]),xlab="",ylab="")

# HISTOGRAM OF EDGES 
sampleEdgeCount <- function(memb, Y, a, b){
  z <- dummy(memb)
  H <- ncol(z)
  V <- dim(Y)[1]
  
  M <- t(z)%*%Y%*%z
  diag(M) <- diag(M)/2
  Tot <- t(z)%*%matrix(1,V,V)%*%z
  diag(Tot) <- (diag(Tot)-table(memb))/2
  Mbar <- Tot - M
  a_n <- lowerTriangle(M, diag=TRUE) + a
  b_bar_n <- lowerTriangle(Mbar, diag=TRUE) + b
  
  theta <- rbeta(length(a_n), a_n, b_bar_n)
  Theta <- matrix(0, H, H)
  Theta[lower.tri(Theta, diag=TRUE)] <- theta
  Theta <- Theta + t(Theta)
  diag(Theta) <- diag(Theta)/2
  
  edge_prob <- z %*% Theta %*% t(z)
  diag(edge_prob) <- NA
  
  # Return total expected number of edges (undirected, no self-loops)
  return(sum(edge_prob[lower.tri(edge_prob)]))
}

expected_edges <- numeric(dim(Z_MFM_WAIC)[2])

for (t in 1:dim(Z_MFM_WAIC)[2]) {
  expected_edges[t] <- sampleEdgeCount(as.factor(Z_MFM_WAIC[,t]), A, a, b)
  if (t %% 1000 == 0) print(paste("Iteration:", t))
}

observed_edges <- sum(A[lower.tri(A)])

# Histogram
hist(expected_edges,
     breaks = 30,
     col = "lightblue",
     border = "white",
     main = "",
     xlab = "Predicted Number of Edges")

abline(v = observed_edges, col = "red", lwd = 2)

# ------------------------------------

# ============= GN ==============

# ------------------------------------

Y = A
N_iter <- 10000
V <- dim(Y)[1]
my_seed <- 123
my_z <- c(1:V)

gamma <- 0.3
probs_gnedin <- HGnedin(V, 1:V, gamma = gamma)
round(sum(1:V*probs_gnedin))

my_prior <- "GN"
Z_GN <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, gamma_GN = 0.3)

# esbm() found in esbm.R (Durante)

Z_mat <- matrix(0, N_iter, N)

for(i in (burn + 1):N_iter){  
  Z_mat[i - burn, ] <- Z_GN[,i] 
  Z_mat[i - burn, ] <- sort(Z_mat[i - burn, ])  
}

result_gn <- round(colMeans(Z_mat))

table(result_gn)

resultgn_reassigned <- recode(result_gn, `1` = 4, `2` = 1, `3` = 3, `4` = 2)
result_gn <- sort(resultgn_reassigned)
table(result_gn)

# cluster confusion matrices to visualize z_hat output
conf_matrix <- as.data.frame(table(sorted_z0, result_gn))

true_labels <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
names(true_labels) <- c(1, 2, 3, 4)

# Enhanced plot
confusion_plot_gn <- ggplot(conf_matrix, aes(x = factor(sorted_z0, levels = c(1, 2, 3, 4)), 
                                             y = factor(result_gn, levels = rev(levels(factor(result_gn)))), 
                                             fill = Freq)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "lightblue") +
  scale_x_discrete(labels = true_labels) +
  labs(
    x = "True Clusters",
    y = "Inferred Clusters"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )


ggsave("confusion_matrix_gn.pdf", 
       plot = confusion_plot_gn, 
       width = 7, 
       height = 5, 
       units = "in")

# ============= Posterior checking and model comparison ==============

# WAIC
set.seed(1)
V <- dim(Y)[1]
burn_in <- 1000
a <- b <- 1
LL <- matrix(nrow=V*(V-1)/2,ncol=9000)

Z_GN_WAIC <- Z_GN[,(burn_in+1):N_iter]

for (t in 1:dim(Z_GN_WAIC)[2]){
  LL[,t]<-sampleLL(Z_GN_WAIC[,t],Y,a,b)
  if (t%%10000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC

quantile(apply(Z_GN[,(burn_in+1):N_iter],2,max))[c(2:4)]

# TRACE PLOT
set.seed(1)
index_traceplot <- sample(c(1:(V*(V-1)/2)),1)
plot(ts(LL[index_traceplot,]),xlab="",ylab="")

# HISTOGRAM OF EDGES 
expected_edges <- numeric(dim(Z_GN_WAIC)[2])

for (t in 1:dim(Z_MFM_WAIC)[2]) {
  expected_edges[t] <- sampleEdgeCount(as.factor(Z_GN_WAIC[,t]), A, a, b)
  if (t %% 1000 == 0) print(paste("Iteration:", t))
}

observed_edges <- sum(A[lower.tri(A)])

# Histogram
hist(expected_edges,
     breaks = 30,
     col = "lightblue",
     border = "white",
     main = "",
     xlab = "Predicted Number of Edges")

abline(v = observed_edges, col = "red", lwd = 2)

# ------------------------------------

# ============= DP ==============

# ------------------------------------

N_iter <- 10000
N = 30
sigma_dp <- 0   
H_dp <- Inf 
alpha_dp <- 8
round(expected_cl_py(V, sigma = sigma_dp, theta = alpha_dp, H = H_dp))

my_prior <- "DP"
Z_DP <- esbm(Y, my_seed, N_iter, my_prior, my_z, a = 1, b = 1, alpha_PY = 8, sigma_PY = 0)

Z_mat <- matrix(0, N_iter, N)

for(i in (burn + 1):N_iter){  
  Z_mat[i - burn, ] <- Z_DP[,i] 
  Z_mat[i - burn, ] <- sort(Z_mat[i - burn, ])  
}

result_dp <- round(colMeans(Z_mat))

table(result_dp)

resultdp_reassigned <- recode(result_dp, `1` = 4, `2` = 1, `3` = 3, `4` = 2, `5` = 5)
result_dp <- sort(resultdp_reassigned)
table(result_dp)

# cluster confusion matrices to visualize z_hat output
conf_matrix <- as.data.frame(table(sorted_z0, result_dp))

true_labels <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Extra Cluster")
names(true_labels) <- c(1, 2, 3, 4, 5)  # Add label for extra cluster

confusion_plot_dp <- ggplot(conf_matrix, aes(x = factor(sorted_z0, levels = c(1, 2, 3, 4, 5)), 
                                             y = factor(result_dp, levels = rev(levels(factor(result_dp)))), 
                                             fill = Freq)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "lightblue") +
  scale_x_discrete(labels = true_labels) +  # Ensure it includes "Extra Cluster"
  labs(
    x = "True Clusters",
    y = "Inferred Clusters"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )


ggsave("confusion_matrix_dp.pdf", 
       plot = confusion_plot_dp, 
       width = 7, 
       height = 5, 
       units = "in")

# ============= Posterior checking and model comparison ==============

# WAIC 
set.seed(1)
V <- dim(Y)[1]
burn_in <- 1000
a <- b <- 1
LL <- matrix(nrow=V*(V-1)/2,ncol=9000)

Z_DP_WAIC <- Z_DP[,(burn_in+1):N_iter]

for (t in 1:dim(Z_DP_WAIC)[2]){
  LL[,t]<-sampleLL(as.factor(Z_DP_WAIC[,t]),Y,a,b)
  if (t%%10000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC

quantile(apply(Z_DP[,(burn_in+1):N_iter],2,max))[c(2:4)]

# HISTOGRAM OF EDGES 
expected_edges <- numeric(dim(Z_DP_WAIC)[2])

for (t in 1:dim(Z_DP_WAIC)[2]) {
  expected_edges[t] <- sampleEdgeCount(as.factor(Z_DP_WAIC[,t]), A, a, b)
  if (t %% 1000 == 0) print(paste("Iteration:", t))
}

observed_edges <- sum(A[lower.tri(A)])

# Histogram
hist(expected_edges,
     breaks = 30,
     col = "lightblue",
     border = "white",
     main = "",
     xlab = "Predicted Number of Edges")

abline(v = observed_edges, col = "red", lwd = 2)


# ------------------------------------

# ======== Traditional SBM =========

# ------------------------------------

set.seed(1)
fit_sbm <- bayesian_sbm_gibbs(A = A, K = 4)
fit_sbm2 <- bayesian_sbm_gibbs(A = A, K = 2)
fit_sbm6 <- bayesian_sbm_gibbs(A = A, K = 6)

table(fit_sbm$z_hat)
result_reassigned <- recode(fit_sbm$z_hat, `1` = 4, `2` = 1, `3` = 2, `4` = 3)
table(result_reassigned)
sorted_z0 <- sort(z_0)
result <- sort(result_reassigned)


# cluster confusion matrices to visualize z_hat output

conf_matrix <- as.data.frame(table(sorted_z0, result))

true_labels <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
names(true_labels) <- c(1, 2, 3, 4)

confusion_plot_sbm <- ggplot(conf_matrix, aes(x = factor(sorted_z0, levels = c(1, 2, 3, 4)), 
                                              y = factor(result, levels = rev(levels(factor(result)))), 
                                              fill = Freq)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = Freq), color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "lightblue") +
  scale_x_discrete(labels = true_labels) +
  labs(
    x = "True Clusters",
    y = "Inferred Clusters"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(family = "Times New Roman"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

ggsave("confusion_matrix_sbm.pdf", 
       plot = confusion_plot_sbm, 
       width = 7, 
       height = 5, 
       units = "in")

# ============= Posterior checking and model comparison ==============

# WAIC 

set.seed(1)
V <- dim(A)[1]
burn_in <- 1000
N_iter <- 10000
a <- b <- 1
LL <- matrix(nrow=V*(V-1)/2,ncol=9000)


Z_SBM_WAIC <- matrix(0, (N_iter - burn_in), V)  
for(i in (burn_in + 1):N_iter){  
  Z_SBM_WAIC[i - burn_in, ] <- fit_sbm6$z_samples[i,]
}

Z_SBM_WAIC <- t(Z_SBM_WAIC)

for (t in 1:dim(Z_SBM_WAIC)[2]){
  LL[,t]<-sampleLL(as.factor(Z_SBM_WAIC[,t]),A,a,b)
  if (t%%10000 == 0){print(paste("Iteration:", t))}
}
WAIC(LL)$WAIC

quantile(apply(t(Z_SBM_WAIC)[,(burn_in+1):N_iter],2,max))[c(2:4)]


# when we misspecify K, i.e. K = 6
result_reassigned <- recode(fit_sbm$z_hat, `1` = 4, `3` = 1, `5` = 3)
table(result_reassigned, z_0)
library(clue)  # For solving label permutations
# Get all possible cluster labels (ensuring full range exists)
all_clusters <- sort(union(unique(fit_sbm$z_hat), unique(z_0)))
# Convert to factor, forcing all clusters to exist
true_labels <- factor(z_0, levels = all_clusters)
predicted_labels <- factor(fit_sbm$z_hat, levels = all_clusters)
# Create confusion matrix (ensuring all clusters exist)
confusion_matrix <- table(predicted_labels, true_labels)
# Solve the assignment problem
best_mapping <- solve_LSAP(confusion_matrix, maximum = TRUE)
# Reassign cluster labels
mapped_labels <- best_mapping[predicted_labels]
# Compute misclassification rate
misclassified <- sum(mapped_labels != as.numeric(true_labels))
misclassification_rate <- misclassified / length(true_labels)
# Print result
print(paste("Misclassification Rate:", round(misclassification_rate, 3)))

# save the output for all four models
save(fit1, Z_DP, Z_GN, fit_sbm, file = "sim0_results.RData")

