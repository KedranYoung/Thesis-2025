bayesian_sbm_gibbs <- function(A, K, alpha_prior = rep(1, K), iter = 20000, burnin = 2000) {
  library(gtools)
  
  n <- nrow(A)
  z <- sample(1:K, n, replace = TRUE)  # Initialize cluster assignments
  
  # Initialize Q matrix
  Q <- matrix(0, K, K)
  # within-group probs
  for (k in 1:K) {
    Q[k, k] <- rbeta(1, 8, 2)
  }
  
  # between-group probs 
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      Q[i, j] <- rbeta(1, 2, 8)
      Q[j, i] <- Q[i, j] 
    }
  }
  
  # Store results
  z_samples <- matrix(0, nrow = iter, ncol = n)
  Q_samples <- array(0, dim = c(iter, K, K))
  
  for (t in 1:iter) {
    
    # Update group proportions pi
    n_k <- table(factor(z, levels = 1:K))
    pi <- rdirichlet(1, alpha_prior + n_k)
    
    # Update Q matrix (edge probabilities)
    for (r in 1:K) {
      for (s in 1:K) {
        A_rs <- sum(A[z == r, z == s])
        n_rs <- sum(z == r) * sum(z == s) - ifelse(r == s, sum(z == r), 0)
        Q[r, s] <- rbeta(1, 1 + A_rs, 1 + n_rs - A_rs)
      }
    }
    
    # Update z_i for each node
    for (i in 1:n) {
      probs <- pi * sapply(1:K, function(l) {
        prod(Q[l, z[-i]]^A[i, -i] * (1 - Q[l, z[-i]])^(1 - A[i, -i]))
      })
      
      z[i] <- sample(1:K, 1, prob = probs / sum(probs))
    }
    
    # Store samples
    z_samples[t, ] <- z
    Q_samples[t, , ] <- Q
  }
  
  # Return posterior means after burn-in
  list(
    z_samples = z_samples,
    z_hat = apply(z_samples[(burnin+1):iter, ], 2, function(x) names(which.max(table(x)))),
    Q_hat = apply(Q_samples[(burnin+1):iter, , ], c(2, 3), mean)
  )
}

# posterior edge count
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