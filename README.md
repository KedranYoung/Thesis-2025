# Nonparametric Methods for Bayesian Community Detection in Complex Networks

This repository is meant to provide the R code used in my thesis to fulfill the requirements of the Masters in Statistics and Analytics degree at the University of Arkansas. I use (1) Legramanti's Extended Stochastic Block Model and (2) Geng's Mixture of Finite Mixtures Stochastic Block Model (MFM-SBM) for community detection in network analysis. (1) https://github.com/danieledurante/ESBM provides the collapsed Gibbs sampler for ESBM, and (2) https://github.com/gengjunxianjohn/MFM-SBM provides the collapsed Gibbs sampler for MFM-SBM. I provide the traditional SBM Gibbs sampler in `Traditional_SBM.R`. 

This thesis modifies the above samplers to both simulated data and three complex real-world networks to determine best model fit for community detection. 

## Traditional SBM

The traditional SBM has the hierarchical structure from which we implement the Gibbs sampler in the .R file `Traditional_SBM.R`. 

   $A_{ij} \mid z_i, z_j, Q \sim \text{Bernoulli}(Q_{z_i z_j})$ \\
    $z_i \mid \pi\sim \text{Categorical}(\pi_1, \dots, \pi_K)$ \\
    $\pi \sim \text{Dirichlet}(\alpha_1, \dots ,\alpha_K)$ \\
    $Q_{rs} \sim \text{Beta}(a, b)$ \\

## DP-SBM 

The Dirichlet Process prior can also be implemented as a nonparametric prior on the community partition structure. The model after which the Gibbs sampler is implemented is given, 

 $A_{ij} \mid z_i, z_j, Q \sim \text{Bernoulli}(Q_{z_i z_j})$ \\
 $z_i \mid P \sim P)$ \\
 $P \sim DP(\alpha P_0)$ \\
 $Q_{rs} \sim \text{Beta}(a, b)$ \\ 

 ## GN-SBM 

The Gnedin Process prior can also be implemented as a nonparametric prior on the community partition structure. The model after which the Gibbs sampler is implemented is given, 

 $A_{ij} \mid z_i, z_j, Q \sim \text{Bernoulli}(Q_{z_i z_j})$ \\
 $z_i \mid P \sim P)$ \\
 $P \sim GN(\gamma P_0)$ \\
 $Q_{rs} \sim \text{Beta}(a, b)$ \\

 ## MFM-SBM

The Mixture of Finite Mixtures SBM also implements a nonparametric prior on the community partition structure where a known prior is placed on the number of communities,


   $A_{ij} \mid - \sim \text{Bernoulli}(Q_{z_iz_j}), 1 \leq i < j \leq n$ \\    
   $\pi \mid K \sim \text{Dirichlet}(\alpha_1, \dots, \alpha_K)$ \\
   $p(z_i = k \mid \pi, K) = \pi_k, i = 1, \dots, n, \,\, k = 1, \dots K$  \\
   $Q_{rs} \sim \text{Beta}(a,b), r,s = 1, \dots, K$ \\
   $K \sim p(\cdot), p(\cdot) \text{ is a PMF}$ \\


 
