# Nonparametric Methods for Bayesian Community Detection in Complex Networks

This repository is meant to provide the R code used in my thesis to fulfill the requirements of the Masters in Statistics and Analytics degree at the University of Arkansas. I use (1) Legramanti's Extended Stochastic Block Model and (2) Geng's Mixture of Finite Mixtures Stochastic Block Model (MFM-SBM) for community detection in network analysis. (1) https://github.com/danieledurante/ESBM provides the collapsed Gibbs sampler for ESBM, and (2) https://github.com/gengjunxianjohn/MFM-SBM provides the collapsed Gibbs sampler for MFM-SBM. 

This thesis modifies the above samplers to both simulated data and three complex real-world networks to determine best model fit for community detection. 

## Traditional SBM

   $ A_{ij} \mid z_i, z_j, Q \sim \text{Bernoulli}(Q_{z_i z_j})$
    $z_i \mid \pi\sim \text{Categorical}(\pi_1, \dots, \pi_K)$
    $\pi \sim \text{Dirichlet}(\alpha_1, \dots ,\alpha_K)$
    $Q_{rs} \sim \text{Beta}(a, b)$

