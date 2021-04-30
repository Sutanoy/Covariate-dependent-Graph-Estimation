# Covariate-dependent-Graph-Estimation

There are two demo codes in this repository which showcases the practical performance of the graph estimator proposed in our paper "An approximate Bayesian approach to covariate dependent graphical modeling". The code discrete_covariate_demo.R considers the case of discrete covariates. The code cont_covariate_attempt2.R considers the toy example presented in the paper with continuous covariates. The cov_vsvb.R and ELBO_calculator are functions called by the demo codes. Specifically, the function cov_vsvb updates the variational parameters and returns the converged estimates for a single graph.

Overview of discrete_covariate_demo.R:
=====================================
One can simply run the demo as is to get some demo examples and some visual results through a heatmap and histograms.

## Data generation
In this file, it is assumed that there are 2 discrete covariate levels. The data are generated from two different covariance matrices as an example, controlled by a \lambda parameter. Depending on whether \lambda_1= \lambda_2, we have the covariate independent model or the covariate dependent model. 

#1. Covariate independent model
`\lambda_1 = \lambda_2 = [15,15,15,15, 0,... 0]`

#2. Covariate dependent model
`\lambda_1 = [15,15,15,15,0,...0]`
`\lambda_2 = [0,...0,15,15,15,15]`

The precision matrix `\Omega_1 = \lambda_1\lambda_1^T + 10 I`, where `I` is the identity matrix. Similarly, `\Omega_2= \lambda_2\lambda_2^T + 10 I`.
Let `Sigma_1= \Omega_1 ^{-1}`, and `Sigma_2 = \Omega_2^{-1}`.
We generate `n/2` samples from `Normal (0, \Sigma_1)` and `n/2` samples from `Normal (0, \Sigma_2)` to form our dataset.

## Generating the covariates
We generate an `n*p` covariate matrix with entries = `-0.1` for each of the p variables for the population with precision matrix `\Omega_1` and entries = `0.1` for each of the p variables for the population with precision matrix `Omega_2`. Thus for each variable, the covariate value for an individual is univariate. As an example, the covariate attached to the FOXC2 protein expression of patient 1 is the univariate FOXC2 RNA expression for the same patient. If instead we used both the RNA expression and CNV expression for FOXC2 gene for the same patient, we would have a two-dimensional covariate attached to the data.

## Overview of the algorithm
1. Fix a variable  `j` as response, and the remaining `p-1` variables as predictor. 
2. From the covariate matrix, define an `n*n` weight matrix where the `i`th row describes the weight vector associated with the n subjects relative to subject `i`. The weights for this model are chosen with an ad-hoc bandwidth value of 0.1. Technically, one can perform a density estimation on the covariate space, but since its basically discrete, we choose a small value of 0.1 
3. Choose the hyperparameter values \pi and \sigma^2 following Carbonetto Stephens and sigma_{\beta}^2 over a grid.
4. Call the cov_vsvb function to update the following variational parameters : The `n* p-1` matrices alpha, mu and S_sq where the `i`th row corresponds to the inclusion probability of the `p-1` predictor variables, mean and standard deviation for the `i`th subject.  
5. Loop over the `p` variables as response to get the `p*(p-1)` matrices corresponding to each of the n subjects in the study.
6. Assume that the diagonal elements in the inclusion probability matrices for each individual is `0`, and apply the post processing `\alpha_{jk}*=(\alpha_{jk} + \alpha_{kj})/2` to symmetrize the matrix.
7. Set the dependence graph to be I{\alpha_{jk}*>0.5}.


## Remarks 
The code calls the cov_vsvb function, which updates the variational parameters and returns the final values of the estimates corresponding to a single graph which corresponds to a fixed individual in the study. By going through a loop, one can calculate the graph estimates corresponding to every individual in the study, and can also be parallelized since the updates are independent. However, the parallelization is yet to be implemented.
The cov_vsvb function itself calls the ELBO_calculator function, which calculates the ELBO corresponding to the current values of the variational parameters for a specific graph corresponding to a single individual. Note the contribution of every individual in the study to the ELBO of the parameters corresponding to a single individual, facilitating the borrowing of information.

#####################################
Overview of cont_covariate_demo.R:

In this file, instead of discrete covariate values, we have three clusters of covariate values, and the individuals in the study have covariate values belonging to one of the three clusters.
