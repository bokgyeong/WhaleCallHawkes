# Supplemental Code for "Analyzing Whale Calls through Novel Multivariate Hawkes Processes"
Authors: Bokgyeong Kang, Erin M. Schliep, Alan E. Gelfand, Tina M. Yack, Christopher W. Clark, and Robert S. Schick

Detailed instructions for implementing the proposed multivariate Hawkes process models, as well as for generating the tables and figures presented in the paper, are provided.

## Required packages
The code has been tested with R version 4.4.0, "Puppy Cup."  The following R packages must be installed before the code will run successfully.

- `Rcpp`
- `RcppArmadillo`
- `foreach`
- `batchmeans`
- `tidyverse`
- `readr`
- `coda`
- `egg`
- `grid`
- `xtable`
- `sf`
- `tigris`
- `suncalc`
- `cowplot`

## Simulate data
- `dataNHPP.R`: Generate a dataset from the model (i) NHPP
- `dataNHPPSE.R`: Generate a dataset from the model (ii) NHPP+SE
- `dataLGCP.R`: Generate a dataset from the model (iii) NHPP+GP
- `dataLGCPSE.R`: Generate a dataset from the model (iv) NHPP+GP+SE
- The simulated datasets are saved in the folder `/data`

## Fit the multivariate Hawkes process models to the simulated datasets 
- `fitNHPP.R`: Fit the model (i) NHPP to the simulated datasets
- `fitNHPPSE.R`: Fit the model (ii) NHPP+SE to the simulated datasets
- `fitLGCP.R`: Fit the model (iii) NHPP+GP to the simulated datasets
- `fitLGCPSE.R`: Fit the model (iv) NHPP+GP+SE to the simulated datasets
- The resulting posterior samples for model parameters are saved in the folder `/fit`

## Compute deviance information criterion (DIC) 

### Obtain posterior samples for loglikelihoods
`loglikNHPP.R`, `loglikNHPPSE.R`, `loglikLGCP.R`, `loglikLGCPSE.R`
- Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the models fitted to the simulated datasets
- The resulting posterior samples for the loglikelihood  are saved in the folder `/loglik`
- We suggest determining the burn-in period by examining the trace plot of the loglikelihood chain.

### Compute DIC


## Obtain posterior samples for random time change theorem (RTCT)
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain.
- `rtctNHPP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (i) NHPP fitted to the simulated datasets
- `rtctNHPPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (ii) NHPP+SE fitted to the simulated datasets
- `rtctLGCP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iii) NHPP+GP fitted to the simulated datasets
- `rtctLGCPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iv) NHPP+GP+SE fitted to the simulated datasets
- The results are saved in the folder `/rtct`
  
## Obtain posterior samples for the expected number of upcalls
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain.
- `numNHPP.R`: Evaluate the expected total number of upcalls, expected number of contact calls, expected number of countercalls for the model (i) NHPP fitted to the simulated datasets
- `numNHPPSE.R`: Evaluate the expected total number of upcalls, expected number of contact calls, expected number of countercalls for the model (ii) NHPP+SE fitted to the simulated datasets
- `numLGCP.R`: Evaluate the expected total number of upcalls, expected number of contact calls, expected number of countercalls for the model (iii) NHPP+GP fitted to the simulated datasets
- `numLGCPSE.R`: Evaluate the expected total number of upcalls, expected number of contact calls, expected number of countercalls for the model (iv) NHPP+GP+SE fitted to the simulated datasets
- The results are saved in the folder `/num`
