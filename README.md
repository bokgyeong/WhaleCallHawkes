# Supplemental Code for "Analyzing Whale Calls through Novel Multivariate Hawkes Processes"
Authors: Bokgyeong Kang, Erin M. Schliep, Alan E. Gelfand, Tina M. Yack, Christopher W. Clark, and Robert S. Schick

Detailed instructions for implementing the proposed multivariate Hawkes process models, as well as for generating the tables and figures presented in the paper, are provided.

## Required packages
The code has been tested with R version 4.4.0, "Puppy Cup."  The following R packages must be installed before the code will run successfully.

`Rcpp`, `RcppArmadillo`, `foreach`, `batchmeans`, `tidyverse`, `readr`, `coda`, `egg`, `grid`, `xtable`, `sf`, `tigris`, `suncalc`, `cowplot`

## Simulation study

### Step 1: Simulate data
- `dataNHPP.R`: Generate a dataset from the model (i) NHPP
- `dataNHPPSE.R`: Generate a dataset from the model (ii) NHPP+E
- `dataLGCP.R`: Generate a dataset from the model (iii) NHPP+GP
- `dataLGCPSE.R`: Generate a dataset from the model (iv) NHPP+GP+E
- The simulated datasets are saved in the directory `/sim/data`

### Step 2: Fit the multivariate Hawkes process models to the simulated datasets 
- `fitNHPP.R`: Fit the model (i) NHPP to the simulated datasets
- `fitNHPPSE.R`: Fit the model (ii) NHPP+E to the simulated datasets
- `fitLGCP.R`: Fit the model (iii) NHPP+GP to the simulated datasets
- `fitLGCPSE.R`: Fit the model (iv) NHPP+GP+E to the simulated datasets
- The resulting posterior samples for model parameters are saved in the directory `/sim/fit`

### Step 3: Compare models via deviance information criterion (DIC) 

#### Compute loglikelihood
- `loglikNHPP.R`, `loglikNHPPSE.R`, `loglikLGCP.R`, `loglikLGCPSE.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for each model fitted to the simulated datasets
- The resulting posterior samples for the loglikelihood  are saved in the directory `/sim/loglik`
- We suggest determining the burn-in period by examining the trace plot of the loglikelihood chain

#### Obtain DIC
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `sumDIC.R`: Evaluate DIC for each model and create Table 2 included in the paper

### Step 4: Assess model adequacy via random time change theorem (RTCT)

#### Compute the expected number of calls occurring between two consecutive events
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `rtctNHPP.R`, `rtctNHPPSE.R`, `rtctLGCP.R`, `rtctLGCPSE.R`: Evaluate $d^{\ast}_{b,i}$ for each model fitted to the simulated datasets
- The results are saved in the directory `/sim/rtct`

#### Obtain Q-Q plot and mean squared difference (MSD) for RTCT
- `sumRTCT.R`: Calculate the posterior mean estimates of the order statistics $\{d_{(i)}^{*}\}$, along with their associated uncertainties and generate Figure 4 and Table 1 as presented in the paper

### Step 5: Perform inference using a compensator

#### Evaluate the expected number of calls
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `numNHPP.R`, `numNHPPSE.R`, `numLGCP.R`, `numLGCPSE.R`: Evaluate the expected total number of calls, expected number of contact calls, expected number of countercalls received at each hydrophone for each model fitted to the simulated datasets
- The results are saved in the directory `/sim/num`

#### Obtain the empirical posterior distributions for the expected number of total calls, contact calls, and countercalls
- `sumNum.R`: Evaluate the expected total number of calls, expected number of contact calls, expected number of countercalls received across hydrophones and create Figure 5 as outlined in the paper








