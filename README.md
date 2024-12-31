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
- `loglikNHPP.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (i) NHPP fitted to the simulated datasets
- `loglikNHPPSE.R`:Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (ii) NHPP+E fitted to the simulated datasets
- `loglikLGCP.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (iii) NHPP+GP fitted to the simulated datasets
- `loglikLGCPSE.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (iv) NHPP+GP+E fitted to the simulated datasets
- The resulting posterior samples for the loglikelihood  are saved in the directory `/sim/loglik`
- We suggest determining the burn-in period by examining the trace plot of the loglikelihood chain

#### Obtain DIC
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `sumDIC.R`: Evaluate DIC for each model and create Table 2 included in the paper

### Step 4: Assess model adequacy via random time change theorem (RTCT)

#### Compute the expected number of calls occurring between two consecutive events
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `rtctNHPP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (i) NHPP fitted to the simulated datasets
- `rtctNHPPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (ii) NHPP+E fitted to the simulated datasets
- `rtctLGCP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iii) NHPP+GP fitted to the simulated datasets
- `rtctLGCPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iv) NHPP+GP+E fitted to the simulated datasets
- The results are saved in the directory `/sim/rtct`

#### Summarize results
- `sumRTCT.R`: Calculate the posterior mean estimates of the order statistics $\{d_{(i)}^{*}\}$, along with their associated uncertainties. Generate Q-Q plots (Figure 4) and calculate the mean squared difference (Table 1) for models (i) to (iv) fitted to each of the simulated datasets

### Step 5: Perform inference using a compensator

#### Evaluate the expected number of calls
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `numNHPP.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (i) NHPP fitted to the simulated datasets
- `numNHPPSE.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (ii) NHPP+E fitted to the simulated datasets
- `numLGCP.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (iii) NHPP+GP fitted to the simulated datasets
- `numLGCPSE.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (iv) NHPP+GP+E fitted to the simulated datasets
- The results are saved in the directory `/sim/num`

#### Summarize results
- `sumNum.R`: Evaluate the expected counts for total calls, contact calls, and countercalls received across hydrophones, and obtain the corresponding empirical posterior distributions (Figure 5)


## Analysis of North Atlantic right whale upcall data

### Step 1: Fit the multivariate Hawkes process models to the dataset 
- `fitNHPP.R`: Fit the model (i) NHPP to the dataset 
- `fitNHPPSE.R`: Fit the model (ii) NHPP+E to the dataset 
- `fitLGCP.R`: Fit the model (iii) NHPP+GP to the dataset 
- `fitLGCPSE.R`: Fit the model (iv) NHPP+GP+E to the dataset 
- The resulting posterior samples for model parameters are saved in the directory `/real/fit`

### Step 3: Compare models via deviance information criterion (DIC) 

#### Compute loglikelihood
- `loglikNHPP.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (i) NHPP fitted to the dataset
- `loglikNHPPSE.R`:Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (ii) NHPP+E fitted to the dataset
- `loglikLGCP.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (iii) NHPP+GP fitted to the dataset
- `loglikLGCPSE.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (iv) NHPP+GP+E fitted to the dataset
- The resulting posterior samples for the loglikelihood  are saved in the directory `/real/loglik`
- We suggest determining the burn-in period by examining the trace plot of the loglikelihood chain

#### Obtain DIC
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `sumDIC.R`: Evaluate DIC for each model and create Table 3 included in the paper

### Step 4: Assess model adequacy via random time change theorem (RTCT)

#### Compute the expected number of calls occurring between two consecutive events
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `rtctNHPP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (i) NHPP fitted to the dataset
- `rtctNHPPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (ii) NHPP+E fitted to the dataset
- `rtctLGCP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iii) NHPP+GP fitted to the dataset
- `rtctLGCPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iv) NHPP+GP+E fitted to the dataset
- The results are saved in the directory `/real/rtct`

#### Summarize results
- `sumRTCT.R`: Calculate the posterior mean estimates of the order statistics $\{d_{(i)}^{*}\}$, along with their associated uncertainties. Generate Q-Q plots (Figure 6) and calculate the mean squared difference (Table 3) for models (i) to (iv) fitted to the dataset

### Step 5: Perform inference using a compensator

#### Evaluate the expected number of calls
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `numNHPP.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (i) NHPP fitted to the dataset
- `numNHPPSE.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (ii) NHPP+E fitted to the dataset
- `numLGCP.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (iii) NHPP+GP fitted to the dataset
- `numLGCPSE.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (iv) NHPP+GP+E fitted to the dataset
- The results are saved in the directory `/real/num`

#### Summarize results
- `sumNum.R`:






