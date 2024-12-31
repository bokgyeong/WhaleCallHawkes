# Supplemental Code for "Analyzing Whale Calls through Novel Multivariate Hawkes Processes"
Authors: Bokgyeong Kang, Erin M. Schliep, Alan E. Gelfand, Tina M. Yack, Christopher W. Clark, and Robert S. Schick

Detailed instructions for implementing the proposed multivariate Hawkes process models, as well as for generating the tables and figures presented in the paper, are provided.

## Required packages
The code has been tested with R version 4.4.0, "Puppy Cup."  The following R packages must be installed before the code will run successfully.

`Rcpp`, `RcppArmadillo`, `foreach`, `batchmeans`, `tidyverse`, `readr`, `coda`, `egg`, `grid`, `xtable`, `sf`, `tigris`, `suncalc`, `cowplot`

## Simulation study

### Step 1. Simulate data
- `/sim/dataNHPP.R`: Generate a dataset from the model (i) NHPP
- `/sim/dataNHPPSE.R`: Generate a dataset from the model (ii) NHPP+E
- `/sim/dataLGCP.R`: Generate a dataset from the model (iii) NHPP+GP
- `/sim/dataLGCPSE.R`: Generate a dataset from the model (iv) NHPP+GP+E
- The simulated datasets are saved in the directory `/sim/data`

### Step 2. Fit the multivariate Hawkes process models to the simulated datasets 
- `/sim/fitNHPP.R`: Fit the model (i) NHPP to the simulated datasets
- `/sim/fitNHPPSE.R`: Fit the model (ii) NHPP+E to the simulated datasets
- `/sim/fitLGCP.R`: Fit the model (iii) NHPP+GP to the simulated datasets
- `/sim/fitLGCPSE.R`: Fit the model (iv) NHPP+GP+E to the simulated datasets
- The resulting posterior samples for model parameters are saved in the directory `/sim/fit`

### Step 3. Compare models via deviance information criterion (DIC) 

#### Compute loglikelihood
- `/sim/loglikNHPP.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (i) NHPP fitted to the simulated datasets
- `/sim/loglikNHPPSE.R`:Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (ii) NHPP+E fitted to the simulated datasets
- `/sim/loglikLGCP.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (iii) NHPP+GP fitted to the simulated datasets
- `/sim/loglikLGCPSE.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (iv) NHPP+GP+E fitted to the simulated datasets
- The resulting posterior samples for the loglikelihood  are saved in the directory `/sim/loglik`
- We suggest determining the burn-in period by examining the trace plot of the loglikelihood chain

#### Obtain DIC
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `/sim/sumDIC.R`: Evaluate DIC for each model and create Table 2 included in the paper

### Step 4. Assess model adequacy via random time change theorem (RTCT)

#### Compute the expected number of calls occurring between two consecutive events
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `/sim/rtctNHPP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (i) NHPP fitted to the simulated datasets
- `/sim/rtctNHPPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (ii) NHPP+E fitted to the simulated datasets
- `/sim/rtctLGCP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iii) NHPP+GP fitted to the simulated datasets
- `/sim/rtctLGCPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iv) NHPP+GP+E fitted to the simulated datasets
- The results are saved in the directory `/sim/rtct`

#### Summarize results
- `/sim/sumRTCT.R`: Calculate the posterior mean estimates of the order statistics $\{d_{(i)}^{*}\}$, along with their associated uncertainties. Generate Q-Q plots (Figure 4) and calculate the mean squared difference (Table 1) for models (i) to (iv) fitted to each of the simulated datasets

### Step 5. Perform inference using a compensator

#### Evaluate the expected number of calls
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `/sim/numNHPP.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (i) NHPP fitted to the simulated datasets
- `/sim/numNHPPSE.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (ii) NHPP+E fitted to the simulated datasets
- `/sim/numLGCP.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (iii) NHPP+GP fitted to the simulated datasets
- `/sim/numLGCPSE.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (iv) NHPP+GP+E fitted to the simulated datasets
- The results are saved in the directory `/sim/num`

#### Summarize results
- `/sim/sumNum.R`: Evaluate the expected counts for total calls, contact calls, and countercalls received across hydrophones, and obtain the corresponding empirical posterior distributions (Figure 5)


## Analysis of North Atlantic right whale upcall data

### Step 1. Fit the multivariate Hawkes process models to the dataset 
- `/real/fitNHPP.R`: Fit the model (i) NHPP to the dataset 
- `/real/fitNHPPSE.R`: Fit the model (ii) NHPP+E to the dataset 
- `/real/fitLGCP.R`: Fit the model (iii) NHPP+GP to the dataset 
- `/real/fitLGCPSE.R`: Fit the model (iv) NHPP+GP+E to the dataset 
- The resulting posterior samples for model parameters are saved in the directory `/real/fit`

### Step 2. Compare models via deviance information criterion (DIC) 

#### Compute loglikelihood
- `/real/loglikNHPP.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (i) NHPP fitted to the dataset
- `/real/loglikNHPPSE.R`:Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (ii) NHPP+E fitted to the dataset
- `/real/loglikLGCP.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (iii) NHPP+GP fitted to the dataset
- `/real/loglikLGCPSE.R`: Evaluate $\log L(\boldsymbol{\theta}_b \mid \mathcal{T})$ for the model (iv) NHPP+GP+E fitted to the dataset
- The resulting posterior samples for the loglikelihood  are saved in the directory `/real/loglik`
- We suggest determining the burn-in period by examining the trace plot of the loglikelihood chain

#### Obtain DIC
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `/real/sumDIC.R`: Evaluate DIC for each model and create Table 3 included in the paper

### Step 3. Assess model adequacy via random time change theorem (RTCT)

#### Compute the expected number of calls occurring between two consecutive events
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `/real/rtctNHPP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (i) NHPP fitted to the dataset
- `/real/rtctNHPPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (ii) NHPP+E fitted to the dataset
- `/real/rtctLGCP.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iii) NHPP+GP fitted to the dataset
- `/real/rtctLGCPSE.R`: Evaluate $d^{\ast}_{b,i}$ for the model (iv) NHPP+GP+E fitted to the dataset
- The results are saved in the directory `/real/rtct`

#### Summarize results
- `/real/sumRTCT.R`: Calculate the posterior mean estimates of the order statistics $\{d_{(i)}^{*}\}$, along with their associated uncertainties. Generate Q-Q plots (Figure 6) and calculate the mean squared difference (Table 3) for models (i) to (iv) fitted to the dataset

### Step 4. Perform inference using a compensator

#### Evaluate the expected number of calls
- We suggest using the burn-in period determined based on the trace plot of the loglikelihood chain
- `/real/numNHPP.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (i) NHPP fitted to the dataset
- `/real/numNHPPSE.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (ii) NHPP+E fitted to the dataset
- `/real/numLGCP.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (iii) NHPP+GP fitted to the dataset
- `/real/numLGCPSE.R`: Evaluate the expected total number of calls, the expected number of contact calls, the expected number of countercalls received at each hydrophone for the model (iv) NHPP+GP+E fitted to the dataset
- The results are saved in the directory `/real/num`

#### Summarize results
- `/real/sumNum.R`: Calculate the posterior mean estimates for the expected number of within-MARU countercalls and cross-MARU countercalls, and create Figure 8 as presented in the paper. Generate the joint posterior distribution of the expected number of contact calls and countercalls received at each MARU (Figure 9)






