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
- All components are saved in the folder `/data`

## Fit the multivariate Hawkes process models to the simulated datasets
- `fitNHPP.R`: Fit the model (i) NHPP to the simulated datasets
- `fitNHPPSE.R`: Fit the model (ii) NHPP+SE to the simulated datasets
- `fitLGCP.R`: Fit the model (iii) NHPP+GP to the simulated datasets
- `fitLGCPSE.R`: Fit the model (iv) NHPP+GP+SE to the simulated datasets
