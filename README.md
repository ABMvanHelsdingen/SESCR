# SESCR
Code and Data for *A Hawkes-type Model extending Spatial Capture-Recapture*, Biometrics, 2024*

## Simulations
Sim-SESCR.R and Sim-OU.R simulate data for Sections 3.1 and 3.2 respectively. Sim-SCR.R simulates data for Appendix C. 
SimulateFunctions.R contains functions used by all scripts. 
The data is saved as .csv files for further use. 

## Fitting
MLE.R fits our SESCR model by MLE as described in Section 2.5.
Bayesian.R fits our SESCR model in a Bayesian framework as described in Section 2.6.
The user can modify lines 9-11 to switch between fitting the output of the three simulation scripts. 
MLEFunction.cpp is the TMB template, which is compiled if necessary by the two fitting files.
FitFunctions.R contains functions for processing and preparing data and the mask, as well as priors and the likelihood for NIMBLE. 


