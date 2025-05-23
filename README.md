# SESCR
Code and Data for *A Spatial Capture-Recapture Model with Hawkes-inspired Detection Rates to account for Dependencies between Traps*, Biometrics, 2025*

## Simulations
Sim-SESCR.R and Sim-OU.R simulate data for Sections 3.1 and 3.2 respectively. Sim-SCR.R simulates data for Web Appendix C. 
SimulateFunctions.R contains functions used by MLE.R and Bayesian.R
The data is saved as .csv files for further use. 

## Case Studies
The folder CS contains raw data for our three case studies found in Section 4 and Web Appendix F.
CaseStudies.R fits SCR models in the *secr* package and our model in both frameworks.
The user must change the variable *N* on line 5 to switch between the Martens (1), Sitka deer (2) and Belarus lynx (3)

## Fitting
MLE.R fits our SESCR model by MLE/semi-complete likelihood as described in Section 2.5.1. 
Bayesian.R fits our SESCR model in a Bayesian data augmentation framework as described in Section 2.5.2.
The user can modify lines 9-12 to switch between fitting the output of the three simulation scripts. 
MLEFunction.cpp is the *TMB* template for MLE, which is compiled if necessary by the two fitting files.
Bayesian model fitting is done in *NIMBLE*.
FitFunctions.R contains functions for processing and preparing data and the mask, as well as priors and the likelihood for data augmentation.


