# SESCR
`R` Code and Data for *A Spatial Capture-Recapture Model with Hawkes-inspired Detection Rates to account for Dependencies between Traps*, 2025*

## Simulations
*Sim-SESCR.R* and *Sim-OU.R* simulate data for Sections 3.1 and 3.2 respectively.<br/> 
*Sim-SCR.R* simulates data for Web Appendix C.<br/>
*SimulateFunctions.R* contains functions used by these three files. 
The data is saved as .csv files for later use. 

## Case Studies
The folder CS contains raw data for our two case studies found in Section 4 and Web Appendix E.
*CaseStudies.R* fits SCR models in the `secr` package and our model in both frameworks.
The user must change the variable *N* on line 5 to switch between the Martens (1) and Sitka deer (2).

## Fitting
*Frequentist.R* fits our SESCR model by semi-complete likelihood as described in Section 2.5.1.<br/> 
*Bayesian.R* fits our SESCR model in a Bayesian data augmentation framework as described in Section 2.5.2.
The user can modify lines 9-12 to switch between fitting the output of the three simulation scripts.
Bayesian model fitting is done in *NIMBLE*.<br/>
*SCLFunction.cpp* is the *TMB* template with the semi-complete likelihood function, which is compiled if necessary by the two fitting files.<br/>
*FitFunctions.R* contains functions for processing and preparing data and the mask, as well as priors and the likelihood for data augmentation.


