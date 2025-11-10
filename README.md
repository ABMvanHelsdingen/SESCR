# SESCR
`R` Code and Data for *A Spatial Capture-Recapture Model with Hawkes-inspired Detection Rates to account for Animal Movement*, 2025*

## Case Studies
The folder CS contains raw data for our two case studies found in Section 4 and Web Appendix E.
*CaseStudies.R* fits SCR models in the `secr` package and our model in both frameworks.
The user must change the variable *N* on line 15 to switch between the marten (1) and the leopard (3) case studies. Case Study 2 is not used by the paper. 

## Simulations
*Sim-SESCR.R* and *Sim-OU.R* simulate data for Sections 3.1 and 3.2 respectively.<br/> 
*Sim-SCR.R* simulates data for Web Appendix C.<br/>
Both files save output as .csv files for later use.
*SimulateFunctions.R* contains functions used by these three files (see top of file for list)

## Simulation Results
Raw results of the simulation studies in .csv format can be found in the *Simulation Results* folder, along with summaries. These summaries show the mean estimates of *N* and *$\sigma$* grouped by true parameter values, the proportion of simulations that SESCR gives a better estimate for each, and mean squared errors.<br/> 
The start of the file names (SESCR, SCR or OUSCR) refers to the simulation mechanism. DA refers to Bayesian Data Augmentation, and F to frequentist MLEs. <br/>

## Figures
Figures of the camera layout and survey region masks used for the simulations and two case studies are included alongside a video representation of our model. In the video, the red dot indicates the location of a detection at time 0, and the black dot is the individual's activity center.<br/>


## Fitting
*Frequentist.R* fits our SESCR model by semi-complete likelihood as described in Section 2.5.1.<br/> 
*Bayesian.R* fits our SESCR model in a Bayesian data augmentation framework as described in Section 2.5.2.
The user can modify lines 4-6 to select the data to which the model is fitted.
Bayesian model fitting is done in `NIMBLE`.<br/>
*SCLFunction.cpp* is the `TMB` template with the semi-complete likelihood function, which is compiled if necessary by the two fitting files.<br/>
*FitFunctions.R* contains functions for processing and preparing data and the mask, as well as priors and the likelihood for data augmentation.
(see top of file for list)

