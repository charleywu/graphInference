# Inference and search on graph-structured spaces

## Table of Contents
* [General info](general-info)
* [Data](#data)
* [Analyses](#analyses)

## General info
This repository contains all the code and (anonymized) data necessary for reproducing the Wu, Schulz, & Gershman (2020). The following text provides a summary of each file and it's function, while there are further comments throughout the code. 


## Data
* `exp1/` contains a csv of the (anonymized) experiment data `experiment1Full.csv` and `graphs.json` defining the graph structures used in the experiment
* `banditTask/` contains a csv of the (anonymized) experiment data `networkBandit.csv` and `network2.json` defining the graph structures used in the experiment


## Analyses
* `brmsModels/` contains Bayesian mixed effects models. Becuase of the file sizes potentially being larger than the maximum allowed limit in github, all `*.brm` files are tracked using `git-lfs`. Please refer to the [git-lfs](https://git-lfs.github.com/) manual to install lfs and pull these files
* `modelResults/` contains the individual cross-validated maximum likelihood estimates from Exp 1 (`Exp1/`) and Exp 2 (`networkBandit/`). It also contains `modelFit.csv` and `paramEstimates.csv` as the compiled dataframes describing the model results in Exp 2 for convenience.  In addition, `Exp1diffevidence.csv` and `Exp2diffevidence.csv` contain the log loss of each model for computing the protected exceedence probabilities, which is then saved as `Exp1PXP.csv` and `Exp2PXP.csv`. Lastly, the model-based analyses of the bonus from from Exp 2 are here for covenience as `BanditBonusRoundmodelDF.Rds`
* `plots/` contains plots for the paper, where the code will save each plot
* `rationalModels/` contains data used in the model simulations for Experiment 2. `parameters/` holds dataframes for the parameter estimates of each model, where the model simulation (`modelSimulations.R`) are saved as csv files
* `utilities.R` contains data preprocessing functions for each experiment and various vector operations that are used across multiple scripts
* `statisticalTests.R` contains code for performing t-tests and correlations, where the output is formatted for Latex and automatically converted to a set number of significant digits for consistency. Contains code from van Doorn et al., (2018) for computing the Bayes factor for Kendall's rank correlation [paper](https://amstat.tandfonline.com/doi/full/10.1080/00031305.2016.1264998)
* `exportImportGraph.R` contains code by Angelo Antonio Salatino for impoirting and exporting igraph objects as json: [original github repo](https://github.com/angelosalatino/graph-importer-R)
* `models.R` contains code defiing models used in both experiments

### Experiment 1: 

* `Exp1Behavior.R` contains all the behavioral analyses
* `Exp1ModelCV.R` contains code for model fitting 
* `Exp1ModelingResults.R` contains code for all model-based analyses

### Experiment 2: 
* `banditExpBehaviorPlots.R` contains all the behavioral analyses
* `banditModelComparisonCV.R` contains code for model fitting 
* `banditExpBehaviorPlots.R` contains code for model-based analyses
* `modelSimulations.R` contains code for computing the simulated learning curves
* `banditBonus.R` contains code used for analyzing the bonus round data

### PXP
* Code for computing the protected probability of exceedence (pxp) is taken from (original github repo)[https://github.com/sjgershm/mfit]
* `bms.m` contains the code for computing pxp, while `pxp.m` defines the input and output files





