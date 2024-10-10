# expo4-torsten-poppk

`bbr.bayes` is an extension to the `bbr` package for traceable and reproducible Bayesian modeling. Currently, `bbr.bayes` supports Stan models (powered by cmdstanr). Torsten extends Stan by providing a suite of functions to facilitate the specification of pharmacometrics models. This includes specification of complex event schedules often encountered when analyzing pharmacometric data.

This repository contains examples of the various tasks associated with using `bbr.bayes` to fit models in Stan.  It is expected that these
will evolve with time.

# Data

- Model estimation data set 
  - location : `data/derived/pk.csv`
  - data specification file : `data/derived/pk.yml`


# Scripts
Current file structure:
- Model management: `script/model-management-demo.Rmd`

The above model management script will be subdivided into multiple scripts corresponding to the online expo structure:
- Initial model creation and submission: `script/initial-model-submission.Rmd`
- Initial model diagnostics: `script/initial-model-diagnostics.Rmd`
- Iterative model development: `script/revised-models.Rmd`
- Standalone generated quantities: `script/generated-quantities.Rmd`
- Summarizing and comparing models: `script/modeling-summary.Rmd`
- Winthin-chain parallel computation: `script/parallel-computation.Rmd`


# Helper functions
- Helper functions for model management: `script/functions-model.R`
- MCMC diagnostics functions: `script/functions-mcmc-diagnostics.R`
- Posterior predictive checking functions:  `script/functions-ppc.R`

# Metworx Version
- metworx-22.09.1

# R Version
- 4.1.3


Copied from internal repo at d296a409bd7a045cc8c9d10257d37b2f88e9d662

