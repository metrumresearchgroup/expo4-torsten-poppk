# expo4-torsten-poppk

This repository contains R code demonstrating the use of `bbr.bayes` and Stan/Torsten in a typical Bayesian (Bayes) population pharmacokinetic (popPK) modeling and simulation (M&S) analysis, which includes model fitting, model evaluation, and model summarization. This demonstration uses the same processes and suite of tools used at Metrum Research Group (MetrumRG) to ensure traceable and reproducible pharmacometrics research; however, it is not meant to be a complete vignette on using all of the features of `bbr.bayes` or the other tools used in the workflow.

This repository contains examples of the various tasks associated with using `bbr.bayes` to fit models in Stan/Torsten.  It is expected that these
will evolve with time.

# Data

- Model estimation data set 
  - location: `data/derived/pk.csv`
  - data specification file: `data/derived/pk.yml`

# Scripts

The results presented in MeRGE Expo 4 were generated by running these scripts in the following order. The specific dependencies among the scripts are shown below. All of the scripts except `script/mcmc-diagnostics.Rmd` and `script/model-diagnostics.Rmd` require installation of Torsten.

- Torsten installation: `script/Torsten-installation.Rmd`
- Intial model creation and submission: `script/initial-model-submission.Rmd`
- MCMC convergence diagnostics: `script/mcmc-diagnostics.Rmd`
  - Requires `ppkexpo1` model and fitting results generated by `script/initial-model-submission.Rmd`.
- Model fitting diagnostics: `script/model-diagnostics.Rmd`
  - Requires `ppkexpo1` model and fitting results generated by `script/initial-model-submission.Rmd`.
- Model revision: `script/update-model.Rmd`
  - Requires `ppkexpo1` model generated by `script/initial-model-submission.Rmd`.
- Generated quantities model creation and submission: `script/generated-quantities.Rmd`
  - Requires `ppkexpo4` model generated by `script/update-model.Rmd`.
- Summarizing and comparing models: `script/modeling-summary.Rmd`
  - Requires `ppkexpo1` model and fitting results generated by `script/initial-model-submission.Rmd`.
  - Requires `ppkexpo2`, `ppkexpo3` and `ppkexpo4` model and fitting results generated by `script/update-model.Rmd`.
- Within-chain parallel computation: `script/parallel-computation.Rmd`
  - Requires `ppkexpo5` model and fitting results generated by `script/generated-quantities.Rmd`.
- Implementing shrinkage priors: `script/shrinkage-priors.Rmd`
  - Requires `ppkexpo4` model generated by `script/update-model.Rmd`.

# Helper functions
- MCMC diagnostics functions: `script/mcmc-diagnostic-functions.R`

# Metworx Version
- metworx-24-04.01.00

# R Version
- 4.3.1


Copied from internal repo at afd81ab5a82aae975ce687fe22f45c8d9f12dc93

