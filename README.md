# BaySSpaTexplicitSurv

**BaySSpaTexplicitSurv** is an R package for simulating and modeling spatial survival data. The package focuses on **right-censored survival data** with spatial correlation and implements **Metropolis-Hastings MCMC** for parameter estimation in STSU Weibull models.

---

## Features

- Simulate spatial survival data with covariates.  
- Handle right-censored data with spatial structure.  
- Implement Metropolis-Hastings MCMC for parameter estimation.  
- Calculate model evaluation metrics: DIC, BIC, MSE.  
- Provide MCMC convergence diagnostics: Gelman-Rubin, trace plots, cumuplot, correlation plots.

---

## Dependencies

The package requires the following R packages:

```r
install.packages(c(
  "MHadaptive", "MLmetrics", "metropolis", "coda", "BayesianTools",
  "geoR", "purrr", "DWreg", "DiscreteWeibull", "invgamma",
  "mvtnorm", "FAmle", "splines", "R2WinBUGS", "mcmcplots",
  "MCMCpack", "psych"
))
