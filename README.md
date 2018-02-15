# bmeta [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/bmeta)](https://cran.r-project.org/package=bmeta)[![CRAN_Download_Badge](http://cranlogs.r-pkg.org/badges/bmeta)](https://cran.r-project.org/package=bmeta)[![CRAN_Download_Badge](http://cranlogs.r-pkg.org/badges/grand-total/bmeta?color=orange)](http://cranlogs.r-pkg.org/badges/grand-total/bmeta?color=orange)
## Bayesian meta-analysis & meta-regression in R

bmeta is a R package that provides a collection of functions for conducting meta-analyses and meta-regressions under a Bayesian context, using JAGS. The package includes functions for computing various effect size or outcome measures (e.g. odds ratios, mean difference and incidence rate ratio) for different types of data (e.g. binary, continuous and count, respectively), based on Markov Chain Monte Carlo (MCMC) simulations. Users are allowed to select fixed- and random-effects models with different prior distributions for the data and the relevant parameters. Meta-regression can be also performed if the effect of additional covariates are considered. Furthermore, the package provides functions for creating posterior distributions and forest plot to display main model output. Traceplots and some other diagnostic plots are also available for assessing model fit and performance.

bmeta works by allowing the user to specify the set of options in a standardised way in the R command terminal. Currently, bmeta implements 22 models; when the user has selected the required configuration (random vs fixed effects; choice of outcome and prior distributions; presence of covariates), bmeta writes a suitable JAGS file in the working directory. This is used to call the package R2Jags and perform the actual analysis, but can also be considered as a sort of "template" - the user can then modify to extend the modelling, modify the priors in a way that is not automatically done by bmeta, or use it for future references.



## Installation
There are two ways of installing `bmeta`. A "stable" version is packaged and available from CRAN. You can install using the following commands.
```R
install.packages("bmeta")
```
NB: `bmeta` requires that JAGS is also installed on your machine. You can find details of JAGS installation [here](mcmc-jags.sourceforge.net/).

The "development" version of `bmeta` is available from this GitHub repository - it will usually be updated more frequently and may be continuously tested. On Windows machines, you need to install a few dependencies, including [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first, e.g. by running
```R
pkgs <- c("R2jags","forestplot","Rtools","devtools")
repos <- c("https://cran.rstudio.com") 
install.packages(pkgs,repos=repos,dependencies = "Depends")
```
before installing the package using `devtools`:
```R
devtools::install_github("giabaio/bmeta")
```
Under Linux or MacOS, it is sufficient to install the package via `devtools`:
```R
install.packages("devtools")
devtools:install_github("giabaio/bmeta")
```
