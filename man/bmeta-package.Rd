\name{bmeta-package}
\alias{bmeta-package}
\docType{package}
\title{
bmeta: A Bayesian Meta-Analysis Package for R
}
\description{
The bmeta package provides a collection of functions for conducting meta-analyses under 
Bayesian context in R. The package includes functions for computing various effect size 
or outcome measures (e.g. odds ratios, mean difference and incidence rate ratio) for 
different types of data based on MCMC simulations. Users are allowed to fit fixed- and 
random-effects models with different priors to the data. Meta-regression can be carried 
out if effects of additional covariates are observed. Furthermore, the package provides 
functions for creating posterior distribution plots and forest plot to display main model 
output. Traceplots and some other diagnostic plots are also available for assessing model 
fit and performance.
}
\details{
\tabular{ll}{
Package: \tab bmeta\cr
Type: \tab Package\cr
Version: \tab 0.1.3\cr
Date: \tab 2016-05-05\cr
License: \tab GPL2 \cr
LazyLoad: \tab yes\cr
}
Bayesian meta-analysis is becoming more frequently accepted as a statistical approach 
for evidence synthesis from multiple studies in health research. The Bayesian methods 
differ inherently from frequentist ones by assuming that model parameters are random 
quantities. Therefore, prior distributions for model parameters can be specified, which 
are normally based on external evidence. The bmeta function provides 22 models with 
commonly used priors for fitting different types of data (i.e. binary, continuous and 
count data). 
}
\author{
Tao Ding, Gianluca Baio

Maintainer: Gianluca Baio <gianluca@stats.ucl.ac.uk>
}
\references{
Alex J Sutton and Keith R Abrams.(2001).Bayesian methods in meta-analysis and evidence 
synthesis. Statistical Methods in Medical Research,10,277-303.

Welton,N.J., Sutton,A.J., Cooper,N., Abrams,K.R.& Ades,A.E.(2012) Evidence synthesis for 
decision making in healthcare. Chichester, UK: John Wiley & Sons, Ltd.
}
\keyword{Bayesian meta-analysis}
