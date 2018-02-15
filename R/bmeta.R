#' Bayesian Meta Analysis/Meta-regression
#' 
#' Function to fit the Bayesian fixed- and random-effects meta-analytic models
#' with or without moderators. Models are designed to include non-informative
#' priors.
#' 
#' Specifying the data
#' 
#' The function can be used to evaluate odds ratios (or log odds ratios), mean
#' difference and incidence rate ratios (or log incidence rate ratios). Users
#' need to specify a list of data to be used in the function. For binary data,
#' events out of case and control arm and sample size of case and control arm
#' need to be listed. For continuous data, mean and standard errors of case and
#' control arm need to be listed if information is available. However, if only
#' mean difference and variance can be retrieved from each study, users need to
#' list mean difference and precision (inverse of variance). Notice that
#' information of all the studies need to be provided in the same format for
#' the function to work properly. For example, the function cannot work if some
#' of the studies provide mean and standard errors of the two arms while the
#' rest studies provide mean difference and variance. For count data, total
#' number of events in the follow-up period of case and control arm, total
#' follow-up person-time in case and control arm should be listed.
#' 
#' If additional impacts of a variable or more than one variable are observed
#' (when meta-regression is expected to be used), users need to provide a
#' matrix with each column either containing a dummy variable or a continuous
#' variable. In case that categorical variables (i.e. ethnicity, age band) are
#' observed and included, users need to first choose a 'baseline' category as
#' reference and then create dummies for each of the rest categories.
#' 
#' Model selection
#' 
#' Apart from 'null' models which apply Bayesian methods to obtain
#' study-specific without pooling-effects, there are 22 models included in this
#' package for pooling study-specific estimates together and producing summary
#' estimate. The number of models designed for binary, continuous and count
#' data are 8, 8 and 6, respectively. The model selection process for binary
#' and count data requires users to specify not only whether meta-analysis or
#' meta-regression is wanted but also the priors to be used.
#' 
#' For binary data, normal and Student t-distribution priors for summary
#' estimates (on log scale) can be selected and it is indicated that Student
#' t-distribution has heavier tails and is therefore more robust to outliers.
#' The argument 'model' here includes 4 options --- \code{std.norm},
#' \code{std.dt}, \code{reg.norm}, \code{reg.dt}.
#' 
#' For continuous data, rather than specifying prior, users need to select
#' whether all studies included report mean and standard errors of two arms
#' separately or only mean difference and variance as discussed above in the
#' 'Specifying the data' section. The argument 'model' here includes 4
#' options--- \code{std.ta}, \code{std.mv}, \code{reg.ta}, \code{reg.mv}
#' ('model' ending with 'ta' represents 'two arms' and ending with 'mv'
#' represents 'mean and variance').
#' 
#' For count data, uniform and half-Cauchy distribution priors for the
#' variability of summary estimates (on log scale) can be selected. It is
#' suggested that half-Cauchy distribution has heavier tails and allows for
#' outliers and accommodates small variances closing to zero. It should be
#' noticed that there is no need to specify priors for fixed-effects models for
#' count data. The argument 'model' here includes 6 options --- \code{std},
#' \code{std.unif}, \code{std.hc}, \code{reg}, \code{reg.unif}, \code{reg.hc}.
#' 
#' In conjunction with the argument 'type'--- \code{fix} or \code{ran}, users
#' can select the specific model wanted for a certain type of data.
#' 
#' @aliases bmeta bmeta.default
#' @param data a data list containing information on observed data (including
#' moderators). See 'details'.
#' @param outcome type of outcome that needs to be specified. For binary,
#' continuous and count data, \code{bin}, \code{ctns} and \code{count} need to
#' be specified, respectively.
#' @param model type of model that needs to be specified. See 'details'.
#' @param type model type---either fixed-effects (\code{fix}) or random-effects
#' model(\code{ran}) needs to be specified.
#' @param n.iter number of iterations to be used in the simulation (default is
#' 10000)
#' @param n.burnin number of burn-in to be used in the simulation (default is
#' 5000)
#' @param n.samples The total number of MCMC simulations saved (including
#' thinning). Default at 1000
#' @param n.chains number of Markov chains to be used in the simulation
#' (default is 2)
#' @param model.file Name of the text file to which the model is saved
#' @param ...  additional (optional) arguments
#' @return \item{mod}{A \code{rjags} object with the results of the model}
#' \item{params}{a list of monitored parameters to be saved} \item{data}{the
#' original dataset} \item{inits}{a list with \code{n.chains} elements, with
#' each element itself being a list of starting values for the model or a
#' function generating initial values} \item{outcome}{selected type of outcome
#' (i.e. bin/ctns/count)} \item{type}{selected type of model (either
#' fixed-/random-effects)} \item{model}{selected model with specific priors}
#' \item{mod0}{independent model without pooling effects}
#' @author Tao Ding Gianluca Baio
#' @references Baio, G.(2012) Bayesian methods in health economics. Chapman
#' Hall, CRC.
#' 
#' Welton, N.J., Sutton, A.J., Cooper, N., Abrams, K.R. & Ades, A.E. (2012)
#' Evidence synthesis for decision making in healthcare. Chichester, UK: John
#' Wiley & Sons, Ltd.
#' @keywords Bayesian meta-analysis
#' @examples
#' 
#' ### Read and format the data (binary)
#' data = read.csv(url("http://www.statistica.it/gianluca/bmeta/Data-bin.csv"))
#' 
#' ### List data for binary outcome (for meta-analysis)
#' d1 <- data.list <- list(y0=data$y0,y1=data$y1,n0=data$n0,n1=data$n1) 
#' 
#' ### List data for binary outcome when there is a covariate (for meta-regression)
#' d1 <- data.list <- list(y0=data$y0,y1=data$y1,n0=data$n0,n1=data$n1,X=cbind(data$X0)) 
#' 
#' ### Select fixed-effects meta-analysis with normal prior for binary data 
#' m1 <- bmeta(d1, outcome="bin", model="std.norm", type="fix",n.iter=100)
#' 
#' ### Select random-effects meta-regression with t-distribution prior for binary
#' ### data
#' m2 <- bmeta(data.list, outcome="bin", model="reg.dt", type="ran",n.iter=100)
#' 
#' 
#' 
#' ### Read and format the data (continuous)
#' data = read.csv(url("http://www.statistica.it/gianluca/bmeta/Data-ctns.csv"))
#' 
#' ### List data for continuous outcome for studies reporting two arms separately
#' ### (for meta-analysis)
#' d1 <- data.list <- list(y0=data$y0,y1=data$y1,se0=data$se0,se1=data$se1) 
#' 
#' ### List data for continuous outcome for studies reporting mean difference and 
#' ### variance with a covariate (for meta-regression)
#' d2 <- data.list2 <- list(y=data$y,prec=data$prec,X=cbind(data$X0))
#' 
#' ### Select fixed-effects meta-analysis with studies reporting information of 
#' ### both arm for continuous data 
#' m1 <- bmeta(data.list, outcome="ctns", model="std.ta", type="fix",n.iter=100)
#' 
#' ### Select random-effects meta-regression with studies reporting mean difference and 
#' ### variance only for continuous data
#' m2 <- bmeta(data.list2, outcome="ctns", model="reg.mv", type="ran",n.iter=100)
#' 
#' 
#' 
#' ### Read and format the data (count)
#' data = read.csv(url("http://www.statistica.it/gianluca/bmeta/Data-count.csv"))  
#' 
#' ### List data for count outcome (for meta-analysis)
#' d1 <- data.list <- list(y0=data$y0,y1=data$y1,p0=data[,6],p1=data[,10])
#' 
#' ### List data for count outcome when there is a covariate (for meta-regression)
#' d2 <- data.list <- list(y0=data$y0,y1=data$y1,p0=data[,6],p1=data[,10],X=cbind(data$X0)) 
#' 
#' ### Select fixed-effects meta-analysis for count data
#' m1 <- bmeta(d1, outcome="count", model="std", type="fix",n.iter=100)
#' 
#' ### Select random-effects meta-analysis with half-Cauchy prior for count data
#' m2 <- bmeta(d1, outcome="count", model="std.hc", type="ran",n.iter=100)
#' 
#' ### Select random-effects meta-regression with uniform prior for count data
#' m3 <- bmeta(d2, outcome="count", model="reg.unif", type="ran",n.iter=100)
#' 
#' 
bmeta<- function(data,outcome=c("bin","ctns","count"),model=c("std.norm","std.dt","reg.norm","reg.dt",
                                                              "std.ta","std.mv","reg.ta","reg.mv",
                                                              "std","std.unif","std.hc","reg","reg.unif","reg.hc"),
                 type=c("fix","ran"),n.iter=10000,n.burnin=5000,n.samples=1000,n.chains=2,model.file="model.txt",...)  UseMethod("bmeta")

