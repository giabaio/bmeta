#' bmeta -- Bayesian Meta Analysis/Meta-regression
#' 
#' Function to fit the Bayesian fixed- and random-effects meta-analytic models
#' with or without moderators. Models are designed to include non-informative
#' priors.
#' 
#' Specifying the data
#' 
#' The function can be used to evaluate odds ratios (or log odds ratios), mean
#' difference and incidence rate ratios (or log incidence rate ratios). Users
#' need to specify a list of data to be used in the function.
#' 
#' For binary data, the meta-analysis models require users to provide a list of
#' values for events out of case arm (y1) and control arm (y0) and the sample
#' size of case arm (n1) and control arm (n0). If there are additional
#' covariates (when meta-regression should be applied), a matrix 'X' needs to
#' be used with each column either containing a dummy variable (either '0' or
#' '1') or a continuous variable. For binary and continuous covariates, the
#' matrix 'X' can be specified by the command X=data$covariate or
#' X=cbind(data$covariate1, data$covariate2). For categorical variables, users
#' need to firstly choose a 'baseline' category as reference and then create
#' dummies for each of the rest categories. For example, if there is a
#' covariate 'ethnicity' with 3 categories ---White, Black and Asian, and
#' suppose we use Asian as our baseline reference group, then dummies for the
#' rest 2 categories, White and Black, should be created, with 1 column
#' containing either '0' or '1' indicating whether a study use White population
#' or not and another one column in exactly the same format indicating whether
#' as study use Black population or not. For the next step, the matrix 'X'
#' should be specified by combining information from these two columns, i.e.
#' X=cbind(data$White,data$Black).
#' 
#' For continuous data, mean (y1,y0) and standard errors (se1,se0) of case and
#' control arm need to be listed, respectively, if information is available.
#' However, if only mean difference and variance can be retrieved from each
#' study, users need to list mean difference (y) and precision
#' (prec)---defining as the inverse of variance. Notice that information of all
#' the studies need to be provided in the same format for the function to work
#' properly. For example, the function cannot work if some of the studies
#' provide mean and standard errors of the two arms while the rest studies
#' provide mean difference and variance. The meta-regression models require
#' exactly the same format of covariates input as illustrated above for the
#' binary data.
#' 
#' For count data, total number of events in the follow-up period of case arm
#' (y1) and control arm (y0), total follow-up person-time in case arm (p1) and
#' control arm (p0) should be listed. The meta-regression models require
#' exactly the same format of input as illustrated above for the binary data.
#' 
#' Notice that for all meta-regression models, moderator effects can only be
#' incorporated into the model through the matrix 'X'.
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
#' The argument 'model' here includes 4 options
#' ---"std.norm","std.dt","reg.norm","reg.dt".
#' 
#' For continuous data, rather than specifying prior, users need to select
#' whether all studies included report mean and standard errors of two arms
#' separately or only mean difference and variance as discussed above in the
#' 'Specifying the data' section. The argument 'model' here includes 4
#' options---"std.ta","std.mv","reg.ta","reg.mv"('model' ending with 'ta'
#' represents 'two arms' and ending with 'mv' represents 'mean and variance').
#' 
#' For count data, uniform and half-Cauchy distribution priors for the
#' variability of summary estimates (on log scale) can be selected. It is
#' suggested that half-Cauchy distribution has heavier tails and allows for
#' outliers and accommodates small variances closing to zero. It should be
#' noticed that there is no need to specify priors for fixed-effects models for
#' count data. The argument 'model' here includes 6 options
#' ---"std","std.unif","std.hc","reg","reg.unif","reg.hc". The meta-analysis
#' models require users to provide a list of values for total number of events
#' in case arm (y1) and control arm (y0) in the follow-up period and the total
#' follow-up person-time in case arm (p1) and control arm (p0).
#' 
#' In conjunction with the argument 'type'---"fix" or "ran", users can select
#' the specific model wanted for a certain type of data.
#' 
#' @param data a data list containing information on observed data (including
#' moderators). See 'details'.
#' @param outcome type of outcome that needs to be specified. For binary,
#' continuous and count data, "bin","ctns" and "count" need to be specified,
#' respectively.
#' @param model type of model that needs to be specified. See 'details'.
#' @param type model type---either fixed-effects("fix") or random-effects
#' model("ran") needs to be specified.
#' @param n.iter number of iterations to be used in the simulation (default is
#' 10000)
#' @param n.burnin number of burn-in to be used in the simulation (default is
#' 5000)
#' @param n.samples number of total iterations to be saved
#' @param n.chains number of Markov chains to be used in the simulation
#' (default is 2)
#' @param model.file file containing the appropriate model selected by user
#' @param ...  additional (optional) arguments
#' @return \item{mod}{A rjags object with the results of the model}
#' \item{params}{a list of monitored parameters to be saved} \item{data}{the
#' original dataset} \item{inits}{a list with n.chains elements, with each
#' element itself being a list of starting values for the model or a function
#' generating initial values} \item{outcome}{selected type of outcome (i.e.
#' bin/ctns/count)} \item{type}{selected type of model (either
#' fixed-/random-effects)} \item{model}{selected model with specific priors}
#' \item{mod0}{independent model without pooling effects}
#' @note %% ~~further notes~~
#' @author Tao Ding Gianluca Baio
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Baio, G.(2012) Bayesian methods in health economics. Chapman
#' Hall, CRC. Welton, N.J., Sutton, A.J., Cooper, N., Abrams, K.R. & Ades, A.E.
#' (2012) Evidence synthesis for decision making in healthcare. Chichester, UK:
#' John Wiley & Sons, Ltd.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' ### Read and format the data (binary)
#' data = read.csv(url("https://gianluca.statistica.it/software/bmeta/Data-bin.csv"))
#' 
#' ### List data for binary outcome (for meta-analysis)
#' data.list <- list(y0=data$y0,y1=data$y1,n0=data$n0,n1=data$n1) 
#' 
#' ### List data for binary outcome when there is a covariate (for meta-regression)
#' data.list <- list(y0=data$y0,y1=data$y1,n0=data$n0,n1=data$n1,X=cbind(data$X0)) 
#' 
#' ### Select fixed-effects meta-analysis with normal prior for binary data 
#' bmeta(data.list, outcome="bin", model="std.norm", type="fix")
#' 
#' ### Select random-effects meta-regression with t-distribution prior for binary
#' ### data
#' bmeta(data.list, outcome="bin", model="reg.dt", type="ran")
#' 
#' 
#' 
#' ### Read and format the data (continuous)
#' data = read.csv(url("https://gianluca.statistica.it/software/bmeta/Data-ctns.csv"))
#' 
#' ### List data for continuous outcome for studies reporting two arms separately
#' ### (for meta-analysis)
#' data.list <- list(y0=data$y0,y1=data$y1,se0=data$se0,se1=data$se1) 
#' 
#' ### List data for continuous outcome for studies reporting mean difference and 
#' ### variance with a covariate (for meta-regression)
#' data.list <- list(y=data$y,prec=data$prec,X=cbind(data$X0))
#' 
#' ### Select fixed-effects meta-analysis with studies reporting information of 
#' ### both arm for continuous data 
#' bmeta(data.list, outcome="ctns", model="std.ta", type="fix")
#' 
#' ### Select random-effects meta-regression with studies reporting mean difference and 
#' ### variance only for continuous data
#' bmeta(data.list, outcome="ctns", model="reg.mv", type="ran")
#' 
#' 
#' 
#' ### Read and format the data (count)
#' data = read.csv(url("https://gianluca.statistica.it/software/bmeta/Data-count.csv"))
#' 
#' ### List data for count outcome (for meta-analysis)
#' data.list <- list(y0=data$y0,y1=data$y1,p0=data$p0,p1=data$p1)
#' 
#' ### List data for count outcome when there is a covariate (for meta-regression)
#' data.list <- list(y0=data$y0,y1=data$y1,p0=data$p0,p1=data$p1,X=cbind(data$X0)) 
#' 
#' ### Select fixed-effects meta-analysis for count data
#' bmeta(data.list, outcome="count", model="std", type="fix")
#' 
#' ### Select random-effects meta-analysis with half-Cauchy prior for count data
#' bmeta(data.list, outcome="count", model="std.hc", type="ran")
#' 
#' ### Select random-effects meta-regression with uniform prior for count data
#' bmeta(data.list, outcome="count", model="reg.unif", type="ran")
#' 
#' 
bmeta<- function(data,outcome=c("bin","ctns","count"),model=c("std.norm","std.dt","reg.norm","reg.dt",
                                                              "std.ta","std.mv","reg.ta","reg.mv",
                                                              "std","std.unif","std.hc","reg","reg.unif","reg.hc"),
                 type=c("fix","ran"),n.iter=10000,n.burnin=5000,n.samples=1000,n.chains=2,model.file="model.txt",...)  UseMethod("bmeta")

