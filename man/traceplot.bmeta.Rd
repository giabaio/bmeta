\name{traceplot.bmeta}
\alias{traceplot.bmeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Traceplot to assess convergence
}
\description{
Function to display a plot of iteration vs. sample values for each variable in the chain
}
\usage{
traceplot.bmeta(x,node,title="",lab="")
}
\arguments{
  \item{x}{
a \code{bmeta} object with results of the model
}
  \item{node}{
variable to be displayed on the traceplot
}
  \item{title}{
title of the plot if specified
}
  \item{lab}{
name of the variable to be displayed on the traceplot
}

}
\author{
Tao Ding
Gianluca Baio
}
\examples{
### Read and format the data (binary)
data = read.csv(url("http://www.statistica.it/gianluca/bmeta/Data-bin.csv"))

### List data for binary outcome 
data.list <- list(y0=data$y0,y1=data$y1,n0=data$n0,n1=data$n1) 

### Select random-effects meta-analysis with t-distribution prior for binary
### data
x <- bmeta(data.list, outcome="bin", model="std.dt", type="ran")


### using output from bmeta to produce traceplot for a specific node
traceplot.bmeta(x,"mu")

### using output from bmeta to produce traceplot and specify the node used
traceplot.bmeta(x,"mu",lab="mu")
}
\keyword{MCMC}
\keyword{Diagnostics}% __ONLY ONE__ keyword per line
