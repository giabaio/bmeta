\name{acf.plot}
\alias{acf.plot}
\title{
Autocorrelation function plot
}
\description{
Function to create autocorrelation function plot to assess convergence
}
\usage{
acf.plot(x,node,title="Autocorrelation function")
}
\arguments{
  \item{x}{
a \code{bmeta} object with results of the model
}
  \item{node}{
variable to be displayed on the plot
}
  \item{title}{
title of the plot, if specified 
}
}
\value{
A plot showing the autocorrelation for the selected node
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

### generate output from bmeta 
x <- bmeta(data=data.list,outcome="bin",model="std.dt",type="ran")

### generate autocorrelation function plot 
acf.plot(x,"alpha[1]")

### generate autocorrelation function plot and specify the title 
acf.plot(x,"alpha[1]",title="Autocorrelation plot")
}
\keyword{Autocorrelation}
\keyword{MCMC}
