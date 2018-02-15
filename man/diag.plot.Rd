\name{diag.plot}
\alias{diag.plot}
\title{
Diagnostic plot to examine model fit
}
\description{
Function to produce plot based on different diagnostic statistics 
}
\usage{
diag.plot(x,diag="Rhat",text=TRUE)
}
\arguments{
  \item{x}{
a \code{bmeta} object with results of the model
}
  \item{diag}{
diagnostic statistics to be used---either the Gelman-Rubin statistic (Rhat) by default or
effective sample size (n.eff)
}
  \item{text}{
should the name of the variables be shown in the graph (TRUE) or not (FALSE)
}

}
\value{
A plot showing the relevant diagnostic stats for each node in the model
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

### generate output using bmeta 
x <- bmeta(data=data.list,outcome="bin",model="std.norm",type="fix")

### run the diagnostic plot to examine the Gelman-Rubin statistic 
diag.plot(x)

### run the diagnostic plot to examine the effective sample size 
diag.plot(x,diag="n.eff")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MCMC}
\keyword{Diagnostics}% __ONLY ONE__ keyword per line
