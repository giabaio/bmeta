\name{funnel.plot}
\alias{funnel.plot}
\title{
Funnel plot to examine publication bias
}
\description{
Function to examine publication bias. For both fixed- and random-effects models, 
estimates from no-pooling effects model are used as study-specific estimates. For 
random-effects models, the corresponding fixed-effects models are implemented at 
background to obtain pooled estimate. For example, if users call \code{bmeta} to 
run random-effects meta-analysis with normal prior, fixed-effects meta-analysis 
with normal prior are implemented at background to obtain pooled estimate for 
graphing. In the absence of publication and heterogeneity, the scatter resembles 
a symmetrical funnel and the triangle area formed by connecting the centred summary 
estimate with its 2.5\% and 97.5\% quantiles on either side includes about 95\% of 
the studies if the fixed-effects model assumption holds (i.e. all the studies 
estimate the same effect).
}
\usage{
funnel.plot(x,xlab=NULL,ylab=NULL,title=NULL,xlim=NULL)
}
\arguments{
  \item{x}{
a \code{bmeta} object with results of the model
}
  \item{xlab}{
title of x-axis. If unspecified, the function sets an appropriate lable by default.
}
  \item{ylab}{
title of x-axis. If unspecified, the function sets an appropriate lable by default.
}
  \item{title}{
title of the plot if specified
}
\item{xlim}{
horozontal limits of the plot region. If unspecified, the function sets the horizontal 
plot limits to (-6,6).
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

### using output from bmeta to produce funnel plot 
funnel.plot(x)

### using output from bmeta and specify title of the plot 
funnel.plot(x,title="funnel plot")

### using output from bmeta and specify the limit of x-axis and title
funnel.plot(x,title="funnel plot",xlim=c(-2,1))
}
\keyword{Funnel plot}
