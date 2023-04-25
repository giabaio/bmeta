#' Traceplot to assess convergence
#' 
#' Function to display a plot of iteration vs. sample values for each variable
#' in the chain
#' 
#' 
#' @param x a \code{bmeta} object with results of the model
#' @param node variable to be displayed on the traceplot
#' @param title title of the plot if specified
#' @param lab name of the variable to be displayed on the traceplot
#' @author Tao Ding Gianluca Baio
#' @keywords MCMC Diagnostics
#' @examples
#' 
#' ### Read and format the data (binary)
#' data = read.csv(url("https://gianluca.statistica.it/software/bmeta/Data-bin.csv"))
#' 
#' ### List data for binary outcome 
#' data.list <- list(y0=data$y0,y1=data$y1,n0=data$n0,n1=data$n1) 
#' 
#' ### Select random-effects meta-analysis with t-distribution prior for binary
#' ### data
#' x <- bmeta(data.list, outcome="bin", model="std.dt", type="ran")
#' 
#' 
#' ### using output from bmeta to produce traceplot for a specific node
#' traceplot.bmeta(x,"mu")
#' 
#' ### using output from bmeta to produce traceplot and specify the node used
#' traceplot.bmeta(x,"mu",lab="mu")
#' 
traceplot.bmeta <- function(x,node,title="",lab=""){
  requireNamespace("R2jags",quietly=TRUE)
  ## node is a string with the name of the node to be plotted
  ## x is the name of the bmeta object containing the MCMC simulations
  xlab <- "Iteration"
  ## this way works with R2jags as well as with R2WinBUGS
  ###  cmd <- ifelse(class(model)=="rjags",mdl <- model$BUGSoutput,mdl <- model)
  mdl <- x$mod$BUGSoutput
  col <- colors()
  if (mdl$n.chains==1) {
    plot(mdl$sims.array[,1,node],t="l",col=col[which(col=="blue")],xlab=xlab,
         ylab=lab,main=title)
  }
  if (mdl$n.chains==2) {
    cols <- c("blue","red")
    plot(mdl$sims.array[,1,node],t="l",col=col[which(col==cols[1])],xlab=xlab,
         ylab=lab,main=title,ylim=range(mdl$sims.array[,1:2,node]))
    points(mdl$sims.array[,2,node],t="l",col=col[which(col==cols[2])])
  }
  if (mdl$n.chains>2) {
    cols <- c("blue","red","green","magenta","orange","brown","azure")
    plot(mdl$sims.array[,1,node],t="l",col=col[which(col==cols[1])],xlab=xlab,
         ylab=lab,main=title,
         ylim=range(mdl$sims.array[,1:mdl$n.chains,node]))
    for (i in 2:mdl$n.chains) {
      points(mdl$sims.array[,i,node],t="l",col=col[which(col==cols[i])])
    }
  }
}
