#' Autocorrelation function plot
#' 
#' Function to create autocorrelation function plot to assess convergence
#' 
#' 
#' @param x a \code{bmeta} object with results of the model
#' @param node variable to be displayed on the plot
#' @param title title of the plot, if specified
#' @return A plot showing the autocorrelation for the selected node
#' @author Tao Ding Gianluca Baio
#' @keywords Autocorrelation MCMC
#' @examples
#' 
#' ### Read and format the data (binary)
#' data = read.csv(url("http://www.statistica.it/gianluca/bmeta/Data-bin.csv"))
#' 
#' ### List data for binary outcome 
#' data.list <- list(y0=data$y0,y1=data$y1,n0=data$n0,n1=data$n1) 
#' 
#' ### generate output from bmeta 
#' x <- bmeta(data=data.list,outcome="bin",model="std.dt",type="ran")
#' 
#' ### generate autocorrelation function plot 
#' acf.plot(x,"alpha[1]")
#' 
#' ### generate autocorrelation function plot and specify the title 
#' acf.plot(x,"alpha[1]",title="Autocorrelation plot")
#' 
acf.plot <- function(x,node,title="Autocorrelation function") {
  requireNamespace("R2jags",quietly=TRUE)
  ##cmd <- ifelse(class(m)=="rjags",mdl <- m$BUGSoutput,mdl <- m)
  mdl <- x$mod$BUGSoutput
  ind2 <- 1:dim(mdl$summary)[1]
  ind <- which(colnames(mdl$sims.matrix)==node)
  if (var(mdl$sims.matrix[,ind])==0) {return(invisible())}
  acf(mdl$sims.matrix[,ind],ylab="Autocorrelation",main=title)
}


