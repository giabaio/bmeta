#' Diagnostic plot to examine model fit
#' 
#' Function to produce plot based on different diagnostic statistics
#' 
#' 
#' @param x a \code{bmeta} object with results of the model
#' @param diag diagnostic statistics to be used---either the Gelman-Rubin
#' statistic (Rhat) by default or effective sample size (n.eff)
#' @param text should the name of the variables be shown in the graph (TRUE) or
#' not (FALSE)
#' @return A plot showing the relevant diagnostic stats for each node in the
#' model
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
#' ### generate output using bmeta 
#' x <- bmeta(data=data.list,outcome="bin",model="std.norm",type="fix")
#' 
#' ### run the diagnostic plot to examine the Gelman-Rubin statistic 
#' diag.plot(x)
#' 
#' ### run the diagnostic plot to examine the effective sample size 
#' diag.plot(x,diag="n.eff")
#' 
#' 
diag.plot <- function(x,diag="Rhat",text=TRUE) {
  tab <- x$mod$BUGSoutput$summary
  
  if(diag=="Rhat") {
    plot(tab[,"Rhat"],xlab="Parameters",ylab="Gelman-Rubin statistic (Rhat)",col="white",
         main="Convergence diagnostics")
    if(text){text(1:dim(tab)[1],tab[,"Rhat"],rownames(tab),cex=.6)} else {
       points(1:dim(tab)[1],tab[,"Rhat"],cex=.5,pch=20)
    }
    abline(h=1.1,lwd=2,lty=2)
  }
  if (diag=="n.eff") {
    plot(tab[,"n.eff"],xlab="Parameters",ylab="Effective sample size",col="white",
         main="Convergence diagnostics")
    if(text){text(1:dim(tab)[1],tab[,"n.eff"],rownames(tab),cex=.6)} else {
       points(1:dim(tab)[1],tab[,"n.eff"],cex=.5,pch=20)
    }
    abline(h=x$mod$BUGSoutput$n.sims,lwd=2,lty=2)    
  }
}

