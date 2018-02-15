##### POSTERIOR PLOT ##### 




#' Posterior distribution plots for summary estimates and between-study
#' standard deviation (measurement of heterogeneity)
#' 
#' Function to create posterior distribution plots for summary estimates and
#' between-study standard deviation based on output from \code{bmeta}
#' 
#' 
#' @param x a \code{bmeta} object with results of the model
#' @param xlim horozontal limits of the plot region. If unspecified, the
#' function sets the horizontal plot limits to (-3,3) for binary and count data
#' and (-5,5) for continuous data.
#' @param xlab title for the x-axis
#' @param main title of the plot. If unspecified, the function sets an
#' appropriate title by default.
#' @param scale logical specifying whether summary estimates need to be
#' displayed on log ("log") or natural scale ("exp"). For continuous data,
#' summary estimates are always displayed on natural scale, therefore, users do
#' not need to specify this option.
#' @param heterogeneity logical specifying whether to resent posterior plot for
#' between-study standard deviation (TRUE) to examine heterogeneity of
#' different studies. If unspecified, FALSE by default.
#' @author Tao Ding Gianluca Baio
#' @references Anzures-Cabrera,J & Higgins,J.P.T.(2010) Graphical displays for
#' meta-analysis: An overview with suggestions for practice.Res Synth
#' Methods,1,66-80.
#' @keywords Bayesian meta-analysis
#' @examples
#' 
#' ### Read and format the data (binary)
#' data = read.csv(url("http://www.statistica.it/gianluca/bmeta/Data-bin.csv"))
#' 
#' ### List data for binary outcome 
#' data.list <- list(y0=data$y0,y1=data$y1,n0=data$n0,n1=data$n1)
#' 
#' ### Select random-effects meta-analysis with t-distribution prior for binary
#' ### data
#' x <- bmeta(data.list, outcome="bin", model="std.dt", type="ran")
#' 
#' ### using output from bmeta to produce posterior plot
#' posterior.plot(x)
#' 
#' ### using output from bmeta and specify the horizontal limits 
#' posterior.plot(x,xlim=c(-2,1))
#' 
#' ### using output from bmeta on natural scale and specify more options 
#' posterior.plot(x,xlim=c(-0.5,2.5),xlab="odds ratio",main="Posterior distribution
#' of pooled odds ratio", scale="exp")
#' 
#' ### examine heterogeneity by producing posterior plot for between-study standard
#' ### deviation
#' posterior.plot(x,heterogeneity=TRUE,xlim=c(0,3),xlab="between-study standard 
#' deviation")
#' 
#' 
posterior.plot<-function(x,xlim=NULL,xlab="",main="Posterior distribution Plot",scale="log",
                         heterogeneity=FALSE){
  requireNamespace("R2jags",quietly=TRUE)
  
  param.sims <- x$mod$BUGSoutput$sims.matrix
  
  if (heterogeneity==FALSE){
  if (x$outcome=="bin" & x$type=="fix"){ 
    if(scale=="log") {param <- param.sims[,"delta"]} 
    if(scale=="exp") {param <- param.sims[,"rho"]}
    
    if(is.null(xlim)){
      hist(param,30,xlim=c(-3,3),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    
    if(scale=="log") {abline(v=0,lwd=2)}
    if(scale=="exp") {abline(v=1,lwd=2)}
  }
  
  
  if (x$outcome=="bin" & x$type=="ran"){ 
    if(scale=="log") {param <- param.sims[,"mu"]} 
    if(scale=="exp") {param <- param.sims[,"rho"]}
    
    if(is.null(xlim)){
      hist(param,30,xlim=c(-3,3),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    if(scale=="log") {abline(v=0,lwd=2)}
    if(scale=="exp") {abline(v=1,lwd=2)} 
  }
  
  
  
  if (x$outcome=="ctns" & x$type=="fix" & (x$model=="std.ta"|x$model=="reg.ta")){
    param <- param.sims[,"delta"]	
    if(is.null(xlim)){  
      hist(param,30,xlim=c(-5,5),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    abline(v=0,lwd=2)
  }
  
  
  if (x$outcome=="ctns" & x$type=="ran" & (x$model=="std.ta"|x$model=="reg.ta")){
    param <- param.sims[,"mu"]
    if(is.null(xlim)){ 
      hist(param,30,xlim=c(-5,5),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
      
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    abline(v=0,lwd=2)
  }
  
  
  if (x$outcome=="ctns" & x$model=="std.mv"){ 
    param <- param.sims[,"mu"] 
    if(is.null(xlim)){ 
      hist(param,30,xlim=c(-5,5),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    abline(v=0,lwd=2)  
  }
  
  
  if (x$outcome=="ctns" & x$model=="reg.mv" & x$type=="fix"){  
    param <- param.sims[,"alpha"]
    if(is.null(xlim)){ 
      hist(param,30,xlim=c(-5,5),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    abline(v=0,lwd=2)  
  }
  
  
  if (x$outcome=="ctns" & x$model=="reg.mv" & x$type=="ran"){  
    param <- param.sims[,"mu"]
    if(is.null(xlim)){ 
      hist(param,30,xlim=c(-5,5),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    abline(v=0,lwd=2)  
  }
  
  
  if (x$outcome=="count" & x$type=="fix"){
    if(scale=="log") {param <- param.sims[,"delta"]} 
    if(scale=="exp") {param <- param.sims[,"IRR"]}
    if(is.null(xlim)){
      hist(param,30,xlim=c(-3,3),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    if(scale=="log") {abline(v=0,lwd=2)}
    if(scale=="exp") {abline(v=1,lwd=2)}
  }
  
  
  
  if (x$outcome=="count" & x$type=="ran"){
    if(scale=="log") {param <- param.sims[,"mu"]} 
    if(scale=="exp") {param <- param.sims[,"IRR"]}
    if(is.null(xlim)){
      hist(param,30,xlim=c(-3,3),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
    points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
    if(scale=="log") {abline(v=0,lwd=2)}
    if(scale=="exp") {abline(v=1,lwd=2)}
  }

  } else {
    
    if(x$type=="fix"){
      print("Heterogeneity test is designed for random-effects models.")
    } else {
    param<-param.sims[,"sigma"]
      if(is.null(xlim)){
      hist(param,30,xlim=c(0,5),main=main,xlab=xlab,ylab="density",freq=FALSE)
    } else {
      hist(param,30,xlim=c(xlim[1],xlim[2]),main=main,
           xlab=xlab,ylab="density",freq=FALSE)
    }
      points(c(quantile(param,.025),quantile(param,.025),quantile(param,.975)),c(0,0,0),lwd=5,type="l")
      
    } 
  }  
}
  
 
