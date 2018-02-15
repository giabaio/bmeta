###### Funnel Plot #####




#' Funnel plot to examine publication bias
#' 
#' Function to examine publication bias. For both fixed- and random-effects
#' models, estimates from no-pooling effects model are used as study-specific
#' estimates. For random-effects models, the corresponding fixed-effects models
#' are implemented at background to obtain pooled estimate. For example, if
#' users call \code{bmeta} to run random-effects meta-analysis with normal
#' prior, fixed-effects meta-analysis with normal prior are implemented at
#' background to obtain pooled estimate for graphing. In the absence of
#' publication and heterogeneity, the scatter resembles a symmetrical funnel
#' and the triangle area formed by connecting the centred summary estimate with
#' its 2.5\% and 97.5\% quantiles on either side includes about 95\% of the
#' studies if the fixed-effects model assumption holds (i.e. all the studies
#' estimate the same effect).
#' 
#' 
#' @param x a \code{bmeta} object with results of the model
#' @param xlab title of x-axis. If unspecified, the function sets an
#' appropriate lable by default.
#' @param ylab title of x-axis. If unspecified, the function sets an
#' appropriate lable by default.
#' @param title title of the plot if specified
#' @param xlim horozontal limits of the plot region. If unspecified, the
#' function sets the horizontal plot limits to (-6,6).
#' @author Tao Ding Gianluca Baio
#' @keywords Funnel plot
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
#' ### using output from bmeta to produce funnel plot 
#' funnel.plot(x)
#' 
#' ### using output from bmeta and specify title of the plot 
#' funnel.plot(x,title="funnel plot")
#' 
#' ### using output from bmeta and specify the limit of x-axis and title
#' funnel.plot(x,title="funnel plot",xlim=c(-2,1))
#' 
funnel.plot<-function(x,xlab=NULL,ylab=NULL,title=NULL,xlim=NULL){
  
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  
  tab0 <- x$mod0$BUGSoutput$summary
  tab <- x$mod$BUGSoutput$summary
   
  if (is.null(xlab)){
    xlab="effect"
  }
  
  if (is.null(ylab)){
    ylab="size"
  }
  
  if (x$type=="fix") {
    
    if(x$outcome=="bin"){
      ind0 <- which((substr(rownames(tab0),1,3)=="lor")==TRUE)
      ind <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      
    } 
    
    if(x$outcome=="ctns") {
      
      if(x$model=="std.ta" | x$model=="reg.ta"){
        ind0 <- which((substr(rownames(tab0),1,4)=="diff")==TRUE)
        ind <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      }
      
      if(x$model=="std.mv"){
        ind0 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)  
        ind <- which((substr(rownames(tab),1,2)=="mu")==TRUE) 
      } 
      
      if(x$model=="reg.mv"){
        ind0 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)  
        ind <- which((substr(rownames(tab),1,5)=="alpha")==TRUE) 
      }
      
    }
    
    if(x$outcome=="count"){
      ind0 <- which((substr(rownames(tab0),1,4)=="lIRR")==TRUE)      
      ind <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
    } 
    
    if(is.null(xlim)){
      plot(tab0[ind0,1],tab0[ind0,2],ylim=c(max(tab0[ind0,2]),0),
           xlab=xlab,ylab=ylab,main=title,xlim=c(-6,6))
    } else {
      plot(tab0[ind0,1],tab0[ind0,2],ylim=c(max(tab0[ind0,2]),0),
           xlab=xlab,ylab=ylab,main=title,xlim=c(xlim[1],xlim[2]))
    }
    abline(v=tab[ind,1])
    segments(tab[ind,1]-(max(tab0[ind0,2])*1.96),max(tab0[ind0,2]),
             tab[ind,1],0,lty=2)
    segments(tab[ind,1]+(max(tab0[ind0,2])*1.96),max(tab0[ind0,2]),
             tab[ind,1],0,lty=2)
    
  }
  
  if(x$type=="ran"){
    
    if(!is.null(x$data$J)) {
      threshold<- 1e-12
      Z0<-x$data$Z0
      for (j in 1:x$data$J) {
        if (mean(x$data$Z0[,j])>threshold) {
          Z0[,j]<-scale(x$data$Z0[,j],scale=FALSE)
        }
      }      
    }
    
    file1 <- "file1.txt"
    
    if (x$outcome=="bin" & x$model=="std.norm") {
      cat("model {   
          for (s in 1:S){
          y0[s]~dbin(pi0[s],n0[s])
          y1[s]~dbin(pi1[s],n1[s])
          logit(pi0[s])<-alpha[s]
          logit(pi1[s])<-alpha[s]+delta
          alpha[s]~dnorm(0,0.0001)
          }
          ###prior###
          delta~dnorm(0,0.0001)
          rho<-exp(delta)
    }",file=file1)
    dataJags1<-list(S=x$data$S,y0=x$data$y0,n0=x$data$n0,y1=x$data$y1,n1=x$data$n1)
    param1 <- c("delta")
    inits1 <- function(){
      list(delta=rnorm(1),alpha=rnorm(x$data$S))
    }
  }
    
    if (x$outcome=="bin" & x$model=="reg.norm") {

      if(x$data$J==1) { 
        cat(" 
        model {  
        for (s in 1:S){
        y0[s]~dbin(pi0[s],n0[s])
        y1[s]~dbin(pi1[s],n1[s])
        logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
        logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta
        alpha[s]~dnorm(0,0.0001)
        }
        ###prior###
        beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])
        delta~dnorm(0,0.0001)
        rho<-exp(delta)
        } ",file=file1)
      
      } else {
        cat(" 
        model {  
        for (s in 1:S){
        y0[s]~dbin(pi0[s],n0[s])
        y1[s]~dbin(pi1[s],n1[s])
        logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
        logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta
        alpha[s]~dnorm(0,0.0001)
        }
        ###prior###
        beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
        delta~dnorm(0,0.0001)
        rho<-exp(delta)
        }",file=file1)   
      
    }         
      
      dataJags1<-list(S=x$data$S,y0=x$data$y0,n0=x$data$n0,y1=x$data$y1,n1=x$data$n1,
                      J=x$data$J,Z0=Z0,m.beta0=rep(0,x$data$J), 
                      tau.beta0=.0001*diag(x$data$J))
      param1 <- c("delta")
      inits1 <- function(){
        list(delta=rnorm(1),beta0=rnorm(x$data$J),alpha=rnorm(x$data$S))
      }
}
    
    if (x$outcome=="bin" & x$model=="std.dt") {
      cat("model {   
          for (s in 1:S){
          y0[s]~dbin(pi0[s],n0[s])
          y1[s]~dbin(pi1[s],n1[s])
          logit(pi0[s])<-alpha[s]
          logit(pi1[s])<-alpha[s]+delta
          alpha[s]~dnorm(0,0.0001)
    }
          ###prior###
          delta~dt(0,0.5,v)
          v~dunif(0,8)
          rho<-exp(delta)
  }",file=file1)
    dataJags1<-list(S=x$data$S,y0=x$data$y0,n0=x$data$n0,y1=x$data$y1,n1=x$data$n1)
    param1 <- c("delta")
    inits1 <- function(){
      list(v=runif(1),alpha=rnorm(x$data$S))
    } 
  }
   
    if (x$outcome=="bin" & x$model=="reg.dt") {
      
      if(x$data$J==1) { 
        cat(" 
            model {  
            for (s in 1:S){
            y0[s]~dbin(pi0[s],n0[s])
            y1[s]~dbin(pi1[s],n1[s])
            logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
            logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta
            alpha[s]~dnorm(0,0.0001)
            }
            ###prior###
            beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])
            delta~dt(0,0.5,v)
            v~dunif(0,8)
            rho<-exp(delta)
            } ",file=file1)
      
      } else {
        cat(" 
            model {  
            for (s in 1:S){
            y0[s]~dbin(pi0[s],n0[s])
            y1[s]~dbin(pi1[s],n1[s])
            logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
            logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta
            alpha[s]~dnorm(0,0.0001)
            }
            ###prior###
            beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
            delta~dt(0,0.5,v)
            v~dunif(0,8)
            rho<-exp(delta)
            }",file=file1)   
      
    }         
      dataJags1<-list(S=x$data$S,y0=x$data$y0,n0=x$data$n0,y1=x$data$y1,n1=x$data$n1,
                      J=x$data$J,Z0=Z0,m.beta0=rep(0,x$data$J), 
                      tau.beta0=.0001*diag(x$data$J))
      param1 <- c("delta")
      inits1 <- function(){
        list(v=runif(1),beta0=rnorm(x$data$J),alpha=rnorm(x$data$S))
      }
    }
    
    
    if (x$outcome=="ctns" & x$model=="std.ta"){
      cat("model{   
          for(s in 1:S) { 
          y0[s]~dnorm(alpha0[s],prec0[s])       
          y1[s]~dnorm(alpha1[s],prec1[s])       
          alpha1[s]<-alpha0[s]+delta
          prec0[s]<-pow(se0[s],-2)
          prec1[s]<-pow(se1[s],-2)
          alpha0[s]~dnorm(0,0.0001)
        }
        ###prior###
        delta~dnorm(0,0.0001)
      }",file=file1) 
      dataJags1<-list(S=x$data$S,y0=x$data$y0,se0=x$data$se0,y1=x$data$y1,se1=x$data$se1)
      param1 <- c("delta")
      inits1 <- function(){
        list(delta=rnorm(1),alpha0=rnorm(x$data$S))
      } 
    }
    
    
    if (x$outcome=="ctns" & x$model=="reg.ta"){
      if(x$data$J==1) {
          cat("model{ 
            for(s in 1:S) { 
              y0[s]~dnorm(alpha0[s],prec0[s])       
              y1[s]~dnorm(alpha1[s],prec1[s])       
              alpha1[s]<-alpha0[s]+delta
              prec0[s]<-pow(se0[s],-2)
              prec1[s]<-pow(se1[s],-2)
              alpha0[s]<-gamma0[s]+Z0[s,]%*%beta0
              gamma0[s]~dnorm(0,0.0001)
            }
            ###prior###
            beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])      
            delta~dnorm(0,0.0001)
          }",file=file1)
    } else { 
          cat("
          model{ 
            for(s in 1:S) { 
              y0[s]~dnorm(alpha0[s],prec0[s])       
              y1[s]~dnorm(alpha1[s],prec1[s])       
              
              alpha1[s]<-alpha0[s]+delta
              prec0[s]<-pow(se0[s],-2)
              prec1[s]<-pow(se1[s],-2)
              alpha0[s]<-gamma0[s]+Z0[s,]%*%beta0
              gamma0[s]~dnorm(0,0.0001)
            }
            ###prior###
            beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])     
            delta~dnorm(0,0.0001)
          }",file=file1) 
      
}
      dataJags1<-list(S=x$data$S,y0=x$data$y0,se0=x$data$se0,y1=x$data$y1,se1=x$data$se1,
                      J=x$data$J,Z0=Z0,m.beta0=rep(0,x$data$J), 
                      tau.beta0=.0001*diag(x$data$J))
      param1 <- c("delta")
      inits1 <- function(){
        list(delta=rnorm(1),beta0=rnorm(x$data$J),gamma0=rnorm(x$data$S))
      } 
}
    
    
    if (x$outcome=="ctns" & x$model=="std.mv"){
      cat("model{  
        for(s in 1:S){
          y[s]~dnorm(mu,prec[s])
        }
        ###prior###
        mu~dnorm(0,0.0001)
      }",file=file1)
      dataJags1<-list(S=x$data$S,y=x$data$y,prec=x$data$prec)
      param1 <- c("mu")
      inits1 <- function(){
        list(mu=rnorm(1))
      } 
    }
    
    
    if (x$outcome=="ctns" & x$model=="reg.mv"){
      if(x$data$J==1){
        cat("
        model{ 
        for(s in 1:S){
        y[s]~dnorm(delta[s],prec[s])        
        delta[s]<-alpha+Z0[s,]%*%beta0
        
        }
        ###prior###
        beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])       
        alpha~dnorm(0,0.0001)
        }",file=file1)
    } else {
      cat("
      model{ 
      for(s in 1:S){ 
      y[s]~dnorm(delta[s],prec[s])        
      delta[s]<-alpha+Z0[s,]%*%beta0
      
      }
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])       
      alpha~dnorm(0,0.0001)
      
      }",file=file1)
    }
      dataJags1<-list(S=x$data$S,y=x$data$y,prec=x$data$prec,J=x$data$J,Z0=Z0,
                      m.beta0=rep(0,x$data$J),tau.beta0=.0001*diag(x$data$J))
      param1 <- c("alpha")
      inits1 <- function(){
        list(alpha=rnorm(1),beta0=rnorm(x$data$J))
      } 
    }
    
    
    if (x$outcome=="count" & (x$model=="std.unif"|x$model=="std.hc")){
       cat("model{  
        for(s in 1:S){
          y0[s]~dpois(lambda0[s])
          y1[s]~dpois(lambda1[s])
          log(lambda0[s])<-xi0[s]+log(p0[s])
          log(lambda1[s])<-xi1[s]+log(p1[s])
          xi0[s]~dunif(-5,5)
          xi1[s]<-xi0[s]+delta      ### delta is the log rate-ratio
        }
        ### Prior ###
        delta~dnorm(0,0.0001)
        IRR<-exp(delta)
      }",file=file1)  
      dataJags1<-list(S=x$data$S,y0=x$data$y0,p0=x$data$p0,y1=x$data$y1,p1=x$data$p1)
      param1 <- c("delta")
      inits1 <- function(){
        list(delta=rnorm(1),xi0=runif(x$data$S))
      } 
    }
    
    
    if (x$outcome=="count" & (x$model=="reg.unif"|x$model=="reg.hc")){
      if(x$data$J==1){
        cat("
        model{ 
        for(s in 1:S){
        y0[s]~dpois(lambda0[s])
        y1[s]~dpois(lambda1[s])
        log(lambda0[s])<-xi0[s]+log(p0[s])
        log(lambda1[s])<-xi1[s]+log(p1[s])
        xi0[s]<-theta[s]+Z0[s,]%*%beta0   
        theta[s]~dunif(-5,5)
        xi1[s]<-xi0[s]+delta
        }
        ### Prior ###
        beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])
        delta~dnorm(0,0.0001)
        IRR<-exp(delta)
        }",file=file1)
    } else {
      cat("
      model{ 
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      log(lambda0[s])<-xi0[s]+log(p0[s])
      log(lambda1[s])<-xi1[s]+log(p1[s])
      xi0[s]<-theta[s]+Z0[s,]%*%beta0   
      theta[s]~dunif(-5,5)
      xi1[s]<-xi0[s]+delta
      }
      ### Prior ###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      delta~dnorm(0,0.0001)
      IRR<-exp(delta)
      }",file=file1)
    }
      dataJags1<-list(S=x$data$S,y0=x$data$y0,p0=x$data$p0,y1=x$data$y1,p1=x$data$p1,
                      J=x$data$J,Z0=Z0,m.beta0=rep(0,x$data$J),
                      tau.beta0=.0001*diag(x$data$J))
      param1 <- c("delta")
      inits1 <- function(){
        list(delta=rnorm(1),beta0=rnorm(x$data$J),theta=runif(x$data$S))
      } 
}
  
  
    # Now also run the fixed-effects model
    working.dir <- getwd()  # defines the current directory to work in
    quiet(mod1 <- R2jags::jags(dataJags1,inits1,param1,model.file="file1.txt",n.chains=2,
                               n.iter=10000,n.burnin=5000,n.thin=1,DIC=TRUE,
                               working.directory=working.dir, progress.bar="text"))
    
    unlink(file1,force=TRUE)  # Now deletes the file "file1.txt" from current directory
    
    
    tab1 <- mod1$BUGSoutput$summary
    
    if(x$outcome=="bin"){
      ind1 <- which((substr(rownames(tab0),1,3)=="lor")==TRUE)
      ind2 <- which((substr(rownames(tab1),1,5)=="delta")==TRUE) 
    } 
    
    
    if(x$outcome=="ctns" & (x$model=="std.ta" | x$model=="reg.ta")){
      ind1 <- which((substr(rownames(tab0),1,4)=="diff")==TRUE)
      ind2 <- which((substr(rownames(tab1),1,5)=="delta")==TRUE)
    }
    
    
    if(x$outcome=="ctns" & x$model=="std.mv"){
      ind1 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)
      ind2 <- which((substr(rownames(tab1),1,2)=="mu")==TRUE)
    }
    
    if(x$outcome=="ctns" & x$model=="reg.mv"){
      ind1 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)
      ind2 <- which((substr(rownames(tab1),1,5)=="alpha")==TRUE)
    }
    
    
    if(x$outcome=="count"){
      ind1 <- which((substr(rownames(tab0),1,4)=="lIRR")==TRUE)
      ind2 <- which((substr(rownames(tab1),1,5)=="delta")==TRUE) 
    } 
    
    
    if(is.null(xlim)){
      plot(tab0[ind1,1],tab0[ind1,2],ylim=c(max(tab0[ind1,2]),0),
         xlab=xlab,ylab=ylab,main=title,xlim=c(-6,6))
    } else {
      plot(tab0[ind1,1],tab0[ind1,2],ylim=c(max(tab0[ind1,2]),0),
           xlab=xlab,ylab=ylab,main=title,xlim=c(xlim[1],xlim[2]))
    }
      abline(v=tab1[ind2,1])
    
      segments(tab1[ind2,1]-(max(tab0[ind1,2])*1.96),max(tab0[ind1,2]),
             tab1[ind2,1],0,lty=2)
      segments(tab1[ind2,1]+(max(tab0[ind1,2])*1.96),max(tab0[ind1,2]),
             tab1[ind2,1],0,lty=2)
    
   } 
  
}

