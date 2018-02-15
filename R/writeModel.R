#### WRITE MODEL #####


#' A function to write a text file encoding the modelling assumptions
#' 
#' The \code{writeModel} function helps to select the proper model to be
#' contained in the 'model.file' for MCMC simulation based on users'
#' specifications.
#' 
#' 
#' @param outcome type of outcome that needs to be specified. For binary,
#' continuous and count data, 'bin', 'ctns' and 'count' need to be specified,
#' respectively.
#' @param model type of model that needs to be specified. There are 14 options:
#' 'std.norm','std.dt',
#' 'reg.norm','reg.dt','std.ta','std.mv','reg.ta','reg.mv','std','std.unif','std.hc','reg',
#' 'reg.unif','reg.hc'.
#' @param type model type---either fixed-effects("fix") or random-effects
#' model("ran") needs to be specified.
#' @param model.file file containing the appropriate model selected by user
#' @param data a data list containing information on observed data (including
#' moderators).
#' @author Tao Ding Gianluca Baio
writeModel<-function(outcome,model,type,model.file,data){
  
  ## Model selection
  ## Selects modules in the model code, according to the distributional assumptions and whether it is a standard meta-analysis or meta-regression
  
  
  
  
  ##############################################
  ## (i) standard meta-analysis for binary data
  
  
  sel.mod.bin.std.norm.fix<-" 
  ## a. binary standard fixed-effects meta-analysis with normal prior
  \n model {   
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
  } "
  
  
  
  
  
  sel.mod.bin.std.dt.fix<-" 
  ## b. binary standard fixed-effects meta-analysis with t-distribution prior
  \n model {  
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
  } "
  
  
  
  
  
  
  sel.mod.bin.std.norm.ran<-"
  ## c. binary standard random-effects meta-analysis with normal prior 
  \n model{  
  for (s in 1:S){
  y0[s]~dbin(pi0[s],n0[s])
  y1[s]~dbin(pi1[s],n1[s])
  logit(pi0[s])<-alpha[s]
  logit(pi1[s])<-alpha[s]+delta[s]
  
  alpha[s]~dnorm(0,0.0001)
  delta[s]~dnorm(mu,tau)
  gamma[s]<-exp(delta[s])
  }
  ###prior###
  mu~dnorm(0,0.0001)
  sigma~dunif(0,5)
  tau<-pow(sigma,-2)
  rho<-exp(mu)
  } "
  
  
  
  
  
  
  sel.mod.bin.std.dt.ran<-"
  ## d. binary standard random-effects meta-analysis with t-distribution prior
  \n model{  
  for (s in 1:S){
  y0[s]~dbin(pi0[s],n0[s])
  y1[s]~dbin(pi1[s],n1[s])
  logit(pi0[s])<-alpha[s]
  logit(pi1[s])<-alpha[s]+delta[s]
  
  delta[s]~dnorm(mu,tau)
  alpha[s]~dnorm(0,0.0001)
  gamma[s]<-exp(delta[s])
  }
  ###Prior###
  sigma~dunif(0,5)
  tau<-pow(sigma,-2)
  mu~dt(0,0.5,v)
  v~dunif(0,8)
  rho<-exp(mu)
  } " 
  
  
  
  
  #########################################
  ## (ii) meta-regression for binary data
  
  if(!is.null(data$J)) {
    
    sel.mod.bin.reg.norm.fix<-if(data$J==1) { " 
      ## a. binary fixed-effects meta-regression with normal prior
      \n model {  
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
      } "
      
    } else {
      " 
      ## a. binary fixed-effects meta-regression with normal prior
      \n model {  
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
      } "   
      
    }   
    
    
    
    
    
    
    
    
    sel.mod.bin.reg.dt.fix<-if(data$J==1) {
      "
      ## b. binary fixed-effects meta-regression with t-distribution prior
      \n model {  
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
      } "
    } else {
      "
      ## b. binary fixed-effects meta-regression with t-distribution prior
      \n model {  
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
    } "
    }
    
    
    
    sel.mod.bin.reg.norm.ran<-if(data$J==1) {
      "
      ## c. binary random-effects meta-regression with normal prior
      \n model {  
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta[s]
      
      alpha[s]~dnorm(0,0.0001)
      delta[s]~dnorm(mu,tau)
      gamma[s]<-exp(delta[s])
      }
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])
      
      mu~dnorm(0,0.0001)
      sigma~dunif(0,5)
      tau<-pow(sigma,-2)
      rho<-exp(mu)
      } "
    } else {
      "
      ## c. binary random-effects meta-regression with normal prior
      \n model {  
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta[s]
      
      alpha[s]~dnorm(0,0.0001)
      delta[s]~dnorm(mu,tau)
      gamma[s]<-exp(delta[s])
      }
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      
      mu~dnorm(0,0.0001)
      sigma~dunif(0,5)
      tau<-pow(sigma,-2)
      rho<-exp(mu)
      } "
    }
    
    
    
    
    sel.mod.bin.reg.dt.ran<-if(data$J==1) {
      " 
      ## d. binary random-effects meta-regression with t-distribution prior
      \n model {  
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta[s]
      
      delta[s]~dnorm(mu,tau)
      alpha[s]~dnorm(0,0.0001)
      gamma[s]<-exp(delta[s])
      }
      ###Prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])
      
      sigma~dunif(0,5)
      tau<-pow(sigma,-2)
      mu~dt(0,0.5,v)
      v~dunif(0,8)
      rho<-exp(mu)
      } "
    } else {
      "
      ## d. binary random-effects meta-regression with t-distribution prior
      \n model {  
      for (s in 1:S){
      y0[s]~dbin(pi0[s],n0[s])
      y1[s]~dbin(pi1[s],n1[s])
      logit(pi0[s])<-alpha[s]+Z0[s,]%*%beta0
      logit(pi1[s])<-alpha[s]+Z0[s,]%*%beta0+delta[s]
      
      delta[s]~dnorm(mu,tau)
      alpha[s]~dnorm(0,0.0001)
      gamma[s]<-exp(delta[s])
      }
      ###Prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      
      sigma~dunif(0,5)
      tau<-pow(sigma,-2)
      mu~dt(0,0.5,v)
      v~dunif(0,8)
      rho<-exp(mu)
      } "
    }
    
  }
  
  
  
  #################################################### 
  ## (iii) standard meta-analysis for continuous data
  
  sel.mod.ctns.std.ta.fix<-"
  ## a. continuous fixed-effects meta-analysis with data available for two arms separately
  \n model{   
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
  
  }"
  
  
  
  
  sel.mod.ctns.std.ta.ran<-"
  ## b. continuous random-effects meta-analysis with data available for two arms separately
  \n model{  
  for(s in 1:S) { 
  y0[s]~dnorm(alpha0[s],prec0[s])       
  y1[s]~dnorm(alpha1[s],prec1[s])       
  
  alpha1[s]<-alpha0[s]+delta[s]
  prec0[s]<-pow(se0[s],-2)
  prec1[s]<-pow(se1[s],-2)
  
  alpha0[s]~dnorm(0,0.0001)
  delta[s]~dnorm(mu,tau)
  }
  
  ###prior###
  mu~dnorm(0,0.0001)
  tau<-pow(sigma,-2)
  sigma~dunif(0,10)
  }"
  
  
  
  
  sel.mod.ctns.std.mv.fix<-"
  ## c. continuous fixed-effects meta-analysis for studies reporting mean difference and pooled variance
  \n model{  
  for(s in 1:S){
  y[s]~dnorm(mu,prec[s])
  
  }
  ###prior###
  mu~dnorm(0,0.0001)
  }"
  
  
  
  
  
  
  sel.mod.ctns.std.mv.ran<-"
  ## d. continuous random-effects meta-analysis for studies reporting mean difference and pooled variance
  \n model{  
  for(s in 1:S){
  y[s]~dnorm(delta[s],prec[s])        
  delta[s]~dnorm(mu,tau) 
  
  }
  ###prior###
  mu~dnorm(0,0.0001)
  tau<-pow(sigma,-2)
  sigma~dunif(0,10)
  }" 
  
  
  
  
  ############################################
  ## (iv) meta-regression for continuous data
  
  if(!is.null(data$J)) {
    sel.mod.ctns.reg.ta.fix<-if(data$J==1){ 
      "
      ## a. continuous fix-effects meta-regression with data available for two arms separately
      \n model{ 
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
      }"
    } else { 
      "
      ## a. continuous fix-effects meta-regression with data available for two arms separately
      \n model{ 
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
      }"
    }
    
    
    
    sel.mod.ctns.reg.ta.ran<-if(data$J==1){
      "
      ## b. continuous random-effects meta-regression with data available for two arms separately
      \n model{ 
      for(s in 1:S) { 
      y0[s]~dnorm(alpha0[s],prec0[s])       
      y1[s]~dnorm(alpha1[s],prec1[s])       
      
      alpha1[s]<-alpha0[s]+delta[s]
      prec0[s]<-pow(se0[s],-2)
      prec1[s]<-pow(se1[s],-2)
      
      alpha0[s]<-gamma0[s]+Z0[s,]%*%beta0
      gamma0[s]~dnorm(0,0.0001)
      delta[s]~dnorm(mu,tau)
      
      }
      
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])     
      mu~dnorm(0,0.0001)
      tau<-pow(sigma,-2)
      sigma~dunif(0,10)
      
      }" 
    } else {
      "
      ## b. continuous random-effects meta-regression with data available for two arms separately
      \n model{ 
      for(s in 1:S) { 
      y0[s]~dnorm(alpha0[s],prec0[s])       
      y1[s]~dnorm(alpha1[s],prec1[s])       
      
      alpha1[s]<-alpha0[s]+delta[s]
      prec0[s]<-pow(se0[s],-2)
      prec1[s]<-pow(se1[s],-2)
      
      alpha0[s]<-gamma0[s]+Z0[s,]%*%beta0
      gamma0[s]~dnorm(0,0.0001)
      delta[s]~dnorm(mu,tau)
      
      }
      
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])      
      mu~dnorm(0,0.0001)
      tau<-pow(sigma,-2)
      sigma~dunif(0,10)
      
    }"   
    }  
    
    
    
    
    sel.mod.ctns.reg.mv.fix<-if(data$J==1){
      "
      ## d. continuous fixed-effects meta-regression for studies reporting mean difference and pooled variance 
      \n model{ 
      for(s in 1:S){
      y[s]~dnorm(delta[s],prec[s])        
      delta[s]<-alpha+Z0[s,]%*%beta0
      
      }
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])       
      alpha~dnorm(0,0.0001)
      }"
    } else {
      "
      ## d. continuous fixed-effects meta-regression for studies reporting mean difference and pooled variance
      \n model{ 
      for(s in 1:S){ 
      y[s]~dnorm(delta[s],prec[s])        
      delta[s]<-alpha+Z0[s,]%*%beta0
      
      }
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])       
      alpha~dnorm(0,0.0001)
      
    }"
    }
    
    
    
    
    sel.mod.ctns.reg.mv.ran<-if(data$J==1){
      "
      ## d. continuous random-effects meta-regression for studies reporting mean difference and pooled variance
      \n model{ 
      for(s in 1:S){
      y[s]~dnorm(delta[s],prec[s])        
      delta[s]<-alpha[s]+Z0[s,]%*%beta0
      alpha[s]~dnorm(mu,tau)
      
      }
      ###prior###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,])       
      mu~dnorm(0,0.0001)
      tau<-pow(sigma,-2)
      sigma~dunif(0,10)
      }"
    } else {
      "
      ## d. continuous random-effects meta-regression for studies reporting mean difference and pooled variance
      \n model{  
      for(s in 1:S){
      y[s]~dnorm(delta[s],prec[s])        
      delta[s]<-alpha[s]+Z0[s,]%*%beta0
      alpha[s]~dnorm(mu,tau)
      
      }
      ###prior###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])       
      mu~dnorm(0,0.0001)
      tau<-pow(sigma,-2)
      sigma~dunif(0,10)
      
    }" 
    }
  }
  
  
  
  
  ##############################################
  ## (v) standard meta-analysis for count data
  
  sel.mod.count.std.fix<-"
  ## a. count fixed-effects meta-analysis
  \n model{  
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
  }"
  
  
  
  
  
  sel.mod.count.std.unif.ran<-"
  ## b. count random-effects meta-analysis with uniform prior
  \n model{ 
  for(s in 1:S){
  y0[s]~dpois(lambda0[s])
  y1[s]~dpois(lambda1[s])
  
  log(lambda0[s])<-xi0[s]+log(p0[s])
  log(lambda1[s])<-xi1[s]+log(p1[s])
  
  xi0[s]~dunif(-5,5)
  xi1[s]<-xi0[s]+delta[s]
  delta[s]~dnorm(mu,tau)
  gamma[s]<-exp(delta[s])
  }
  ### Prior ###
  mu~dnorm(0,0.0001)
  tau<-pow(sigma,-2)
  sigma~dunif(0,10)
  IRR<-exp(mu)
  }"
  
  
  
  
  sel.mod.count.std.hc.ran<-"
  ## c. count random-effects meta-analysis with halfcauchy prior  
  \n model{ 
  for(s in 1:S){
  y0[s]~dpois(lambda0[s])
  y1[s]~dpois(lambda1[s])
  
  log(lambda0[s])<-xi0[s]+log(p0[s])
  log(lambda1[s])<-xi1[s]+log(p1[s])
  
  xi0[s]~dunif(-5,5)
  xi1[s]<-xi0[s]+delta[s]
  delta[s]~dnorm(mu,tau)
  gamma[s]<-exp(delta[s])
  }
  ### Prior ###
  mu~dnorm(0,0.0001)
  tau<-pow(sigma,-2)
  
  sigma<-abs(z.sigma)/pow(epsilon.sigma,0.5)
  z.sigma~dnorm(0,tau.z.sigma)
  epsilon.sigma~dgamma(0.5,0.5)
  tau.z.sigma<-pow(B.sigma,-2)
  B.sigma~dunif(0,0.5)
  
  IRR<-exp(mu) 
  }"
  
  
  
  
  
  
  ##########################################
  ## (vi) meta-regression for count data 
  
  if(!is.null(data$J)) {
    sel.mod.count.reg.fix<-if(data$J==1){
      "
      ## a. count fixed-effects meta-regression
      \n model{ 
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
      }"
    } else {
      "
      ## a. count fixed-effects meta-regression
      \n model{ 
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
      
      }" 
    }
    
    
    
    
    
    sel.mod.count.reg.unif.ran<-if(data$J==1){
      "
      ## b. count random-effects meta-regression with uniform prior
      \n model{ 
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      
      log(lambda0[s])<-xi0[s]+log(p0[s])
      log(lambda1[s])<-xi1[s]+log(p1[s])
      
      xi0[s]<-theta[s]+Z0[s,]%*%beta0   
      theta[s]~dunif(-5,5)
      xi1[s]<-xi0[s]+delta[s]
      delta[s]~dnorm(mu,tau)
      gamma[s]<-exp(delta[s])
      }
      ### Prior ###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,]) 
      mu~dnorm(0,0.0001)
      tau<-pow(sigma,-2)
      sigma~dunif(0,10)
      IRR<-exp(mu)
      }"
    } else {
      "
      ## b. count random-effects meta-regression with uniform prior
      \n model{ 
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      
      log(lambda0[s])<-xi0[s]+log(p0[s])
      log(lambda1[s])<-xi1[s]+log(p1[s])
      
      xi0[s]<-theta[s]+Z0[s,]%*%beta0   
      theta[s]~dunif(-5,5)
      xi1[s]<-xi0[s]+delta[s]
      delta[s]~dnorm(mu,tau)
      gamma[s]<-exp(delta[s])
      }
      ### Prior ###
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      mu~dnorm(0,0.0001)
      tau<-pow(sigma,-2)
      sigma~dunif(0,10)
      IRR<-exp(mu)
    }" 
    }
    
    
    
    
    
    sel.mod.count.reg.hc.ran<-if(data$J==1){"
      ## c. count random-effects meta-regression with halfcauchy prior
      \n model{ 
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      
      log(lambda0[s])<-xi0[s]+log(p0[s])
      log(lambda1[s])<-xi1[s]+log(p1[s])
      
      xi0[s]<-theta[s]+Z0[s,]%*%beta0   
      theta[s]~dunif(-5,5)
      xi1[s]<-xi0[s]+delta[s]
      delta[s]~dnorm(mu,tau)
      gamma[s]<-exp(delta[s])                                                                                       
      }
      ### Prior ###
      beta0[1:J]~dnorm(m.beta0[],tau.beta0[,]) 
      mu~dnorm(0,0.0001)
      tau<-pow(sigma,-2)
      
      sigma<-abs(z.sigma)/pow(epsilon.sigma,0.5)
      z.sigma~dnorm(0,tau.z.sigma)
      epsilon.sigma~dgamma(0.5,0.5)
      tau.z.sigma<-pow(B.sigma,-2)
      B.sigma~dunif(0,0.5)
      IRR<-exp(mu)
      
      }"
    } else {
      "
      ## c. count random-effects meta-regression with halfcauchy prior
      \n model{ 
      for(s in 1:S){
      y0[s]~dpois(lambda0[s])
      y1[s]~dpois(lambda1[s])
      
      log(lambda0[s])<-xi0[s]+log(p0[s])
      log(lambda1[s])<-xi1[s]+log(p1[s])
      
      xi0[s]<-theta[s]+Z0[s,]%*%beta0   
      theta[s]~dunif(-5,5)
      xi1[s]<-xi0[s]+delta[s]
      delta[s]~dnorm(mu,tau)
      gamma[s]<-exp(delta[s])
      }
      ### Prior ###
      
      beta0[1:J]~dmnorm(m.beta0[],tau.beta0[,])
      mu~dnorm(0,0.0001)
      tau<-pow(sigma,-2)
      
      sigma<-abs(z.sigma)/pow(epsilon.sigma,0.5)
      z.sigma~dnorm(0,tau.z.sigma)
      epsilon.sigma~dgamma(0.5,0.5)
      tau.z.sigma<-pow(B.sigma,-2)
      B.sigma~dunif(0,0.5)
      
      IRR<-exp(mu)         
      
      }"  
    }
    
  }
  
  
  
  
  
  
  ###############################################################
  ## Combines the required modules in a singel model code ##
  mod.outcome<-c("bin","ctns","count")
  mod.sel<-c("std.norm","std.dt","reg.norm","reg.dt",
             "std.ta","std.mv","reg.ta","reg.mv",
             "std","std.unif","std.hc","reg","reg.unif","reg.hc")
  mod.type<-c("fix","ran")
  
  lab.outcome<-c("Model for binary data ", "Model for continuous data ", "Model for count data ")
  lab.sel<-c("Meta-analysis with normal prior ","Meta-analysis with t-distribution prior ",
             "Meta-regression with normal prior ", "Meta-regression with t-distribution prior ",
             "Meta-analysis with data avalaible for two arms separately ", 
             "Meta-analysis for studies reporting mean difference and pooled variance ",
             "Meta-regression with data avalaible for two arms separately ", 
             "Meta-regression for studies reporting mean difference and pooled variance ",
             "Meta-analysis ","Meta-analysis with uniform prior ","Meta-analysis with halfcauchy prior ",
             "Meta-regression ","Meta-regression with uniform prior ","Meta-regression with halfcauchy prior ")
  lab.type<-c("Fixed-effects ", "Random-effects ")
  time <- Sys.time()
  
  for (mo in 1:length(mod.outcome)){
    for (ms in 1:length(mod.sel)){
      for (mt in 1:length(mod.type)){
        if (outcome==mod.outcome[mo] & model==mod.sel[ms] & type==mod.type[mt]){
          summary.text<- paste0("#",lab.outcome[mo],lab.type[mt],lab.sel[ms])
          print.time <- paste0("#",time)
          txt<- paste0("summary.text, print.time, sel.mod.", mod.outcome[mo], ".", mod.sel[ms], ".", mod.type[mt])  
          cmd<- paste0("model<- paste(",txt,")")
          eval(parse(text=cmd))
        }
      }
    }
  }
  
  filein <- model.file
  writeLines(model, con=filein)
  
}
