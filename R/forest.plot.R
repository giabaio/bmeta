##### FOREST PLOT #####




#' Function to create forest plot
#' 
#' A function to call package \code{forestplot} from R library and produce
#' forest plot using results from \code{bmeta}. The posterior estimate and
#' credible interval for each study are given by a square and a horizontal
#' line, respectively. The summary estimate is drawn as a diamond.
#' 
#' 
#' @param x a \code{bmeta} object with results of the model
#' @param title title of the plot
#' @param xlab title of the x-axis label
#' @param log estimates on natural scale is displayed by default. If TRUE, log
#' scale is used (i.e.  log odds ratio, log incidence rate ratio). For
#' continuous data, estimates are always presented on natural scale and users
#' do not need to specify this argument.
#' @param study.label label for each study and the summary estimate. See
#' details.
#' @param clip lower and upper limits for clipping credible intervals to arrows
#' @param lines selects the colour for the lines of the intervals. If the extra
#' option \code{add.null} is set to \code{TRUE}, then \code{lines} should be
#' specified as a two-element vector. If the user fails to do so, \code{bmeta}
#' will overwrite this setting and select suitable values.
#' @param box selects the colour for mean study-specific estimates. If the
#' extra option \code{add.null} is set to \code{TRUE}, then \code{box} should
#' be specified as a two-element vector. If the user fails to do so,
#' \code{bmeta} will overwrite this setting and select suitable values.
#' @param summary selects the colour for the pooled estimate
#' @param box.symb selects the symbol used to plot the mean. Options are "box"
#' (default) or "circle"
#' @param label.cex defines the size of the text for the label. Defaults at .8
#' of normal size
#' @param xlab.cex defines the size of the text for the x-label. Defaults at 1
#' of the normal size
#' @param ticks.cex defines the size of the text for the x-axis ticks. Defaults
#' at .8 of the normal size
#' @param ...  Additional arguments. Includes
#' 
#' - \code{add.null} = TRUE/FALSE. If set to true, adds a plot of the null
#' (no-pooling model) - \code{line.margin} = the distance between lines in case
#' multiple graphs are shown on the same plot - \code{box.size} = the size of
#' the summary box - \code{new.page} = TRUE/FALSE. If set to true, then a new
#' graph overwrite the existing one - \code{zero} (x-axis coordinate for zero
#' line. If you provide a vector of length 2 it will print a rectangle instead
#' of just a line. Default at 0 or 1 depending on log scale) - \code{legend} =
#' a legend for the multi-graph plot (including the null/no-pooling model)
#' @author Tao Ding Gianluca Baio
#' @keywords Forest plot
#' @examples
#' 
#' ### Read and format the data (binary)
#' data = read.csv(url("http://www.statistica.it/gianluca/bmeta/Data-bin.csv"))
#' 
#' ### List data for binary outcome 
#' data.list <- list(y0=data$y0,y1=data$y1,n0=data$n0,n1=data$n1) 
#' 
#' ### Select fixed-effects meta-analysis with normal prior for binary data 
#' x <- bmeta(data.list, outcome="bin", model="std.norm", type="fix")
#' 
#' ### Plot forest plot 
#' forest.plot(x)
#' 
#' ### Plot forest plot on log scale 
#' forest.plot(x,log=TRUE)
#' 
#' ### Select random-effects meta-analysis with t-distribution prior for binary
#' ### data
#' x <- bmeta(data.list, outcome="bin", model="std.dt", type="ran")
#' 
#' ### Plot 'two-line' forest plot showing estimates from both randome-effects 
#' ### model and no-pooling effects model for comparison
#' forest.plot(x,add.null=TRUE,title="Two-line forestplot for comparison")
#' 
#' 
#' ### Read and format the data (continuous)
#' data = read.csv(url("http://www.statistica.it/gianluca/bmeta/Data-ctns.csv"))
#' 
#' ### List data for continuous outcome
#' data.list <- list(y0=data$y0,y1=data$y1,se0=data$se0,se1=data$se1)  
#' 
#' ### Select fix-effects meta-analysis for studies reporting two arms separately
#' x <- bmeta(data=data.list,outcome="ctns",model="std.ta",type="fix")
#' 
#' ### Define for individual studies
#' study.label <- c(paste0(data$study,", ",data$year),"Summary estimate")
#' 
#' ### Produce forest plot with label for each study and control the lower and upper 
#' ### limits for clipping credible intervals to arrows   
#' forest.plot(x,study.label=study.label,clip=c(-7,4))
#' 
#' 
#' 
forest.plot <- function(x,title=NULL,xlab=NULL,log=FALSE,study.label=NULL,clip=c(-3,3),
                        lines="black",box="blue",summary="orange",box.symb="box",label.cex=.8,
                        xlab.cex=1,ticks.cex=.8,...) {
  
  # x = a bmeta object including the MCMC simulations for the results
  # title = possibly a string including the title for the forest plot
  # log = indicates whether to report the analysis on the log (TRUE) or the natural (FALSE, default) scale
  # study.label = a vector of study labels. If NULL then forest.plot will create one 
  # clip = graphical parameter to determine the range of values showed for the x-axis
  # lines = selects the colour for the lines of the intervals
  # box = selects the colour for the mean estimate
  # summary = selects the colour for the pooled estimate
  # box.symb = selects the symbol used to plot the mean. Options are "box" (default) or "circle"
  # label.cex = defines the size of the text for the label. Defaults at .8 of normal size
  # xlab.cex = defines the size of the text for the x-label. Defaults at 1 of the normal size
  # ticks.cex = defines the size of the text for the x-axis ticks. Defaults at .8 of the normal size
  # ... = additional arguments, including
  #  - add.null = TRUE/FALSE. If set to true, adds a plot of the null (no-pooling model)
  #  - line.margin = the distance between lines in case multiple graphs are shown on the same plot
  #  - box.size = the size of the summary box
  #  - new.page = TRUE/FALSE. If set to true, then a new graph overwrite the existing one
  #  - zero (x-axis coordinate for zero line. If you provide a vector of length 2 it will print 
  #         a rectangle instead of just a line. Default at 0 or 1 depending on log scale)
  #  - legend = a legend for the multi-graph plot
  
  exArgs <- list(...)  
  
  tab0 <- x$mod0$BUGSoutput$summary 
  tab <- x$mod$BUGSoutput$summary 
  
  if (is.null(study.label)) {
    study.label <- c(paste0("Study",1:x$data$S),"Summary")
  }
  
  is.sum=c(rep(FALSE,length(study.label)-1),TRUE)
  
  # Defines the default values for the original forestplot options
  if(exists("line.margin",where=exArgs)){line.margin=exArgs$line.margin} else {line.margin=0.1}
  if(exists("boxsize",where=exArgs)){boxsize=exArgs$boxsize} else {boxsize=0.4}
  if(exists("new_page",where=exArgs)){new_page=exArgs$new_page} else {new_page=TRUE}
  
  
  # If no xlabel has been specified, selects a suitable one
  if(is.null(xlab)) {
    cond <- c((x$outcome=="bin" & log==FALSE),(x$outcome=="bin" & log==TRUE),
              x$outcome=="ctns",(x$outcome=="count" & log==FALSE),
              (x$outcome=="count" & log==TRUE))
    labs <- c("OR","log(OR)","Mean","IRR","log(IRR)")
    xlab <- labs[which(cond==TRUE)]
  }
  
  
  ### different color/box symbol option for simple plots and plots showing results from two different models 
  # Checks whether the null (no-pooling) model should be also plotted & sets the relevant parameters for either cases
  if(exists("add.null",where=exArgs)) {add.null=exArgs$add.null} else {add.null=FALSE}
  if (add.null==FALSE){
    col <- fpColors(lines=lines,box=box,summary=summary)
    fn.ci_norm <- fpDrawNormalCI
    if(box.symb=="circle") {fn.ci_norm <- fpDrawCircleCI}
  } else {
    if(length(lines)!=2){lines <- c("grey","black")}
    if(length(box)!=2){box <-c("blue","darkred")}
    col <- fpColors(lines=c(lines[1],lines[2]),box=c(box[1],box[2]),summary=summary)
    fn.ci_norm <- c(fpDrawNormalCI, fpDrawCircleCI)     
  }
  
  txt_gp=fpTxtGp(xlab=grid::gpar(cex=xlab.cex),ticks=grid::gpar(cex=ticks.cex),label=grid::gpar(cex=label.cex))
  
  
  ### forest plot for fixed-effects model ###
  if (x$type=="fix") {
    
    if(x$outcome=="bin"){
      ind <- which((substr(rownames(tab0),1,3)=="lor")==TRUE)
      ind0 <- which((substr(rownames(tab0),1,2)=="or")==TRUE)
      ind1 <- which((substr(rownames(tab),1,3)=="rho")==TRUE)
      ind2 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      
      if(log==FALSE){
        xlog <- FALSE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=1}
        
        forestplot::forestplot(labeltext=study.label,mean=c(tab0[ind0,1],tab[ind1,1]),lower=c(tab0[ind0,3],tab[ind1,3]),
                               upper=c(tab0[ind0,7],tab[ind1,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,xlog=xlog,col=col,zero=zero,fn.ci_norm=fn.ci_norm)  
      } else {
        xlog <- TRUE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=0}
        
        forestplot::forestplot(labeltext=study.label,mean=c(tab0[ind,1],tab[ind2,1]),lower=c(tab0[ind,3],tab[ind2,3]),
                               upper=c(tab0[ind,7],tab[ind2,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,clip=clip,col=col,zero=zero,fn.ci_norm=fn.ci_norm)  
        
      }
      
    } 
    
    if(x$outcome=="ctns") {
      log==FALSE
      
      if(x$model=="std.ta" | x$model=="reg.ta"){
        ind0 <- which((substr(rownames(tab0),1,4)=="diff")==TRUE)
        ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)    
      }
      
      
      if(x$model=="std.mv"){
        ind0 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)
        ind1 <- which((substr(rownames(tab),1,2)=="mu")==TRUE)
      }
      
      if(x$model=="reg.mv"){
        ind0 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)
        ind1 <- which((substr(rownames(tab),1,5)=="alpha")==TRUE)
      }
      
      if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=0}
      
      forestplot::forestplot(labeltext=study.label,mean=c(tab0[ind0,1],tab[ind1,1]),lower=c(tab0[ind0,3],tab[ind1,3]),
                             upper=c(tab0[ind0,7],tab[ind1,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                             txt_gp=txt_gp,xlab=xlab,col=col,clip=clip,zero=zero,fn.ci_norm=fn.ci_norm)  
      
      
    }
    
    if(x$outcome=="count"){
      ind <- which((substr(rownames(tab0),1,4)=="lIRR")==TRUE)
      ind0 <- which((substr(rownames(tab0),1,3)=="IRR")==TRUE)
      ind1 <- which((substr(rownames(tab),1,3)=="IRR")==TRUE)
      ind2 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      
      if(log==FALSE){
        xlog <- FALSE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=1}
        
        forestplot::forestplot(labeltext=study.label,mean=c(tab0[ind0,1],tab[ind1,1]),lower=c(tab0[ind0,3],tab[ind1,3]),
                               upper=c(tab0[ind0,7],tab[ind1,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,xlog=xlog,col=col,zero=zero,fn.ci_norm=fn.ci_norm)  
      } else {
        xlog <- TRUE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=0}
        
        forestplot::forestplot(labeltext=study.label,mean=c(tab0[ind,1],tab[ind2,1]),lower=c(tab0[ind,3],tab[ind2,3]),
                               upper=c(tab0[ind,7],tab[ind2,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,clip=clip,col=col,zero=zero,fn.ci_norm=fn.ci_norm)  
        
      }
      
      
    }
    
    
  }
  
  
  ### forest plot for random-effects models ###   
  if(x$type=="ran"){
    if(x$outcome=="bin" & log==TRUE){
      ind0 <- which((substr(rownames(tab0),1,3)=="lor")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE) 
    } 
    
    if(x$outcome=="bin" & log==FALSE){
      ind0 <- which((substr(rownames(tab0),1,2)=="or")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="gamma")==TRUE)
      ind2 <- which((substr(rownames(tab),1,3)=="rho")==TRUE)
    }
    
    if(x$outcome=="ctns" & (x$model=="std.ta"|x$model=="reg.ta")){
      log=TRUE
      ind0 <- which((substr(rownames(tab0),1,4)=="diff")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE)
    }
    
    if(x$outcome=="ctns" & x$model=="std.mv"){
      log=TRUE
      ind0 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE)
    }  
    
    if(x$outcome=="ctns" & x$model=="reg.mv"){
      log==TRUE
      ind0 <- which((substr(rownames(tab0),1,2)=="mu")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="alpha")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE)
    }
    
    if(x$outcome=="count" & log==TRUE){
      ind0 <- which((substr(rownames(tab0),1,4)=="lIRR")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="delta")==TRUE)
      ind2 <- which((substr(rownames(tab),1,2)=="mu")==TRUE) 
    } 
    
    if(x$outcome=="count" & log==FALSE){
      ind0 <- which((substr(rownames(tab0),1,3)=="IRR")==TRUE)
      ind1 <- which((substr(rownames(tab),1,5)=="gamma")==TRUE)
      ind2 <- which((substr(rownames(tab),1,3)=="IRR")==TRUE)
    }
    
    
    if (add.null==FALSE){
      
      if((x$outcome=="bin" & log==FALSE)|(x$outcome=="count" & log==FALSE)){
        xlog <- FALSE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=1}
        
        forestplot::forestplot(labeltext=study.label,mean=c(tab[ind1,1],tab[ind2,1]),lower=c(tab[ind1,3],tab[ind2,3]),
                               upper=c(tab[ind1,7],tab[ind2,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,xlog=xlog,col=col,zero=zero,fn.ci_norm=fn.ci_norm)
      } else {
        xlog <- TRUE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=0}
        
        forestplot::forestplot(labeltext=study.label,mean=c(tab[ind1,1],tab[ind2,1]),lower=c(tab[ind1,3],tab[ind2,3]),
                               upper=c(tab[ind1,7],tab[ind2,7]),is.summary=is.sum,new_page=new_page,title=title,boxsize=boxsize,
                               txt_gp=txt_gp,xlab=xlab,clip=clip,col=col,zero=zero,fn.ci_norm=fn.ci_norm)  
        
      }
      
      
    } else {
      
      ### adding results from model (study-specific estimates and summary estimate) together
      comb1<-c(tab[ind1,1],tab[ind2,1]) 
      comb2<-c(tab[ind1,3],tab[ind2,3])
      comb3<-c(tab[ind1,7],tab[ind2,7])
      
      ### results from null model but how to avoid the summary estimate appear twice???   
      comb4<-c(tab0[ind0,1],500) 
      comb5<-c(tab0[ind0,3],500)
      comb6<-c(tab0[ind0,7],500)
      
      mean0<-as.matrix(comb4)
      lower0<-as.matrix(comb5)
      upper0<-as.matrix(comb6)
      
      mean1<-as.matrix(comb1)
      lower1<-as.matrix(comb2)
      upper1<-as.matrix(comb3)
      
      if((x$outcome=="bin" & log==FALSE)|(x$outcome=="count" & log==FALSE)){
        
        xlog <- FALSE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=1}
        if(exists("legend",where=exArgs)){
          legend=exArgs$legend
        } else {
          legend=c("Estimates from random-effects model","Estimates from no-pooling effects model")
        }
        
        forestplot::forestplot(labeltext=study.label,mean=cbind(mean1[,1],mean0[,1]),
                               lower=cbind(lower1[,1],lower0[,1]),upper=cbind(upper1[,1],upper0[,1]),
                               is.summary=is.sum,new_page=new_page,boxsize=boxsize,line.margin=line.margin,xlog=xlog,
                               fn.ci_norm=fn.ci_norm,col=col,title=title,legend=legend,txt_gp=txt_gp,
                               xlab=xlab,clip=clip,zero=zero)  
      } else{
        
        xlog <- TRUE
        if(exists("zero",where=exArgs)){zero=exArgs$zero} else {zero=0}
        if(exists("legend",where=exArgs)){
          legend=exArgs$legend
        } else {
          legend=c("Estimates from random-effects model","Estimates from no-pooling effects model")
        }
        
        forestplot::forestplot(labeltext=study.label,mean=cbind(mean1[,1],mean0[,1]),
                               lower=cbind(lower1[,1],lower0[,1]),upper=cbind(upper1[,1],upper0[,1]),
                               is.summary=is.sum,new_page=new_page,boxsize=boxsize,line.margin=line.margin,
                               fn.ci_norm=fn.ci_norm,col=col,title=title,legend=legend,txt_gp=txt_gp,
                               xlab=xlab,clip=clip,zero=zero)  
        
      }
      
      
      
      
    }
    
  }
  
}   

