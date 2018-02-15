########################################


#' Print method for \code{bmeta} objects
#' 
#' Function to print output from function \code{bmeta}
#' 
#' 
#' @param x a \code{bmeta} object with results of the model
#' @param \dots other arguments
#' @author Tao Ding Gianluca Baio
print.bmeta<-function(x,...){
  print(x$mod,interval=c(0.025,0.975),digit=3)
}

