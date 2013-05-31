print.acifit <- function(x,...){
  
  cat("Result of fitaci.\n\n")
  
  cat("Data and predictions:\n")
  print(x$data)
  
  cat("\nParameters:\n")
  
  PARS <- summary(x$nlsfit)$coefficients[,1:2]
  print(PARS)
  cat("\n")
  
}


summary.acifit <- function(x,...){
  
  print.acifit(x, ...)
  
}