
#' @S3method print acifit
#' @rdname fitaci
print.acifit <- function(x,...){
  
  cat("Result of fitaci.\n\n")
  
  cat("Data and predictions:\n")
  print(x$df)
  
  cat("\nParameters (at 25C):\n")
  
  print(x$pars)
  cat("\n")
  
}

#' @S3method summary acifit
#' @rdname fitaci
summary.acifit <- function(x,...){
  
  print.acifit(x, ...)
  
}