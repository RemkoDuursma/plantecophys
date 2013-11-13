
#' Predict from a non-linear regression object
#' @description Predicts from a fitted non-linear regression object, as returned by \code{\link{nls}}. 
#' Also calculates a confidence interval using a case bootstrap resampling method (using \code{\link{bootCase}}) from the \code{car} package.
#' @param object An object of class 'nls'
#' @param from 

#'@examples
#'
#' f <- function(t, k=0.5)1 - exp(-k*t)
#' x <- seq(0,12,length=101)
#' y <- f(x) + rnorm(101,sd=0.1)
#' dfr <- data.frame(x=x,y=y)
#' nls1 <- nls(y ~ b0 + A*(1-exp(-k*x)), start=list(b0=0,A=1,k=1), data=dfr)
#' 
#' p <- predict_nls(nls1, from=min(x), to=max(x), interval="confidence")
#' 
#' plot(x,y)
#' with(p, {
#'   lines(x, pred, lty=5)
#'   lines(x, lwr, lty=4, col="red")
#'   lines(x, upr, lty=4, col="red")
#' })
#' P <- predictNLS(nls1, newdata = data.frame(x=p$x))
#' lines(p$x, P$summary[,5], col="blue", lty=5)
#' lines(p$x, P$summary[,6], col="blue", lty=5)
# 


predict_nls <- function(object, from=NULL, to=NULL, x=NULL,interval = c("none", "confidence"), 
                        level=0.95, 
                        n=101, nboot=999, add=TRUE, ...){
  
  
  interval <- match.arg(interval)
  
  # find out name of X variable in nls object
  getx <- function(object){
    coefs <- names(coef(object))
    yvar <- as.character(object$m$formula()[[2]])
    els <- ls(object$m$getEnv())
    xvar <- setdiff(els, c(coefs,yvar))
    return(xvar)
  }
  
  if(is.null(x)){
    xi <- seq(from,to, length=n)
  } else {
    xi <- x
  }
  
  f <- object$m$formula()[[3]]
  
  
  dfr <- as.data.frame(xi)
  names(dfr) <- getx(object)
  pred <- predict(object, dfr)
  l <- list()
  l$x <- xi
  l$pred <- pred
  class(l) <- "nlspred"
  
  if(interval == "confidence"){
    d <- (1-level)/2
    b <- bootCase(object, B=nboot)
    parnames <- names(coef(object))
    npars <- length(parnames)
    preds <- list()
    
    for(i in seq_len(nboot)){
      
      p <- c(as.list(b[i,]),list(xi))
      names(p)[length(p)] <- getx(object)
      preds[[i]] <- eval(f,p)
    }
    preds <- do.call(rbind,preds)
    
    l$lwr <- apply(preds,2,function(x)quantile(x, d))
    l$upr <- apply(preds,2,function(x)quantile(x, level+d))
    
    l$boot <- b
  }
  
  return(l)
  
  
}


