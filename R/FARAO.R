FARAO <- function(ca=380, ...){
	
	fx <- function(ca,...)optimize(OPTfun, interval=c(0,ca), maximum=TRUE,ca=ca,...)$maximum
	optimalcis <- mapply(fx,ca=ca,...)
	
	res <- as.data.frame(OPTfun(ci=optimalcis, retobjfun=FALSE, ca=ca, ...))
	
return(res)
}

FARABBT <- function(ca=380, ...){
	
	fx <- function(ca,...)optimize(BBobjfun, interval=c(0,ca), maximum=TRUE, ca=ca,...)$maximum
	optimalcis <- mapply(fx,ca=ca,...)
	
	res <- as.data.frame(t(mapply(BBobjfun, ci=optimalcis, retobjfun=FALSE, ca=ca, ...)))
	
return(res)
}

