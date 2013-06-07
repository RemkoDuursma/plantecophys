

  
fitaci <- function(dat, nlsmethod="default"){
  
  if(!"PAR" %in% names(dat)){
    dat$PAR <- 1800
    warning("PAR not in dataset; assumed PAR = 1800.")
  }
    
  acifun_wrap <- function(...){
    r <- Aci(...)
    r$ALEAF
  }
  
  nlsfit <- nls2(Photo ~ acifun_wrap(Ci, PPFD=PAR, Vcmax=Vcmax, Jmax=Jmax, Rd=Rd, Tleaf=Tleaf),
                data=dat, algorithm=nlsmethod,
                start=list(Vcmax=80, Jmax=160, Rd=1))

  p <- coef(nlsfit)
  acirun <- with(dat, Aci(Ci, PPFD=PAR, Vcmax=p[[1]], Jmax=p[[2]], Rd=p[[3]], Tleaf=Tleaf))
  
#   dfr <- cbind(dat, acirun[,2:4])  
    
  l <- list()  
  l$data <- dfr
  l$pars <- p
  l$nlsfit <- nlsfit
  
  class(l) <- "acifit"
  
  
return(l)
}
