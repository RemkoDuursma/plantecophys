
#' Fit A-Ci curve
#' 
#' @description Fits the Farquhar model of photosynthesis to measurements of photosynthesis and intercellular CO2 concentration (Ci). Estimates Jmax, Vcmax, Rd and their (approximate) standard errors. 
#' @details Fits the entire model at once (no splicing)... details to follow. nls2.
#' @export
#' @param dat Dataset with Ci, Photo, Tleaf, PPFD (the last two are optional).
#' @param varnames List of names of variables (see Details).
#' @param nlsmethod Method passed to nls2 (!!not implemented yet)  
fitaci <- function(dat, varnames=list(ALEAF="Photo", Tleaf="Tleaf", Ci="Ci", PPFD="PPFD"),
                   nlsmethod="default"){
  
  if(!varnames$PPFD %in% names(dat)){
    dat$PPFD <- 1800
    warning("PPFD not in dataset; assumed PPFD = 1800.")
  } else dat$PPFD <- dat[,varnames$PPFD]
  
  
  if(!varnames$Tleaf %in% names(dat)){
    dat$Tleaf <- 25
    warning("Tleaf not in dataset; assumed Tleaf = 25.")
  } else dat$Tleaf <- dat[,varnames$Tleaf]
  
  dat$Ci <- dat[,varnames$Ci]
  dat$ALEAF <- dat[,varnames$ALEAF]
  
  acifun_wrap <- function(...){
    r <- Aci(...)
    r$ALEAF
  }
  
  nlsfit <- nls2(Photo ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, Jmax=Jmax, Rd=Rd, Tleaf=Tleaf),
                data=dat, algorithm=nlsmethod,
                start=list(Vcmax=80, Jmax=160, Rd=1))

  p <- coef(nlsfit)
  acirun <- with(dat, Aci(Ci, PPFD=PPFD, Vcmax=p[[1]], Jmax=p[[2]], Rd=p[[3]], Tleaf=Tleaf))
  acirun$Ameas <- dat$ALEAF
  acirun$ELEAF <- NULL
  acirun$GS <- NULL
  acirun$Ca <- NULL
  names(acirun)[names(acirun) == "ALEAF"] <- "Amodel"
  
  # shuffle
  avars <- match(c("Ci","Ameas","Amodel"),names(acirun))
  acirun <- acirun[,c(avars, setdiff(1:ncol(acirun), avars))]
  
  l <- list()  
  l$df <- acirun
  l$pars <- summary(nlsfit)$coefficients[,1:2]
  l$nlsfit <- nlsfit
  
  class(l) <- "acifit"
  
return(l)
}
