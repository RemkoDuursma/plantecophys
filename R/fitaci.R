
#' Fit A-Ci curve
#' 
#' @description Fits the Farquhar model of photosynthesis to measurements of photosynthesis and intercellular \eqn{CO_2}{CO2} concentration (Ci). Estimates Jmax, Vcmax, Rd and their (approximate) standard errors. 
#' @param dat Dataset with Ci, Photo, Tleaf, PPFD (the last two are optional).
#' @param varnames List of names of variables (see Details).
#' @param nlsmethod Method passed to nls2 ('algorithm' in nls2).
#' @param quiet If TRUE, no messages are written to the screen.
#' @param group For batch analysis using \code{fitacis}, the name of the grouping variable in the dataframe.
#' @details Uses non-linear regression to fit an A-Ci curve. No assumptions are made on which part of the curve is Vcmax or Jmax limited. Three parameters are estimated, Jmax, Vcmax (both at 25deg C) and Rd (at the measurement temperature).
#' 
#' When plotting the fit, the A-Ci curve is simulated using the \code{\link{Aci}} function, with leaf temperature (Tleaf) and PPFD set to the mean value for the dataset. If PPFD is not provided in the dataset, it is assumed to equal 1800 mu mol m-2 s-1.
#' @export
#' @rdname fitaci
fitaci <- function(dat, varnames=list(ALEAF="Photo", Tleaf="Tleaf", Ci="Ci", PPFD="PPFD"),
                   nlsmethod="default", quiet=FALSE){
  
  if(!varnames$PPFD %in% names(dat)){
    dat$PPFD <- 1800
    if(!quiet)warning("PPFD not in dataset; assumed PPFD = 1800.")
  } else dat$PPFD <- dat[,varnames$PPFD]
  
  
  if(!varnames$Tleaf %in% names(dat)){
    dat$Tleaf <- 25
    if(!quiet)warning("Tleaf not in dataset; assumed Tleaf = 25.")
  } else dat$Tleaf <- dat[,varnames$Tleaf]
  
  dat$Ci <- dat[,varnames$Ci]
  dat$ALEAF <- dat[,varnames$ALEAF]
  
  acifun_wrap <- function(...){
    r <- Aci(...)
    r$ALEAF
  }
  
  # Guess some initial values.
  maxCi <- max(dat$Ci)
  maxPhoto <- dat$ALEAF[which.max(dat$Ci)]
  VJ <- maxPhoto / ((maxCi - 45) / (maxCi + 2*45))
  Jmax_guess <- VJ*4
  Vcmax_guess <- Jmax_guess/1.9
  Rd_guess <- 0.015*Vcmax_guess
  
  nlsfit <- nls2(Photo ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, Jmax=Jmax, Rd=Rd, Tleaf=Tleaf),
                data=dat, algorithm=nlsmethod,
                start=list(Vcmax=Vcmax_guess, Jmax=Jmax_guess, Rd=Rd_guess))

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
