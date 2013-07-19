
#' Fit the Farquhar Berry von Caemmerer model of photosynthesis
#' 
#' @description Fits the Farquhar model of photosynthesis to measurements of photosynthesis and intercellular \eqn{CO_2}{CO2} concentration (Ci). Estimates Jmax, Vcmax, Rd and their (approximate) standard errors. Temperature dependencies are taken into account, see \code{\link{Photosyn}}.
#' @param dat Dataset with Ci, Photo, Tleaf, PPFD (the last two are optional).
#' @param varnames List of names of variables (see Details).
#' @param nlsmethod Method passed to nls2 ('algorithm' in nls2).
#' @param quiet If TRUE, no messages are written to the screen.
#' @param group For batch analysis using \code{fitacis}, the name of the grouping variable in the dataframe.
#' @details Uses non-linear regression to fit an A-Ci curve. No assumptions are made on which part of the curve is Vcmax or Jmax limited. Three parameters are estimated, Jmax, Vcmax (both at 25deg C) and Rd (at the measurement temperature).
#' 
#' When plotting the fit, the A-Ci curve is simulated using the \code{\link{Aci}} function, with leaf temperature (Tleaf) and PPFD set to the mean value for the dataset. If PPFD is not provided in the dataset, it is assumed to equal 1800 mu mol m-2 s-1.
#' @return A list of class 'acifit', with three components:
#' 'df' is a dataframe with the original data, and the fitted photosynthetic rate (Amodel). 'pars' contains the parameter estimates and their approximate standard errors. 'nlsfit' is the object returned by \code{\link{nls2}}, and contains more detail on the quality of the fit.
#' @examples
#' # Fit an A-Ci curve on a dataframe that contains Ci, Photo and optionally Tleaf and PPFD. Here, we use the built-in example dataset 'acidata1'.
#' f <- fitaci(acidata1)
#' 
#' # Make a standard plot
#' plot(f)
#' 
#' # Look at a summary of the fit
#' summary(f)
#' 
#' # Extract coefficients only
#' coef(f)
#' 
#' # The object 'f' also contains the original data with predictions.
#' # Here, Amodel are the modelled (fitted) values, Ameas are the measured values.
#' with(f$df, plot(Amodel, Ameas))
#' abline(0,1)
#' 
#' # The fitted values can also be extracted with the fitted() function:
#' fitted(f)
#' 
#' # The non-linear regression (nls) fit is stored as well,
#' summary(f$nlsfit)
#' 
#' @export
#' @rdname fitaci
fitaci <- function(dat, varnames=list(ALEAF="Photo", Tleaf="Tleaf", Ci="Ci", PPFD="PPFD"),
                   nlsmethod="default", quiet=FALSE, ...){
  
  # Set extra parameters if provided
  m <- as.list(match.call())
  a <- as.list(formals(fitaci))
  f <- names(formals(Photosyn))
  
  extrapars <- setdiff(names(m), c(names(a),""))
  for(i in seq_along(extrapars)){
    if(extrapars[i] %in% f)
      formals(Photosyn)[extrapars[i]] <- m[[extrapars[i]]]
    else
      warning("Parameter ",extrapars[i]," not recognized.")
  }
  photpars <- formals(Photosyn)
  removevars <- c("whichA")
  photpars <- photpars[-which(names(photpars) %in% removevars)]
    
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
  
  acifun_wrap <- function(Ci,...){
    r <- Photosyn(Ci=Ci,...)
    r$ALEAF
  }
  
  # Guess some initial values.
  # Here I assume, Jmax/Vcmax = 1.9, GammaStar=45, Rd/Vcmax=0.015.
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
  
  # Save function itself, the formals contain the parameters used to fit the A-Ci curve.
  l$Photosyn <- Photosyn
  
  class(l) <- "acifit"
  
return(l)
}

#' @S3method print acifit
#' @method print acifit
print.acifit <- function(x,...){
  
  cat("Result of fitaci.\n\n")
  
  cat("Data and predictions:\n")
  print(x$df)
  
  cat("\nEstimated parameters:\n")
  
  print(x$pars)
  cat("Note: Vcmax, Jmax are at 25C, Rd is at measurement T.")
  cat("\n\n")
  
  cat("Parameter settings:\n")
  fm <- formals(x$Photosyn)
  pars <- c("alpha","theta","EaV","EdVC","delsC","EaJ","EdVJ","delsJ")
  fm <- fm[pars]
  cat(paste0(names(fm)," = ", unlist(fm),"\n"))
  
}

#' @S3method summary acifit
#' @method summary acifit
summary.acifit <- function(x,...){
  
  print.acifit(x, ...)
  
}


#' @S3method coef acifit
#' @method coef acifit
coef.acifit <- function(object, ...){
 v <- unname(object$pars[,1])
 names(v) <- rownames(object$pars)
return(v)
}

#' @S3method fitted acifit
#' @method fitted acifit
fitted.acifit <- function(object,...){
  
  object$df$Amodel
  
}




