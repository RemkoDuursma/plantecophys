
#' Fit the Farquhar Berry von Caemmerer model of photosynthesis
#' 
#' @description Fits the Farquhar model of photosynthesis to measurements of photosynthesis and intercellular \eqn{CO_2}{CO2} concentration (Ci). Estimates Jmax, Vcmax, Rd and their (approximate) standard errors. Temperature dependencies are taken into account, see \code{\link{Photosyn}}.
#' @param dat Dataset with Ci, Photo, Tleaf, PPFD (the last two are optional).
#' @param varnames List of names of variables (see Details).
#' @param Tcorrect If TRUE, Vcmax and Jmax are corrected to 25C. Otherwise, Vcmax and Jmax are estimated at measurement temperature.
#' @param quiet If TRUE, no messages are written to the screen.
#' @param group For batch analysis using \code{fitacis}, the name of the grouping variable in the dataframe.
#' @details Uses non-linear regression to fit an A-Ci curve. No assumptions are made on which part of the curve is Vcmax or Jmax limited. Three parameters are estimated, Jmax, Vcmax (both at 25deg C) and Rd (at the measurement temperature).
#' 
#' When plotting the fit, the A-Ci curve is simulated using the \code{\link{Aci}} function, with leaf temperature (Tleaf) and PPFD set to the mean value for the dataset. If PPFD is not provided in the dataset, it is assumed to equal 1800 mu mol m-2 s-1.
#' @return A list of class 'acifit', with five components:
#' \describe{
#' \item{df}{A dataframe with the original data, the fitted photosynthetic rate (Amodel), Jmax and Vcmax-limited gross rates (Aj, Ac)}
#' \item{pars}{Contains the parameter estimates and their approximate standard errors}
#' \item{nlsfit}{The object returned by \code{\link{nls2}}, and contains more detail on the quality of the fit}
#' \item{Photosyn}{A copy of the \code{\link{Photosyn}} function with the arguments adjusted for the current fit. That is, Vcmax, Jmax and Rd are set to those estimated in the fit, and Tleaf and PPFD are set to the mean value in the dataset.}
#' \item{Ci_transition}{The Ci at which photosynthesis transitions from Vcmax to Jmax limited photosynthesis.}
#' }
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
#' # The curve generator is stored as f$Photosyn:
#' # Calculate photosynthesis at some value for Ci, using estimated parameters and mean Tleaf, PPFD for the dataset
#' f$Photosyn(Ci=820)
#' 
#' # Photosynthetic rate at the transition point:
#' f$Photosyn(Ci=f$Ci_transition)$ALEAF
#' 
#' @export
#' @rdname fitaci
fitaci <- function(dat, varnames=list(ALEAF="Photo", Tleaf="Tleaf", Ci="Ci", PPFD="PARi"),
                   Tcorrect=TRUE, quiet=FALSE, ...){
  
  # Set extra parameters if provided
  m <- as.list(match.call())
  a <- as.list(formals(fitaci))
  f <- names(formals(Photosyn))
  
  extrapars <- setdiff(names(m), c(names(a),""))
  for(i in seq_along(extrapars)){
    if(extrapars[i] %in% f){
      val <- m[[extrapars[i]]]
      formals(Photosyn)[extrapars[i]] <- val
    } else {
      warning("Parameter ",extrapars[i]," not recognized.")
    }
  }
  photpars <- formals(Photosyn)
  removevars <- c("whichA")
  photpars <- photpars[-which(names(photpars) %in% removevars)]
  
  # Check if PAR is provided
  if(!varnames$PPFD %in% names(dat)){
    dat$PPFD <- 1800
    if(!quiet)warning("PARi not in dataset; assumed PARi = 1800.")
  } else dat$PPFD <- dat[,varnames$PPFD]
  
  # Check if Tleaf is provided
  if(!varnames$Tleaf %in% names(dat)){
    dat$Tleaf <- 25
    if(!quiet)warning("Tleaf not in dataset; assumed Tleaf = 25.")
  } else dat$Tleaf <- dat[,varnames$Tleaf]
  
  dat$Ci <- dat[,varnames$Ci]
  dat$ALEAF <- dat[,varnames$ALEAF]
  
  
  # Needed to avoid apparent recursion below.
  TcorrectVJ <- Tcorrect
  
  # Wrapper around Photosyn; this wrapper will be sent to nls2. 
  acifun_wrap <- function(Ci,...){
    r <- Photosyn(Ci=Ci,Tcorrect=TcorrectVJ,...)
    r$ALEAF
  }
  
#   # Guess some initial values.
#   # Here I assume, Jmax/Vcmax = 1.9, GammaStar=45, Rd/Vcmax=0.015.
#   maxCi <- max(dat$Ci)
#   maxPhoto <- dat$ALEAF[which.max(dat$Ci)]
#   VJ <- maxPhoto / ((maxCi - 45) / (maxCi + 2*45))
#   Jmax_guess <- VJ*4
#   Vcmax_guess <- Jmax_guess/1.9
#   Rd_guess <- 0.015*Vcmax_guess
#  
#   
#   nlsfit <- nls2(Photo ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, Jmax=Jmax, Rd=Rd, Tleaf=Tleaf),
#                 data=dat, algorithm=nlsmethod,
#                 start=list(Vcmax=Vcmax_guess, Jmax=Jmax_guess, Rd=Rd_guess))
# 
#   browser()
#   
#   # Pre-fit ; this finds the best starting values
#   n <- 50
#   nlsfit_pre <- nls2(Photo ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, Jmax=Jmax, Rd=Rd, Tleaf=Tleaf),
#                  data=dat, algorithm="plinear-brute",
#                  start=data.frame(Vcmax=seq(5,350,length=n),
#                                    Jmax=seq(5,350,length=n),
#                                    Rd=seq(0.01, 11, length=n)))
#   
  
  aciSS <- function(Vcmax, Jmax, Rd){
    Photo_mod <- acifun_wrap(dat$Ci, PPFD=dat$PPFD, 
                             Vcmax=Vcmax, Jmax=Jmax, 
                             Rd=Rd, Tleaf=dat$Tleaf)
    SS <- sum((dat$Photo - Photo_mod)^2)
  return(SS)
  }
  
  # Guess some initial values.
  # Here I assume, Jmax/Vcmax = 1.9, GammaStar=45, Rd/Vcmax=0.015.
  maxCi <- max(dat$Ci)
  mi <- which.max(dat$Ci)
  maxPhoto <- dat$ALEAF[mi]
  Tl <- dat$Tleaf[mi]
  gammastar <- TGammaStar(Tl)
  
  VJ <- maxPhoto / ((maxCi - gammastar) / (maxCi + 2*gammastar))
  Jmax_guess <- VJ*4
  
  dato <- dat[dat$Ci < 150 & dat$Ci > 60,]
  if(nrow(dato) > 0){
    Km <- TKm(dato$Tleaf)
    gammastar <- TGammaStar(dato$Tleaf)
    vcmax <- with(dato, Photo / ((Ci - gammastar)/(Ci + Km)))
    Vcmax_guess <- mean(vcmax)
  } else {
    Vcmax_guess <- Jmax_guess/1.8 
  }
  
  Rd_guess <- 0.03*Vcmax_guess
  
  d <- 0.3
  n <- 20
  gg <- expand.grid(Vcmax=seq(Vcmax_guess*(1-d),Vcmax_guess*(1+d),length=n),
                    #Jmax=seq(Jmax_guess*(1-d),Jmax_guess*(1+d),length=n),
                    Rd=seq(Rd_guess*(1-d),Rd_guess*(1+d),length=n))

  m <- with(gg, mapply(aciSS, Vcmax=Vcmax, Jmax=Jmax_guess, Rd=Rd))
  ii <- which.min(m)
  browser()
  # Now fit with optimized starting values
  nlsfit <- nls(Photo ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, 
                                    Jmax=Jmax, Rd=Rd, Tleaf=Tleaf),
                  data=dat, control=nls.control(maxiter=500),
                  start=list(Vcmax=gg$Vcmax[ii], Jmax=Jmax_guess, Rd=gg$Rd[ii]))
  
  # Using fitted coefficients, get predictions from model.
  p <- coef(nlsfit)
  acirun <- Photosyn(Ci=dat$Ci, 
                     Vcmax=p[[1]], Jmax=p[[2]], Rd=p[[3]], 
                     PPFD=dat$PPFD, 
                     Tleaf=dat$Tleaf,
                     Tcorrect=Tcorrect)
  acirun$Ameas <- dat$ALEAF
  acirun$ELEAF <- NULL
  acirun$GS <- NULL
  acirun$Ca <- NULL
  names(acirun)[names(acirun) == "ALEAF"] <- "Amodel"
  
  # shuffle
  avars <- match(c("Ci","Ameas","Amodel"),names(acirun))
  acirun <- acirun[,c(avars, setdiff(1:ncol(acirun), avars))]
  
  # Organize output
  l <- list()  
  l$df <- acirun[order(acirun$Ci),]
  l$pars <- summary(nlsfit)$coefficients[,1:2]
  l$nlsfit <- nlsfit
  l$Tcorrect <- Tcorrect
  
  # Save function itself, the formals contain the parameters used to fit the A-Ci curve.
  # First save Tleaf, PPFD in the formals (as the mean of the dataset)
  formals(Photosyn)$Tleaf <- mean(dat$Tleaf)
  formals(Photosyn)$PPFD <- mean(dat$PPFD)
  formals(Photosyn)$Vcmax <- l$pars[1]
  formals(Photosyn)$Jmax <- l$pars[2]
  formals(Photosyn)$Rd <- l$pars[3]
  l$Photosyn <- Photosyn
  
  # Store Ci at which photosynthesis transitions from Jmax to Vcmax limitation
  l$Ci_transition <- findCiTransition(l$Photosyn)
  
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
  if(x$Tcorrect)
    cat("Note: Vcmax, Jmax are at 25C, Rd is at measurement T.")
  else
    cat("Note: Vcmax, Jmax, Rd are at measurement T.")
  
  cat("\n\n")
  
  cat("Parameter settings:\n")
  fm <- formals(x$Photosyn)
  pars <- c("alpha","theta","EaV","EdVC","delsC","EaJ","EdVJ","delsJ")
  fm <- fm[pars]
  cat(paste0(names(fm)," = ", unlist(fm),"\n"))
  
}

#' @S3method summary acifit
#' @method summary acifit
summary.acifit <- function(object,...){
  
  print.acifit(object, ...)
  
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



#' @S3method plot acifit
#' @method plot acifit
#' @param what The default is to plot both the data and the model fit, or specify 'data' or 'model' to plot one of them.
#' @param add If TRUE, adds to the current plot
#' @param pch The plotting symbol for the data
#' @param xlim Limits for the X axis, if left blank estimated from data
#' @param ylim Limits for the Y axis, if left blank estimated from data
#' @param whichA By default all three photosynthetic rates are plotted (Aj=Jmax-limited (blue), Ac=Vcmax-limited (red), Hyperbolic minimum (black)). Or, specify one or two of them. 
#' @param addzeroline If TRUE, the default, adds a dashed line at y=0
#' @param addlegend If TRUE, adds a legend (by default does not add a legend if add=TRUE)
#' @rdname fitaci
plot.acifit <- function(x, what=c("data","model"), xlim=NULL, ylim=NULL, whichA=c("Ac","Aj","Amin"), add=FALSE, pch=19, addzeroline=TRUE, addlegend=!add, transitionpoint=TRUE, ...){
  
  if(is.null(ylim))ylim <- with(x$df, c(min(Ameas), 1.1*max(Ameas)))
  if(is.null(xlim))xlim <- with(x$df,c(0, max(Ci)))
  
  Ci <- with(x$df, seq(min(Ci), max(Ci), length=101))
  
  # Exact model used to fit the A-Ci curve was saved in the object.
  pred <- x$Photosyn(Ci=Ci)
  
  if(!add){
    with(x$df, plot(Ci, Ameas, type='n',
                    ylim=ylim,
                    xlim=xlim,
                    xlab=expression(italic(C)[i]~~(ppm)),
                    ylab=expression(italic(A)[net]~~(mu*mol~m^-2~s^-1)),
                    ...
    ))
  }
  if("data" %in% what)with(x$df, points(Ci, Ameas, pch=pch,...))
  
  if("model" %in% what){
    if("Aj" %in% whichA)with(pred, points(Ci, Aj-Rd, type='l', col="blue"))
    if("Ac" %in% whichA)with(pred, points(Ci, Ac-Rd, type='l', col="red"))
    if("Amin" %in% whichA)with(pred, points(Ci, ALEAF, type='l', col="black", lwd=2))
  }
  
  if(transitionpoint)
    points(x$Ci_transition, x$Photosyn(Ci=x$Ci_transition)$ALEAF, pch=21, bg="lightgrey", cex=0.8)
  
  if(addzeroline)
    abline(h=0, lty=3)
  
  if(addlegend){
    legend("bottomright", c(expression(italic(A)[c]),
                            expression(italic(A)[j]),
                            "Limiting rate"), lty=1, lwd=c(1,1,2), col=c("red","blue","black"))
  }
  
}
