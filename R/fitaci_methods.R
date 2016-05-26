

#' @export print.acifit
#' @S3method print acifit
#' @rdname fitaci
print.acifit <- function(x,...){
  
  cat("Result of fitaci.\n\n")
  
  cat("Data and predictions:\n")
  print(x$df)
  
  cat("\nEstimated parameters:\n")
  
  print(x$pars)
  if(x$Tcorrect)
    cat("Note: Vcmax, Jmax are at 25C, Rd is at measurement T.\n")
  else
    cat("Note: Vcmax, Jmax, Rd are at measurement T.\n")
  
  if(!is.na(x$gmeso)){
    cat("Note: Mesophyll conductance was input, Vcmax and Jmax are Cc-based rates.\n")
  }
  
  if(x$Rd_measured)
    cat("Note: measured Rd was provided, only Vcmax and Jmax were fit.\n")
  
  cat("\nCurve was fit using method: ", x$fitmethod, "\n")
  
  if(!is.na(x$citransition)){
    cat("\nCi transition was constrained to be: ", x$citransition, "\n")
    cat("Actual fitted Ci transition: ", x$Ci_transition,"\n")
  }

  cat("\nParameter settings:\n")
  fm <- formals(x$Photosyn)
  pars <- c("Patm","alpha","theta","EaV","EdVC","delsC","EaJ","EdVJ","delsJ")
  fm <- unlist(fm[pars])
  cat(paste0(pars," = ", fm,"\n"))
  
  if(!x$gstarinput | !x$kminput){
    cat("\nEstimated from Tleaf (shown at mean Tleaf):\n")
    if(!x$gstarinput)cat("GammaStar = ",x$GammaStar,"\n")
    if(!x$kminput)cat("Km = ",x$Km,"\n")
  }
  
  if(x$gstarinput | x$kminput){
    cat("\nSet by user:\n")
    if(x$gstarinput)cat("GammaStar = ",x$GammaStar,"\n")
    if(x$kminput)cat("Km = ",x$Km,"\n")
  }
  
}

#' @export summary.acifit
#' @S3method summary acifit
#' @rdname fitaci
summary.acifit <- function(object,...){
  
  print.acifit(object, ...)
  
}


#' @export coef.acifit
#' @S3method coef acifit
#' @rdname fitaci
coef.acifit <- function(object, ...){
  v <- unname(object$pars[,1])
  names(v) <- rownames(object$pars)
  return(v)
}

#' @export fitted.acifit
#' @S3method fitted acifit
#' @rdname fitaci
fitted.acifit <- function(object,...){
  
  object$df$Amodel
  
}


#' @export plot.acifit
#' @S3method plot acifit
#' @param x For plot.acifit, an object returned by \code{fitaci}
#' @param xlim Limits for the X axis, if left blank estimated from data
#' @param ylim Limits for the Y axis, if left blank estimated from data
#' @param whichA By default all three photosynthetic rates are plotted (Aj=Jmax-limited (blue), Ac=Vcmax-limited (red), Hyperbolic minimum (black)). Or, specify one or two of them. 
#' @param what The default is to plot both the data and the model fit, or specify 'data' or 'model' to plot one of them, or 'none' for neither (only the plot region is set up)
#' @param add If TRUE, adds to the current plot
#' @param pch The plotting symbol for the data
#' @param addzeroline If TRUE, the default, adds a dashed line at y=0
#' @param addlegend If TRUE, adds a legend (by default does not add a legend if add=TRUE)
#' @param legendbty Box type for the legend, passed to argument bty in \code{\link{legend}}.
#' @param transitionpoint For plot.acifit, whether to plot a symbol at the transition point.
#' @param linecols Vector of three colours for the lines (limiting rate, Ac, Aj), if one value provided it is used for all three.
#' @param lwd Line widths, can be a vector of length 2 (first element for Ac and Aj, second one for the limiting rate).
#' @rdname fitaci
#' @importFrom graphics points
#' @importFrom graphics abline
#' @importFrom graphics legend
plot.acifit <- function(x, what=c("data","model","none"), xlim=NULL, ylim=NULL, 
                        whichA=c("Ac","Aj","Amin","Ap"), add=FALSE, pch=19, 
                        addzeroline=TRUE, addlegend=!add, legendbty='o',
                        transitionpoint=TRUE, linecols=c("black","blue","red"),
                        lwd=c(1,2),
                        ...){
  
  # Note that Ci on the X-axis is in molar units!
  if(is.null(ylim))ylim <- with(x$df, c(min(Ameas), 1.1*max(Ameas)))
  if(is.null(xlim))xlim <- with(x$df,c(0, max(Ci_original)))
  if(length(lwd)==1)lwd <- c(lwd,lwd)
  if(length(linecols)==1)linecols <- rep(linecols,3)
  
  # Vector of Ci values at which to evaluate fitted ACi curve.
  Ci <- with(x$df, seq(min(Ci_original), max(Ci_original), length=101))
  
  # Exact model used to fit the A-Ci curve was saved in the object.
  # (parameter settings etc. are preserved)
  pcor <- mean(x$df$Patm)/100
  pred <- x$Photosyn(Ci=Ci * pcor)
  pred$Ci_original <- pred$Ci / pcor
  
  # Is there a TPU limitation?
  TPUlimit <- any(pred$Ap < pred$Aj)
  
  if(!add){
    with(x$df, plot(Ci_original, Ameas, type='n',
                    ylim=ylim,
                    xlim=xlim,
                    xlab=expression(italic(C)[i]~~(ppm)),
                    ylab=expression(italic(A)[n]~~(mu*mol~m^-2~s^-1)),
                    ...
    ))
  }
  if("data" %in% what)with(x$df, points(Ci_original, Ameas, pch=pch,...))
  
  if("model" %in% what){
    if("Aj" %in% whichA)with(pred, lines(Ci_original, Aj-Rd, col=linecols[2],lwd=lwd[1]))
    if("Ac" %in% whichA)with(pred, lines(Ci_original, Ac-Rd, col=linecols[3],lwd=lwd[1]))
    if("Ap" %in% whichA & TPUlimit)with(pred, lines(Ci_original, Ap-Rd, col="darkgrey", lty=5, lwd=lwd[1]))
    if("Amin" %in% whichA)with(pred, lines(Ci_original, ALEAF, col=linecols[1], lwd=lwd[2]))
  }
  
  if(transitionpoint && "model" %in% what)
    points(x$Ci_transition / pcor, x$Photosyn(Ci=x$Ci_transition)$ALEAF, pch=21, bg="lightgrey", cex=0.8)
  
  if(addzeroline)
    abline(h=0, lty=3)
  
  if(addlegend & ! TPUlimit){
    legend("bottomright", c(expression(italic(A)[c]),
                            expression(italic(A)[j]),
                            "Limiting rate"), lty=1, lwd=c(lwd[1],lwd[1],lwd[2]), 
           col=linecols[3:1], bty=legendbty, bg="white")
  }
  if(addlegend & TPUlimit){
    legend("bottomright", c(expression(italic(A)[c]),
                            expression(italic(A)[j]),
                            expression(italic(A)[p]),
                            "Limiting rate"), lty=1, lwd=c(rep(lwd[1],3), lwd[2]), 
           col=c(linecols[3:2],"darkgrey",linecols[1]), bty=legendbty, bg="white")
  }
  
}
