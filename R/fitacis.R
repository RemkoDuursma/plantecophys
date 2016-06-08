#' Fit multiple A-Ci curves at once
#' 
#' @description A convenient function to fit many curves at once, by calling \code{\link{fitaci}} for every group in the dataset. The data provided must include a variable that uniquely identifies each A-Ci curve.
#' 
#' @param data Dataframe with Ci, Photo, Tleaf, PPFD (the last two are optional). For \code{fitacis}, also requires a grouping variable.
#' @param group The name of the grouping variable in the dataframe (an A-Ci curve will be fit for each group separately).
#' @param fitmethod Method to fit the A-Ci curve. Either 'default' (Duursma 2015), or 'bilinear'. See Details.
#' @param progressbar Display a progress bar (default is TRUE).
#' @param quiet If TRUE, no messages are written to the screen.
#' @param x For \code{plot.acifits}, an object returned from \code{fitacis}
#' @param xlim,ylim The X and Y axis limits.
#' @param add If TRUE, adds the plots to a current plot.
#' @param how If 'manyplots', produces a single plot for each A-Ci curve. If 'oneplot' overlays all of them.
#' @param highlight If a name of a curve is given (check names(object), where object is returned by acifits), all curves are plotted in grey, with the highlighted one on top.
#' @param what What to plot, either 'model' (the fitted curve), 'data' or 'none'. See examples.
#' @param object For \code{coef.acifits}, the object returned by \code{fitacis}.
#' @param \dots Further arguments passed to \code{\link{fitaci}} (in the case of \code{fitacis}), or \code{\link{plot.acifit}} (in the case of \code{plot.acifits}).
#' 
#' @examples
#' 
#' # Fit many curves (using an example dataset)
#' # The bilinear method is much faster, but compare using default!
#' fits <- fitacis(manyacidat, "Curve", fitmethod="bilinear")
#' with(coef(fits), plot(Vcmax, Jmax))
#' 
#' # The resulting object is a list, with each component an object as returned by fitaci
#' # So, we can extract one curve:
#' fits[[1]]
#' plot(fits[[1]])
#' 
#' # Plot all curves in separate figures with plot(fits)
#' # Or, in one plot:
#' plot(fits, how="oneplot")
#' 
#' # Note that parameters can be passed to plot.acifit. For example,
#' plot(fits, how="oneplot", what="data", col="blue")
#' plot(fits, how="oneplot", add=TRUE, what="model", lwd=c(1,1))
#' 
#' # Other elements can be summarized with sapply. For example, look at the RMSE:
#' rmses <- sapply(fits, "[[", "RMSE")
#' plot(rmses, type='h', ylab="RMSE", xlab="Curve nr")
#' 
#' # And plot the worst-fitting curve:
#' plot(fits[[which.max(rmses)]])
#' 
#' 
#' @export
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
fitacis <- function(data, group, fitmethod=c("default","bilinear"),
                    progressbar=TRUE, quiet=FALSE, ...){
  
  fitmethod <- match.arg(fitmethod)
  
  if(!group %in% names(data))
    stop("group variable must be in the dataframe.")
  
  if(quiet)progressbar <- FALSE
  
  data$group <- data[,group]
  tb <- table(data$group)
  if(any(tb == 0))
    stop("Some levels of your group variable have zero observations.\nUse droplevels() or fix data otherwise!")
  
  d <- split(data, data[,"group"])  
  ng <- length(d)
  fits <- do_fit_bygroup(d, 1:ng, progressbar, fitmethod, ...)
  
  if(any(!fits$success)){
    if(!quiet){
      group_fail <- names(d)[!fits$success]
      message("The following groups could not be fit:")
      message(paste(group_fail,collapse="\n"))
    }
    
    # Refit bad curves using the 'bilinear' method
    if(fitmethod == "default"){
      if(!quiet)message("Fitting remaining curves with fitmethod='bilinear'.")
      refits <- do_fit_bygroup(d, which(!fits$success), progressbar=FALSE, fitmethod="bilinear", ...)
      fits$fits[!fits$success] <- refits$fits
    }
  }
  
  l <- fits$fits
  class(l) <- "acifits"
  attributes(l)$groupname <- group
  
return(l)
}


do_fit_bygroup <- function(d, which=NULL, progressbar, fitmethod, ...){
  
  ng <- length(d)
  if(is.null(which))which <- 1:ng
  success <- vector("logical", length(which))
  
  if(progressbar){
    wp <- txtProgressBar(title = "Fitting A-Ci curves", 
                         label = "", min = 0, max = ng, initial = 0, 
                         width = 50, style=3)
  }
  
  fits <- list()
  for(i in which){
    f <- try(fitaci(d[[i]], quiet=TRUE, fitmethod=fitmethod, ...), silent=TRUE)
    success[i] <- !inherits(f, "try-error")
    
    fits[[i]] <- if(success[i]) f else NA
    if(progressbar)setTxtProgressBar(wp, i)
  }
  if(progressbar)close(wp)
  
  names(fits) <- names(d)[which]
  
  l <- list(fits=fits, success=success)
}


#' @export plot.acifits
#' @S3method plot acifits
#' @rdname fitacis
plot.acifits <- function(x, how=c("manyplots","oneplot"),
                         highlight=NULL, ylim=NULL,xlim=NULL,
                         add=FALSE, what=c("model","data","none"),
                         ...){
  
  how <- match.arg(how)
  what <- match.arg(what)
  
  if(is.null(ylim)){
    amax <- max(sapply(x, function(x)max(x$df$Amodel)))
    amin <- max(sapply(x, function(x)min(x$df$Amodel)))
    ylim <- c(amin,amax)
  }
  if(is.null(xlim)){
    cimax <- max(sapply(x, function(x)max(x$df$Ci)))
    cimin <- min(sapply(x, function(x)min(x$df$Ci)))
    xlim <- c(cimin,cimax)
  }
  
  if(how == "manyplots"){
    if(add)Warning("Argument 'add' ignored when making multiple plots.")  
    
    for(i in seq_along(x)){
      plot.acifit(x[[i]],main=names(x)[i],xlim=xlim,ylim=ylim,...)
    }
  }
  
  if(how == "oneplot"){
    
    if(!is.null(highlight)){
      if(!highlight %in% names(x))
          stop("Curve ID not found.")
      
      hi <- which(names(x) == highlight)
      
      if(!add){
        plot.acifit(x[[1]], what="none", ylim=ylim, xlim=xlim, whichA="Amin", ...)
      }
      
      for(i in seq_along(x)){
        plot.acifit(x[[i]], what=what, whichA="Amin", add=TRUE,
                    linecols="grey",...)  
      }
      plot.acifit(x[[hi]], what=what, whichA="Amin", add=TRUE,
                  linecols="black",...)  
      
    } else {
      if(!add)plot.acifit(x[[1]], what="none",ylim=ylim, xlim=xlim, whichA="Amin", ...)
      for(i in seq_along(x))
        plot.acifit(x[[i]], what=what, whichA="Amin", add=TRUE,...)  
    }
    
    
    
  }
}


#' @export coef.acifits
#' @S3method coef acifits
#' @rdname fitacis
coef.acifits <- function(object,...){
  
  f <- lapply(object, function(x)c(x$pars))
  pars <- as.data.frame(do.call(rbind,f))
  rn <- rownames(object[[1]]$pars)
  nm <- c(rn, paste0(rn,"_SE"))
  names(pars) <- nm
  
  d <- data.frame(group=names(object))
  names(d) <- attr(object,"group")
  pars <- cbind(d,pars)
  rownames(pars) <- NULL
  
return(pars)
}




