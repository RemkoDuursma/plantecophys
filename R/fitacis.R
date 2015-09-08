#' @export
#' @rdname fitaci
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
fitacis <- function(data, group, progressbar=TRUE, quiet=FALSE, ...){
  
  if(!group %in% names(data))
    stop("group variable must be in the dataframe.")
  
  if(quiet)progressbar <- FALSE
  
  data$group <- data[,group]
  tb <- table(data$group)
  if(any(tb == 0))
    stop("Some levels of your group variable have zero observations.\nUse droplevels() or fix data otherwise!")
  
  d <- split(data, data[,"group"])  
  ng <- length(d)
  success <- vector("logical", ng)
  
  if(progressbar){
    wp <- txtProgressBar(title = "Fitting A-Ci curves", 
                       label = "", min = 0, max = ng, initial = 0, 
                       width = 50, style=3)
  }
  
  fits <- list()
  for(i in 1:ng){
    f <- try(fitaci(d[[i]], quiet=TRUE, ...), silent=TRUE)
    success[i] <- !inherits(f, "try-error")
    
    fits[[i]] <- if(success[i]) f else NA
    if(progressbar)setTxtProgressBar(wp, i)
  }
  if(progressbar)close(wp)

  names(fits) <- names(d)
  
  if(any(!success)){
    if(!quiet){
      message("The following groups could not be fit:")
      print(names(d)[!success])
    }
  }
  
  # toss unfitted ones
  fits <- fits[success]

  class(fits) <- "acifits"
  attributes(fits)$groupname <- group
  
return(fits)
}


#' @export plot.acifits
#' @S3method plot acifits
#' @param how If 'manyplots', produces a single plot for each A-Ci curve. If 'oneplot' overlays all of them.
#' @param highlight If a name of a curve is given (check names(object), where object is returned by acifits), all curves are plotted in grey, with the highlighted one on top.
#' @rdname fitaci
plot.acifits <- function(x, how=c("manyplots","oneplot"),
                         highlight=NULL,
                         ...){
  
  how <- match.arg(how)
  
  if(how == "manyplots"){
  for(i in seq_along(x))
    plot.acifit(x[[i]],main=names(x)[i],...)
  }
  
  if(how == "oneplot"){
    
    amax <- max(sapply(x, function(x)max(x$df$Amodel)))
    amin <- max(sapply(x, function(x)min(x$df$Amodel)))
    
    
    if(!is.null(highlight)){
      if(!highlight %in% names(x))
          stop("Curve ID not found.")
      
      hi <- which(names(x) == highlight)
      
      plot.acifit(x[[1]], what="none",ylim=c(amin,amax), whichA="Amin", ...)
      for(i in seq_along(x)){
        plot.acifit(x[[i]], what="model", whichA="Amin", add=TRUE,
                    linecols="grey",...)  
      }
      plot.acifit(x[[hi]], what="model", whichA="Amin", add=TRUE,
                  linecols="black",...)  
      
    } else {
      plot.acifit(x[[1]], what="none",ylim=c(amin,amax), whichA="Amin", ...)
      for(i in seq_along(x))
        plot.acifit(x[[i]], what="model", whichA="Amin", add=TRUE,...)  
    }
    
    
    
  }
}


#' @export coef.acifits
#' @S3method coef acifits
#' @rdname fitaci
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




