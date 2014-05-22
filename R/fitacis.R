#' @export
#' @rdname fitaci
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
#' @param how If 'manyplots', produces a single plot for each A-Ci curve. If 'oneplot' overlays all of them.
#' @rdname fitaci
plot.acifits <- function(x, how=c("manyplots","oneplot"), ...){
  
  how <- match.arg(how)
  
  if(how == "manyplots"){
  for(i in seq_along(x))
    plot.acifit(x[[i]],main=names(x)[i],...)
  }
  
  if(how == "oneplot"){
    amax <- max(sapply(fits, function(x)max(x$df$Amodel)))
    amin <- max(sapply(fits, function(x)min(x$df$Amodel)))
    
    plot.acifit(x[[1]], what="model",ylim=c(amin,amax), whichA="Amin")
    for(i in seq_along(x))plot.acifit(x[[i]], what="model", whichA="Amin", add=TRUE,...)
  }
}


#' @export coef.acifits
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




