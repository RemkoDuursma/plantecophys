#' @export
#' @rdname fitaci
fitacis <- function(data, group,...){
  
  if(!group %in% names(data))
    stop("group variable must be in the dataframe.")
  
  data$group <- data[,group]
  tb <- table(data$group)
  if(any(tb == 0))
    stop("Some levels of your group variable have zero observations.\nUse droplevels() or fix data otherwise!")
  
  d <- split(data, data[,"group"])  
  ng <- length(d)
  success <- vector("logical", ng)
  
  wp <- txtProgressBar(title = "Fittin A-Ci curves", 
                       label = "", min = 0, max = ng, initial = 0, width = 50,style=3)
  
  fits <- list()
  for(i in 1:ng){
    f <- try(fitaci(d[[i]],quiet=TRUE,...), silent=TRUE)
    success[i] <- !inherits(f, "try-error")
    
    fits[[i]] <- if(success[i]) f else NA
    setTxtProgressBar(wp, i)
  }
  close(wp)

  names(fits) <- names(d)
  
  if(any(!success)){
    message("The following groups could not be fit:")
    print(names(d)[!success])
  }
  
  # toss unfitted ones
  fits <- fits[success]

  class(fits) <- "acifits"

return(fits)
}

#' @method plot acifits
#' @S3method plot acifits
#' @param how If 'manyplots', produces a single plot for each A-Ci curve. If 'oneplot' overlays all of them.
#' @rdname fitaci
#' @export
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


#' @method coef acifits
#' @S3method coef acifits
coef.acifits <- function(object,...){
  
  f <- lapply(object, function(x)c(x$pars))
  pars <- as.data.frame(do.call(rbind,f))
  rn <- rownames(object[[1]]$pars)
  nm <- c(rn, paste0(rn,"_SE"))
  names(pars) <- nm
  
return(pars)
}



