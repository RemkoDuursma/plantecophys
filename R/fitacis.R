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
plot.acifits <- function(x,...){
  
  for(i in seq_along(x))
    plot.acifit(x[[i]],main=names(x)[i],...)
    
}






