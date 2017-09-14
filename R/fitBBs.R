#' Fit Ball-Berry type models of stomatal conductance to many groups at once
#'@description A batch utility for the \code{\link{fitBB}} function, to fit the model for each group in a dataframe. 
#'@param data Input dataframe, containing all variables needed to fit the model.
#'@param group Name of the grouping variable in the dataframe (quoted), the model will be fit for each group defined by this variable.
#'@param \dots Further parameters passed to \code{\link{fitBB}}, see there for a full description.
#'@export
#'@examples
#'\dontrun{
#'# If you have a factor variable in your dataset called 'species', and you
#'# want to fit the Ball-Berry model for each of the species:
#'myfits <- fitBBs(mydata, "species", model="BallBerry")
#'
#'# A dataframe with coefficients is returned by coef()
#'coef(myfits)
#'
#'}
fitBBs <- function(data, group, ...){
  
  
  datasp <- split(data, data[,group])
  
  fits <- lapply(datasp, fitBB, ...)
  
  class(fits) <- "BBfits"
  
return(fits)
}

#' @method coef BBfits
#' @export
coef.BBfits <- function(object, ...){
  
  p <- do.call(rbind,lapply(object, "[[", "coef"))
  dfrout <- cbind(data.frame(group=rownames(p)), p)
  rownames(dfrout) <- NULL
  
return(dfrout)
}

#' @method print BBfits
#' @export
print.BBfits <- function(x, ...){
  
  cat("Result of fitBBs.\n")
  cat("Fitted", x[[1]]$gsmodel, "model to", length(x), "groups\n")
  
  if(x[[1]]$fitg0){
    cat("Both g0 and g1 were estimated.\n\n")
  } else {
    cat("Only g1 was estimated (g0 = 0).\n\n")
  }
  
  cat("To return dataframe with coefficients, do: coef(myfit).\n")
  cat("(where myfit is the name of the object returned by fitBBs)\n")
}
