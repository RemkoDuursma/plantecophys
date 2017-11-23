#' Calculate transition points for fitted A-Ci curves
#'@description Calculates the Ci at the transition points between Ac & Aj (point 1), and Aj and Ap (point 2). The latter is not NA only when TPU was estimated (and estimable), see \code{\link{fitaci}}, argument \code{fitTPU}.
#'@param object Either an object returned by \code{fitaci}, or a copy of the \code{Photosyn} function.
#'@param \dots Further arguments passed to the \code{Photosyn} function.
#'@details This function is also used by \code{fitaci}, the results are stored in elements \code{Ci_transition} and \code{Ci_transition2}.
#'@export
findCiTransition <- function(object, ...){
  
  if(is.function(object)){
    photofun <- object
  } else {
    photofun <- object$Photosyn
  }
  
  O <- function(Ci, photofun, point, ...){
    x <- photofun(Ci=Ci, ...)
    if(point == 1){
      Ac <- x$Ac - x$Rd
      Aj <- x$Aj - x$Rd
      return((Ac - Aj)^2)
    } 
    if(point == 2){
      Aj <- x$Aj - x$Rd
      Ap <- x$Ap - x$Rd
      return((Aj - Ap)^2)
    }
  }
  p1 <- optimize(O, c(25,2000), photofun=photofun, point=1, ...)$minimum
  
  r <- photofun(Ci=500, ...)
  
  if(r$Ap >= 1000 | is.na(r$Ap)){
    p2 <- NA
  } else {
    p2 <- optimize(O, c(400,2000), photofun=photofun, point=2, ...)$minimum
  }
  
  return(c(Ac_Aj=p1, Aj_Ap=p2))
}
