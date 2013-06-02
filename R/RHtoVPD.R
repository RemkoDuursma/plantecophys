#' Conversions between relative humidity and vapour pressure deficit
#' 
#' @description Converts relative humidity (%) to vapour pressure deficit (kPa), or vice versa. To convert from relative humidity, use the \code{RHtoVPD} function, use \code{VPDtoRH} for the other way around. The water vapor saturation pressure is calculated with \code{esat}, which is from Jones 1992.
#' @param RH Relative humidity (%)
#' @param TdegC Temperature (degrees C)
#' @param VPD Vapour pressure deficit (kPa)
#' @param Pa Atmospheric pressure
#' @export
#' @aliases RHtoVPD VPDtoRH esat
#' @references Jones, H.G. 1992. Plants and microclimate: a quantitative approach to environmental plant physiology. 2nd Edition., 2nd Edn. Cambridge University Press, Cambridge. 428 p.
#' @author Remko Duursma
RHtoVPD <- function(RH, TdegC, Pa=101){
	esatval <- esat(TdegC)
	e <- (RH/100) * esatval
	VPD <- (esatval - e)/1000
return(VPD)
}

VPDtoRH <- function(VPD, TdegC, Pa=101){
  esatval <- esat(TdegC)
  e <- pmax(0, esatval - VPD*1000)
  RH <- 100 * e/esatval
  return(RH)
}

esat <- function(TdegC, Pa=101){  
  a <- 611.21
  b <- 17.502
  c <- 240.97
  f <- 1.0007 + 3.46 * 10^-8 * Pa * 1000
  esatval <- f * a * (exp(b * TdegC/(c + TdegC)))
  return(esatval)
}