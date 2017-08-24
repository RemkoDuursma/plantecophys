#' Conversions between relative humidity, vapour pressure deficit and dewpoint
#' 
#' @description A collection of functions to convert between relative humidity (RH) (\%), 
#' vapour pressure deficit (VPD) (kPa), 
#' dew point temperature, and leaf- or air temperature-based VPD or RH. To convert from 
#' relative humidity to VPD, 
#' use the \code{RHtoVPD} function,  use \code{VPDtoRH} for the other way around. The water 
#' vapor saturation pressure is 
#' calculated with \code{esat}. Use \code{DewtoVPD} to 
#' convert from dewpoint temperature to VPD. The functions \code{VPDleafToAir} and \code{VPDairToLeaf}
#' convert VPD from a leaf temperature to an air-temperature basis and vice versa. The 
#' functions \code{RHleafToAir} a \code{RHairToLeaf} do the same for relative humidity.
#' @details The function describing saturated vapor pressure with temperature is taken from 
#' Jones (1992). All other calculations follow directly from the standard definitions, for
#'  which Jones (1992) may also be consulted.
#' @param RH Relative humidity (\%)
#' @param TdegC Temperature (degrees C) (either leaf or air)
#' @param Tair Air temperature (degrees C)
#' @param Tleaf Leaf temperature (degrees C)
#' @param VPD Vapour pressure deficit (kPa)
#' @param Pa Atmospheric pressure (kPa)
#' @param Tdew Dewpoint temperature (degrees C)
#' @export RHtoVPD VPDtoRH esat DewtoVPD VPDtoDew VPDleafToAir VPDairToLeaf
#' @rdname Conversions
#' @references Jones, H.G. 1992. Plants and microclimate: a quantitative approach to 
#' environmental plant physiology. 2nd Edition., 2nd Edn. Cambridge University Press, Cambridge. 428 p.
#' @author Remko Duursma
RHtoVPD <- function(RH, TdegC, Pa=101){
	esatval <- esat(TdegC, Pa)
	e <- (RH/100) * esatval
	VPD <- (esatval - e)/1000
return(VPD)
}
#' @rdname Conversions
VPDtoRH <- function(VPD, TdegC, Pa=101){
  esatval <- esat(TdegC, Pa)
  e <- pmax(0, esatval - VPD*1000)
  RH <- 100 * e/esatval
  return(RH)
}
#' @rdname Conversions
esat <- function(TdegC, Pa=101){  
  a <- 611.21
  b <- 17.502
  c <- 240.97
  f <- 1.0007 + 3.46 * 10^-8 * Pa * 1000
  esatval <- f * a * (exp(b * TdegC/(c + TdegC)))
  return(esatval)
}

# inverse of esat (calc T given a saturation vapor pressure)
T_esat <- function(sat, Pa=101){
  a <- 611.21
  b <- 17.502
  c <- 240.97
  f <- 1.0007 + 3.46 * 10^-8 * Pa * 1000
  
  phi <- log(sat/(f*a))
  (c*phi)/(b-phi)
}
#' @rdname Conversions
VPDtoDew <- function(VPD, TdegC, Pa=101){
  
  esatval <- esat(TdegC, Pa)
  e <- pmax(0, esatval - VPD*1000)
  T_esat(e, Pa)
}

#' @rdname Conversions
DewtoVPD <- function(Tdew, TdegC, Pa=101){
  
  # Actual vapor pressure.
  e <- esat(Tdew, Pa)
  
  # saturated:
  esatval <- esat(TdegC)
  
  return((esatval - e)/1000)
}
#' @rdname Conversions
VPDleafToAir <- function(VPD, Tleaf, Tair, Pa=101){
  
  e <- esat(Tleaf, Pa) - VPD*1000
  vpd <- esat(Tair, Pa) - e
  
  return(vpd/1000)
}
#' @rdname Conversions
VPDairToLeaf <- function(VPD, Tair, Tleaf, Pa=101){
  
  e <- esat(Tair, Pa) - VPD*1000
  vpd <- esat(Tleaf, Pa) - e
  
  return(vpd/1000)
}
#' @rdname Conversions
#' @export
RHleafToAir <- function(RH, Tleaf, Tair, Pa=101){
  
  e <- (RH/100)*esat(Tleaf, Pa)
  rh <- e/esat(Tair, Pa)
  
  return(rh*100)
}
#' @rdname Conversions
#' @export
RHairToLeaf <- function(RH, Tair, Tleaf, Pa=101){
  
  e <- (RH/100)*esat(Tair, Pa)
  rh <- e/esat(Tleaf, Pa)
  
  return(rh*100)
}


