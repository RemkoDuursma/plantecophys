#' Convert relative humidity to vapour pressure deficit
#' 
#' @description Converts relative humidity (%) to vapour pressure deficit (kPa) using the formula from Jones 1992.
#' @param RH Relative humidity (%)
#' @param TdegC Temperature (degrees C)
#' @param Pa Atmospheric pressure
#' @export
#' @references Jones, H.G. 1992. Plants and microclimate: a quantitative approach to environmental plant physiology. 2nd Edition., 2nd Edn. Cambridge University Press, Cambridge. 428 p.
RHtoVPD <- function(RH, TdegC, Pa=101){
	esatval <- esat(TdegC)
	e <- (RH/100) * esatval
	VPD <- (esatval - e)/1000
return(VPD)
}
