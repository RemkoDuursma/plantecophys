DewtoVPD <- function(Tdew, TdegC, Pa=101){

	# Actual vapor pressure.
	e <- esat(Tdew)
	
	# saturated:
	esatval <- esat(TdegC)
	
	return((esatval - e)/1000)
}