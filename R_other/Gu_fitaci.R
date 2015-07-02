fitaci_Gu <- function(data, varnames=list(ALEAF="Photo", Tleaf="Tleaf", Ci="Ci", PPFD="PARi", Rd="Rd"), ...){
  
  data$Ci <- data[,varnames$Ci]
  data$ALEAF <- data[,varnames$ALEAF]
  
  # Transition points to be tested
  # Make sure we have at least 3 on either side.
  # (Note, transition points are placed in between Ci points, fitaci will re-estimate them).
  ci <- sort(data$Ci)
  nci <- length(ci)
  
  trans <- ci[2:(nci-2)] + diff(ci[2:(nci-1)])/2
  ntrans <- length(trans)
  
  # Loop through, fitaci every time and collect results.
  l <- list()
  for(i in 1:ntrans){
    
    r <- fitaci(data, varnames=varnames, citransition=trans[i], ...)
    r$npoint_vcmax <- length(ci[ci < trans[i]])
    r$npoint_jmax <- length(ci[ci > trans[i]])
    r$trans <- trans[i]
    l[[i]] <- r
  }

return(l)
}

f <- fitaci_Gu(c3)
