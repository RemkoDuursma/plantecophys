

acifun_uv <- function(Ci, PPFD=1500, Tleaf=25, alpha=0.24, theta=0.85, 
                      Jmax=100, Vcmax=50, Rd=1, ...){
  
  
  Jfun <- function(PPFD, alpha, Jmax, theta){
    (alpha*PPFD + Jmax - sqrt((alpha*PPFD + Jmax)^2 - 4*alpha*theta*PPFD*Jmax))/(2*theta)
  }
  
  Km <- exp(38.05-79430/(8.314*(Tleaf+273)))*(1+210/exp(20.3-36380/(8.314*(Tleaf+273))))
  GammaStar <- 42.75*exp(37830*(Tleaf-25)/(8.314*298.15*(Tleaf+273.15)))
  
  Rgas <- 8.314
  
  Ea <- 82620.87
  EaJ <- 39676.89
  EdJ <- 200000
  delsJ <- 641.3615
  delsC <- 645.1013
  
  # Vcmax <- Vcmax*exp(Ea/Rgas*(1/298.15 - 1/(Tleaf+273.15)))  
  Vcmax <- Vcmax * exp((Tleaf-25)*Ea/(Rgas*(Tleaf+273.15)*(25 + 273.15))) * 
    (1+exp((delsC*(25 + 273.15)-EdJ)/(Rgas*(25 + 273.15)))) / 
    (1+exp((delsC*(Tleaf+273.15)-EdJ)/(Rgas*(Tleaf+273.15))))
  
  J1 <- 1+exp((298.15*delsJ-EdJ)/Rgas/298.15)
  J2 <- 1+exp(((Tleaf+273.15)*delsJ-EdJ)/Rgas/(Tleaf+273.15))
  Jmax <- Jmax*exp(EaJ/Rgas*(1/298.15 - 1/(Tleaf+273.15)))*J1/J2
  J <- Jfun(PPFD, alpha, Jmax, theta)
  VJ <- J/4
  
  Ac <- Vcmax*(Ci - GammaStar)/(Ci + Km) - Rd
  Aj <- VJ * (Ci-GammaStar)/(Ci + 2*GammaStar) - Rd
  

  Am <- pmin(Ac,Aj)
  
  
#   hmshape <- 0.9999
#   Am <- (Ac+Aj - sqrt((Ac+Aj)^2-4*hmshape*Ac*Aj))/(2*hmshape) - Rd
  
  vec <- c(Ci,Ac,Aj,Am)
  names(vec) <- c("Ci","Ac","Aj", "Am")
  return(vec)
}



acifun <- function(...){
  
  x <- mapply(acifun_uv, ...)
  return(as.data.frame(t(x)))
  
}





