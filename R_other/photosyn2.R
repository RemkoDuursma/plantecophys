


photosynF <- function(VPD=1.5, 
                      Ca=380, 
                      g1=4,
                      g0=0, 
                      gk=0.5,
                      vpdmin=0.5,
                      PPFD=1500, 
                      Tleaf=25, 
                      alpha=0.24, 
                      theta=0.85, 
                      Jmax=100, 
                      Vcmax=50, 
                      Rd0=0.92, 
                      EaV = 82620.87,
                      EaJ = 39676.89,
                      EdJ = 200000,
                      delsJ = 641.3615,
                      delsC = 645.1013,
                      whichA=c("Ah","Amin","Ac","Aj"),
                      ...){
  
  
  
  
  whichA <- match.arg(whichA)
  
  # g1 and g0 are ALWAYS IN UNITS OF H20
  # G0 must be converted (but no G1, see below)
  g0 <- g0/1.6
  
  Jfun <- function(PPFD, alpha, Jmax, theta){
    (alpha*PPFD + Jmax - sqrt((alpha*PPFD + Jmax)^2 - 4*alpha*theta*PPFD*Jmax))/(2*theta)
  }
  
  
  Km <- exp(38.05-79430/(8.314*(Tleaf+273)))*(1+210/exp(20.3-36380/(8.314*(Tleaf+273))))
  GammaStar <- 42.75*exp(37830*(Tleaf-25)/(8.314*298.15*(Tleaf+273.15)))
  Qten <- 1.95 # as in MAESPA (exp(0.067*10))
  Rd <- Rd0*Qten^((Tleaf-25)/10)
  
  
  # Note that if PPFD =0, GS=0 not g0, to be consistent with MAESPA.
  if(PPFD == 0){
    vec <- c(Ca,-Rd,0,0,-Rd,-Rd,Rd,VPD,Tleaf)
    names(vec) <- c("CI","ALEAF","GS","ELEAF","Ac","Aj","Rd","VPD","TLEAF")
    return(vec)
  }
  
  # From Belinda's spreadsheet.
  # 
  # Ea <- 47590 #51560
  # EaJ <- 37259 #43790
  # EdJ <- 200000
  # delsJ <- 640.02 #644.434
  
  Rgas <- 8.314
  
  # updated, aug 30 (Yan-Shih)
#   EaV <- 82620.87
#   EaJ <- 39676.89
#   EdJ <- 200000
#   delsJ <- 641.3615
#   delsC <- 645.1013
  
  
  V1 <- (1+exp((delsC*(25 + 273.15)-EdJ)/(Rgas*(25 + 273.15))))
  V2 <- (1+exp((delsC*(Tleaf+273.15)-EdJ)/(Rgas*(Tleaf+273.15))))
  Vcmax <- Vcmax * exp((Tleaf-25)*EaV/(Rgas*(Tleaf+273.15)*(25 + 273.15))) * V1/V2    
  
  J1 <- 1+exp((298.15*delsJ-EdJ)/Rgas/298.15)
  J2 <- 1+exp(((Tleaf+273.15)*delsJ-EdJ)/Rgas/(Tleaf+273.15))
  Jmax <- Jmax*exp(EaJ/Rgas*(1/298.15 - 1/(Tleaf+273.15)))*J1/J2
  
  # Electron transport rate
  J <- Jfun(PPFD, alpha, Jmax, theta)
  
  
  # full model. NOTE: 1.6 not here because we need GCO2!
  vpduse <- VPD
  vpduse[vpduse < vpdmin] <- vpdmin
  GSDIVA <- (1 + g1/(vpduse^(1-gk)))/Ca
  
  
  # Following calculations are used for both BB & BBL models.
  # Solution when Rubisco activity is limiting
  A = g0 + GSDIVA * (Vcmax - Rd)
  B = (1. - Ca*GSDIVA) * (Vcmax - Rd) + g0 * (Km - Ca)- GSDIVA * (Vcmax*GammaStar + Km*Rd)
  C = -(1. - Ca*GSDIVA) * (Vcmax*GammaStar + Km*Rd) - g0*Km*Ca
  
  CIC = (- B + sqrt(B*B - 4*A*C)) / (2*A)
  
  # Solution when electron transport rate is limiting
  VJ <- J/4
  A = g0 + GSDIVA * (VJ - Rd)
  B = (1 - Ca*GSDIVA) * (VJ - Rd) + g0 * (2.*GammaStar - Ca)- 
    GSDIVA * (VJ*GammaStar + 2.*GammaStar*Rd)
  C = -(1 - Ca*GSDIVA) * GammaStar * (VJ + 2*Rd) - 
    g0*2*GammaStar*Ca
  
  CIJ = (- B + sqrt(B*B - 4*A*C)) / (2*A)
  
  Ac <- Vcmax*(CIC - GammaStar)/(CIC + Km)
  Aj <- VJ * (CIJ-GammaStar)/(CIJ + 2*GammaStar)
  
  Ci <- if(Ac < Aj)CIC else CIJ
  
  hmshape <- 0.9999
  Am <- (Ac+Aj - sqrt((Ac+Aj)^2-4*hmshape*Ac*Aj))/(2*hmshape) - Rd
  
  # Conductance to CO2
  if(whichA == "Ah")GS <- g0 + GSDIVA*Am
  if(whichA == "Aj")GS <- g0 + GSDIVA*(Aj-Rd)
  if(whichA == "Ac")GS <- g0 + GSDIVA*(Ac-Rd)
  
  # H2O
  GS <- GS*1.6
  
  E <- 1000*GS*VPD/101
  
  vec <- c(Ci,Am,GS,E,Ac,Aj,Rd,VPD,Tleaf)
  names(vec) <- c("CI","ALEAF","GS","ELEAF","Ac","Aj","Rd","VPD","TLEAF")
  return(vec)
}

photosyn2 <- function(...){
  
  x <- mapply(photosynF, ...)
  return(as.data.frame(t(x)))
  
}


# 
# # If gm=NULL, it is actually infinite!
# ACi <- function(Ci, gm=NULL, PPFD=1500, Tleaf=25, alpha=0.24, theta=0.85, Jmax=100, Vcmax=50, Rd0=0.92,
#                 what=c("Am","Ac","Aj")){
#   
#   
#   what <- match.arg(what)
#   
#   Jfun <- function(PPFD, alpha, Jmax, theta){
#     (alpha*PPFD + Jmax - sqrt((alpha*PPFD + Jmax)^2 - 4*alpha*theta*PPFD*Jmax))/(2*theta)
#   }
#   
#   # Ac and Aj functions when gm is not infinite.
#   # Belinda, pers com (August 14 2012)
#   Ac_fun <- function(Ci, gm, Vcmax, Rd, Km, GammaStar){
#     
#     a <- -1/gm
#     b <- (Vcmax-Rd)/gm + Ci + Km
#     c <- Rd*(Ci + Km) - Vcmax*(Ci - GammaStar)
#     
#     (-b+sqrt(b^2-4*a*c))/(2*a)
#   }
#   Aj_fun <- function(Ci, gm, J, Rd, Km, GammaStar){
#     
#     a <- -1/gm
#     b <- (J/4 - Rd)/gm +Ci+2*GammaStar
#     c <- Rd*(Ci + 2*GammaStar) - (J/4)*(Ci-GammaStar)
#     
#     (-b+sqrt(b^2-4*a*c))/(2*a)
#   }
#   
#   
#   
#   Km <- exp(38.05-79430/(8.314*(Tleaf+273)))*(1+210/exp(20.3-36380/(8.314*(Tleaf+273))))
#   GammaStar <- 42.75*exp(37830*(Tleaf-25)/(8.314*298.15*(Tleaf+273.15)))
#   Qten <- 1.95 # as in MAESPA (exp(0.067*10))
#   Rd <- Rd0*Qten^((Tleaf-25)/10)
#   
#   # From Belinda's spreadsheet.
#   # 
#   # Ea <- 47590 #51560
#   # EaJ <- 37259 #43790
#   # EdJ <- 200000
#   # delsJ <- 640.02 #644.434
#   
#   Rgas <- 8.314
#   
#   # updated, aug 30 (Yan-Shih)
#   Ea <- 82620.87
#   EaJ <- 39676.89
#   EdJ <- 200000
#   delsJ <- 641.3615
#   delsC <- 645.1013
#   
#   
#   # Vcmax <- Vcmax*exp(Ea/Rgas*(1/298.15 - 1/(Tleaf+273.15)))	
#   Vcmax <- Vcmax * exp((Tleaf-25)*Ea/(Rgas*(Tleaf+273.15)*(25 + 273.15))) * 
#     (1+exp((delsC*(25 + 273.15)-EdJ)/(Rgas*(25 + 273.15)))) / 
#     (1+exp((delsC*(Tleaf+273.15)-EdJ)/(Rgas*(Tleaf+273.15))))
#   
#   J1 <- 1+exp((298.15*delsJ-EdJ)/Rgas/298.15)
#   J2 <- 1+exp(((Tleaf+273.15)*delsJ-EdJ)/Rgas/(Tleaf+273.15))
#   Jmax <- Jmax*exp(EaJ/Rgas*(1/298.15 - 1/(Tleaf+273.15)))*J1/J2
#   J <- Jfun(PPFD, alpha, Jmax, theta)
#   
#   # gm not infinite
#   if(!is.null(gm)){
#     Ac <- Ac_fun(Ci, gm, Vcmax, Rd, Km, GammaStar)
#     Aj <- Aj_fun(Ci, gm, J, Rd, Km, GammaStar)
#   } else {
#     Ac <- Vcmax*(Ci - GammaStar)/(Ci + Km)
#     Aj <- (J/4) * (Ci-GammaStar)/(Ci + 2*GammaStar)
#   }
#   
#   
#   hmshape <- 0.999
#   Am <- (Ac+Aj - sqrt((Ac+Aj)^2-4*hmshape*Ac*Aj))/(2*hmshape) - Rd
#   
#   if(what=="Am")return(Am-Rd)
#   if(what=="Ac")return(Ac-Rd)
#   if(what=="Aj")return(Aj)
# }	




# curve(ACi(x, gm=0.5), from=40, to=1500)
# curve(ACi(x, gm=0.1), add=TRUE, lty=5)
# curve(ACi(x, gm=0.05), add=TRUE, lty=5)
# curve(ACi(x, gm=0.025), add=TRUE, lty=5)


# 
# photosynF <- function(VPD=1.5, Ca=380, g1=4, g0=0, gm=0.2, gm_max=gm,  ...){
#   
#   
#   O <- function(Ci,...){
#     A <- ACi(Ci, gm=gm_max, ...)
#     gs <- g0 + 1.6*(1 + g1/sqrt(VPD))*A/Ca
#     Ci2 <- Ca - A/(gs/1.6)
#     return((Ci-Ci2)^2)
#   }
#   
#   # Ci at optimal solution. What stomata are going for.
#   Ci <- optimize(O, interval=c(0,Ca),...)$minimum
#   A <- ACi(Ci, gm=gm_max, ...)
#   # gs as if A is operating at max gm.
#   gs <- g0 + 1.6*(1 + g1/sqrt(VPD))*A/Ca
#   
#   # Actual photosynthesis - effects of low gm.
#   A <- ACi(Ci, gm=gm, ...)
#   
#   # Actual Ci, given gs and A.
#   Ci <- Ca - A/(gs/1.6)
#   
#   vec <- c(Ci,A,gs,VPD,Ca)
#   names(vec) <- c("Ci","A","gs","VPD","Ca")
#   return(vec)
# }
# 
# 
# photosyn2 <- function(...){
#   
#   x <- mapply(photosynF, ...)
#   return(as.data.frame(t(x)))
#   
# }
