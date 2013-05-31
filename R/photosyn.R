
# !! need new name, to avoid confusion with GasExchangeR::photosyn
Photosyn <- function(VPD=1.5, 
                     Ca=380, 
                     g1=4,
                     g0=0, 
                     gk=0.5,
                     vpdmin=0.5,
                     PPFD=1500, 
                     alpha=0.24, 
                     theta=0.85, 
                     Jmax=100, 
                     Vcmax=50, 
                      
                     Tleaf=25, 
                     Rd0 = 0.92,
                     Q10 = 1.92,
                     EaV = 82620.87,
                     EaJ = 39676.89,
                     EdJ = 200000,
                     delsJ = 641.3615,
                     delsC = 645.1013,
                     
                     whichA=c("Ah","Amin","Ac","Aj")){

  
  whichA <- match.arg(whichA)
  
  
  # Functions
  
  # Non-rectangular hyperbola
  Jfun <- function(PPFD, alpha, Jmax, theta){
    (alpha*PPFD + Jmax - 
       sqrt((alpha*PPFD + Jmax)^2 - 4*alpha*theta*PPFD*Jmax))/(2*theta)
  }
  
  # Hard-wired parameters.
  Kmfun <- function(Tleaf){
    exp(38.05-79430/(8.314*(Tleaf+273)))*
      (1+210/exp(20.3-36380/(8.314*(Tleaf+273))))
  }
  
  # Hard-wired parameters.
  GammaStarfun <- function(Tleaf){
    42.75*exp(37830*(Tleaf-25)/(8.314*298.15*(Tleaf+273.15)))
  }
  
  # Hard-wired parameters.
  Vcmaxfun <- function(Tleaf){
    V1 <- (1+exp((delsC*(25 + 273.15)-EdJ)/(Rgas*(25 + 273.15))))
    V2 <- (1+exp((delsC*(Tleaf+273.15)-EdJ)/(Rgas*(Tleaf+273.15))))
    exp((Tleaf-25)*EaV/(Rgas*(Tleaf+273.15)*(25 + 273.15))) * V1/V2 
  }
  
  # Hard-wired parameters.
  Jmaxfun <- function(Jmax){
    J1 <- 1+exp((298.15*delsJ-EdJ)/Rgas/298.15)
    J2 <- 1+exp(((Tleaf+273.15)*delsJ-EdJ)/Rgas/(Tleaf+273.15))
    exp(EaJ/Rgas*(1/298.15 - 1/(Tleaf+273.15)))*J1/J2
  }
  
  
  # Do all calculations that can be vectorized, send as argument
  # to non-vectorized bits ('photosynF')
  
  # Constants; hard-wired parameters.
  Rgas <- 8.314
  
  
  # Pre-calculations
  
  # g1 and g0 are ALWAYS IN UNITS OF H20
  # G0 must be converted (but no G1, see below)
  g0 <- g0/1.6
  
  # Leaf respiration
  #! need dayfrac
  #! acclimation
  Rd <- Rd0*Q10^((Tleaf-25)/10)
  
  # Km, GammaStar
  Km <- Kmfun(Tleaf)
  GammaStar <- GammaStarfun(Tleaf)
  
  #-- Vcmax
  # Need function flexibility here.
  Vcmax <- Vcmax * Vcmaxfun(Tleaf)
  Jmax <- Jmax * Jmaxfun(Tleaf)
  
  # Electron transport rate
  J <- Jfun(PPFD, alpha, Jmax, theta)
  VJ <- J/4
  
  
  #--- non-vectorized workhorse
  #! This can be vectorized, if we exclude Zero PPFD first,
  #! move whichA to somewhere else (at the end), and don't take minimum
  #! of Ac, Aj (calc both, return pmin at the end).
  photosynF <- function(VJ,PPFD,Ca,Rd,Tleaf,vpdmin,g0,g1,gk,
                        Vcmax,Jmax,Km,GammaStar,whichA){
    

    # Note that if PPFD =0, GS=0 not g0, to be consistent with MAESPA.
    #! Keep it this way as a reminder that night-time not covered by this.
    #! This needs to be updated when we change output; annoying
    if(PPFD == 0){
      vec <- c(Ca,-Rd,0,0,-Rd,-Rd,Rd,VPD,Tleaf)
      return(vec)
    }
        
    # full model. NOTE: 1.6 not here because we need GCO2!
    vpduse <- VPD
    vpduse[vpduse < vpdmin] <- vpdmin
    GSDIVA <- (1 + g1/(vpduse^(1-gk)))/Ca
    
    # Taken from MAESTRA.
    # Following calculations are used for both BB & BBL models.
    # Solution when Rubisco activity is limiting
    A <- g0 + GSDIVA * (Vcmax - Rd)
    B <- (1. - Ca*GSDIVA) * (Vcmax - Rd) + g0 * 
      (Km - Ca)- GSDIVA * (Vcmax*GammaStar + Km*Rd)
    C <- -(1. - Ca*GSDIVA) * (Vcmax*GammaStar + Km*Rd) - g0*Km*Ca
    
    CIC <- (- B + sqrt(B*B - 4*A*C)) / (2*A)
    
    # Solution when electron transport rate is limiting
    A <- g0 + GSDIVA * (VJ - Rd)
    B <- (1 - Ca*GSDIVA) * (VJ - Rd) + g0 * (2.*GammaStar - Ca)- 
      GSDIVA * (VJ*GammaStar + 2.*GammaStar*Rd)
    C <- -(1 - Ca*GSDIVA) * GammaStar * (VJ + 2*Rd) - 
      g0*2*GammaStar*Ca
    
    CIJ <- (- B + sqrt(B*B - 4*A*C)) / (2*A)
    
    Ac <- Vcmax*(CIC - GammaStar)/(CIC + Km)
    Aj <- VJ * (CIJ-GammaStar)/(CIJ + 2*GammaStar)
    
    
    Ci <- if(Ac < Aj)CIC else CIJ
    
    # Hyperbolic minimum.
    hmshape <- 0.9999
    Am <- (Ac+Aj - sqrt((Ac+Aj)^2-4*hmshape*Ac*Aj))/(2*hmshape) - Rd
    
    # Conductance to CO2
    if(whichA == "Ah")GS <- g0 + GSDIVA*Am
    if(whichA == "Aj")GS <- g0 + GSDIVA*(Aj-Rd)
    if(whichA == "Ac")GS <- g0 + GSDIVA*(Ac-Rd)
    
    # H2O
    GS <- GS*1.6
    
    # Transpiration rate; perfect coupling.
    E <- 1000*GS*VPD/101
    
    vec <- c(Ci,Am,GS,E,Ac,Aj,Rd,VPD,Tleaf)
    return(vec)
  }

  
  # Do runs
  x <- mapply(photosynF, 
              VJ=VJ,
              PPFD=PPFD,
              Ca=Ca,
              Rd=Rd,
              Tleaf=Tleaf,
              vpdmin=vpdmin,
              g0=g0,
              g1=g1,
              gk=gk,
              Vcmax=Vcmax,
              Jmax=Jmax,
              Km=Km,
              GammaStar=GammaStar,
              whichA=whichA)
  
  df <- as.data.frame(t(x))
  names(df) <- c("CI","ALEAF","GS","ELEAF","Ac","Aj","Rd","VPD","TLEAF")
return(df)
}

