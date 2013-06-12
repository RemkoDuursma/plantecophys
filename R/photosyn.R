#' @title Coupled leaf gas exchange model
#' 
#' @export
#' @rdname Photosyn
#' @description A coupled photosynthesis - stomatal conductance model, based on the Farquhar model of photosynthesis, and a Ball-Berry type model of stomatatal conductance. Includes temperature sensitivity of photosynthetic parameters, dark respiration.
#' @details The coupled photosynthesis - stomatal conductance model finds the intersection between the supply
#' of CO2 by diffusion, and the demand for CO2 by photosynthesis. See Farquhar and Sharkey (1982) for basic description of this type of model.
#'  \if{html}{\figure{supplydemand.jpg}{Supply-demand}}
#'  \if{latex}{\figure{supplydemand}{Supply-demand}} 
#'  
#'  Figure 1. \code{Photosyn} calculates the operating point from the intersection of the supply and demand curves, here shown as the red dot.
#'  
#'  The model of Farquhar et al. (1980) is used to estimate the dependence of photosynthesis rate on Ci.
#'  
#'  The temperature response of photosynthetic parameters, including Vcmax, Jmax, Gammastar, and Km follow Medlyn et al. 2002. 
#'  
#'  At the moment, only one stomatal conductance model is implemented, which is a slightly more general form of the model of Medlyn et al. 2011 (see Duursma et al. 2013). It is given by (in notation of the parameters and output variables of \code{Photosyn}),
#'  
#'  GS = G0 + 1.6*(1 + G1/D^GK)*ALEAF/CA
#'  
#'  where GK = 0.5 if stomata behave optimally (cf. Medlyn et al. 2011).
#' @param VPD Vapour pressure deficit (kPa)
#' @param Ca Atmospheric CO2 concentration (ppm)
#' @param PPFD Photosynthetic photon flux density ('PAR') (mu mol m-2 s-1)
#' @param Tleaf Leaf temperature (degrees C)
#' @param Patm Atmospheric pressure (kPa)
#' @param g0,g1 Parameters of Medlyn et al. (2011) model of stomatal conductance (gs).
#' @param gk Optional, exponent of VPD in gs model (Duursma et al. 2013)
#' @param vpdmin Below vpdmin, VPD=vpdmin, to avoid very high gs.
#' @param alpha Quantum yield of electron transport (mol mol-1)
#' @param theta Shape of light response curve.
#' @param Jmax Maximum rate of electron transport at 25 degrees C (mu mol m-2 s-1)
#' @param Vcmax Maximum carboxylation rate at 25 degrees C (mu mol m-2 s-1)
#' @param Rd Dark respiration rate, optional (if not provided, calculated from Tleaf)
#' @param Rd0 Dark respiration rata at reference temperature (\code{TrefR})
#' @param Q10 Temperature sensitivity of Rd.
#' @param TrefR Reference temperature for Rd.
#' @param Rdayfrac Ratio of Rd in the light vs. in the dark.
#' @param EaV, EdVC, delsC Vcmax temperature response parameters
#' @param EaJ, EdVJ, delsJ Jmax temperature response parameters
#' @param Ci Optional, intercellular CO2 concentration
#' @param whichA Which assimilation rate does gs respond to?
#' @references 
#'
#' Farquhar, G.D., S. Caemmerer and J.A. Berry. 1980. A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta. 149:78-90.
#' 
#' Farquhar, G. D., & Sharkey, T. D. (1982). Stomatal conductance and photosynthesis. Annual review of plant physiology, 33(1), 317-345.
#' 
#' Medlyn, B.E., E. Dreyer, D. Ellsworth, M. Forstreuter, P.C. Harley, M.U.F. Kirschbaum, X. Le Roux, P. Montpied, J. Strassemeyer, A. Walcroft, K. Wang and D. Loustau. 2002. Temperature response of parameters of a biochemically based model of photosynthesis. II. A review of experimental data. Plant Cell and Environment. 25:1167-1179.
#' 
#' Medlyn, B.E., R.A. Duursma, D. Eamus, D.S. Ellsworth, I.C. Prentice, C.V.M. Barton, K.Y. Crous, P. De Angelis, M. Freeman and L. Wingate. 2011. Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology. 17:2134-2144.
#' 
#' Duursma, R.A., Payton, P., Bange, M.P., Broughton, K.J., Smith, R.A., Medlyn, B.E., Tissue, D. T., 2013,  Near-optimal response of instantaneous transpiration efficiency to vapour pressure deficit, temperature and [CO2] in cotton (Gossypium hirsutum L.). Agricultural and Forest Meteorology 168 : 168 - 176.
#' 
#' @aliases Photosyn Aci
#' @details If Ci is provided as an input, this function calculates an A-Ci curve. Otherwise, Ci is calculated from the intersection between the 'supply' and 'demand' relationships, using the stomatal conductance model of Medlyn et al. (2011). 
#' @return A dataframe.
Photosyn <- function(VPD=1.5, 
                     Ca=400, 
                     PPFD=1500,
                     Tleaf=25,
                     Patm=101,
                     
                     g1=4,
                     g0=0, 
                     gk=0.5,
                     vpdmin=0.5,
                      
                     alpha=0.24, 
                     theta=0.85, 
                     Jmax=100, 
                     Vcmax=50, 
                       
                     Rd0 = 0.92,
                     Q10 = 1.92,
                     Rd=NULL,
                     TrefR = 25,
                     Rdayfrac = 1.0,
                     
                     EaV = 82620.87,
                     EdVC = 0,
                     delsC = 645.1013,
                     
                     EaJ = 39676.89,
                     EdVJ = 200000,
                     delsJ = 641.3615,
                     
                     Ci = NULL,         
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
    if(EdVC > 0){
      V1 <- (1+exp((delsC*(25 + 273.15)-EdVC)/(Rgas*(25 + 273.15))))
      V2 <- (1+exp((delsC*(Tleaf+273.15)-EdVC)/(Rgas*(Tleaf+273.15))))
      f <- V1/V2
    } else f <- 1
    
    exp((Tleaf-25)*EaV/(Rgas*(Tleaf+273.15)*(25 + 273.15))) * f
  }
  
  # Hard-wired parameters.
  Jmaxfun <- function(Jmax){
    J1 <- 1+exp((298.15*delsJ-EdVJ)/Rgas/298.15)
    J2 <- 1+exp(((Tleaf+273.15)*delsJ-EdVJ)/Rgas/(Tleaf+273.15))
    exp(EaJ/Rgas*(1/298.15 - 1/(Tleaf+273.15)))*J1/J2
  }
  
  
  # Do all calculations that can be vectorized, send as argument
  # to non-vectorized bits ('photosynF')
  
  # Constants; hard-wired parameters.
  Rgas <- 8.314
  
  
  # Pre-calculations
  
  # g1 and g0 are ALWAYS IN UNITS OF H20
  # G0 must be converted (but no G1, see below)
  g0 <- g0/1.57
  
  # Leaf respiration
  #! acclimation
  if(is.null(Rd))
    Rd <- Rdayfrac*Rd0*Q10^((Tleaf-TrefR)/10)
  
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
  
  # Medlyn et al. 2011 model gs/A. NOTE: 1.6 not here because we need GCO2!
  vpduse <- VPD
  vpduse[vpduse < vpdmin] <- vpdmin
  GSDIVA <- (1 + g1/(vpduse^(1-gk)))/Ca
  
  # If CI not provided
  if(is.null(Ci)){
  
    #--- non-vectorized workhorse
    #! This can be vectorized, if we exclude Zero PPFD first,
    #! move whichA to somewhere else (at the end), and don't take minimum
      #! of Ac, Aj (calc both, return pmin at the end).
      getCI <- function(VJ,GSDIVA,PPFD,VPD,Ca,Tleaf,vpdmin,g0,Rd,
                            Vcmax,Jmax,Km,GammaStar){
        
    
        if(PPFD == 0){
          vec <- c(Ca,Ca)
          return(vec)
        }
            
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
        
        if(A == 0)
          CIJ <- -C/B
        else
          CIJ <- (- B + sqrt(B*B - 4*A*C)) / (2*A)
        
        return(c(CIJ,CIC))
      }
    
    # get Ci
    x <- mapply(getCI, 
                VJ=VJ,
                GSDIVA = GSDIVA,
                PPFD=PPFD,
                VPD=VPD,
                Ca=Ca,
                Tleaf=Tleaf,
                vpdmin=vpdmin,
                g0=g0,
                Rd=Rd,
                Vcmax=Vcmax,
                Jmax=Jmax,
                Km=Km,
                GammaStar=GammaStar)
      
      CIJ <- x[1,]
      CIC <- x[2,]
    } else {
      
      # Ci provided (A-Ci function mode)
      CIJ <- Ci
      CIC <- Ci
    }
    
    # Get photosynthetic rate  
    Ac <- Vcmax*(CIC - GammaStar)/(CIC + Km)
    Aj <- VJ * (CIJ-GammaStar)/(CIJ + 2*GammaStar)
    
    # When below light-compensation points, Ci=Ca.
    lesslcp <- Aj-Rd < 0
    
    if(length(Ca) == 1)Ca <- rep(Ca, length(CIJ))
    if(length(GammaStar) == 1)GammaStar <- rep(GammaStar, length(CIJ))
    if(length(VJ) == 1)VJ <- rep(VJ, length(CIJ))
    
    CIJ[lesslcp] <- Ca[lesslcp]
    Aj[lesslcp] <- VJ[lesslcp] * (CIJ[lesslcp] - GammaStar[lesslcp]) / (CIJ[lesslcp] + 2*GammaStar[lesslcp])
    
      
    # Ci
    Ci <- vector("numeric",length(CIC))
    Ci <- ifelse(Aj < Ac, CIJ, CIC)
    
    # Hyperbolic minimum.
    hmshape <- 0.9999
    Am <- (Ac+Aj - sqrt((Ac+Aj)^2-4*hmshape*Ac*Aj))/(2*hmshape) - Rd
    
    # Conductance to CO2
    if(whichA == "Ah")GS <- g0 + GSDIVA*Am
    if(whichA == "Aj")GS <- g0 + GSDIVA*(Aj-Rd)
    if(whichA == "Ac")GS <- g0 + GSDIVA*(Ac-Rd)
    
    # H2O
    GS <- GS*1.57
    
    # Transpiration rate; perfect coupling.
    E <- 1000*GS*VPD/Patm
    
    df <- data.frame( Ci=Ci,
                      ALEAF=Am,
                      GS=GS,
                      ELEAF=E,
                      Ac=Ac,
                      Aj=Aj,
                      Rd=Rd,
                      VPD=VPD,
                      Tleaf=Tleaf,
                      Ca=Ca,
                      PPFD=PPFD)
  
return(df)
}

#'@rdname Photosyn
#'@export
Aci <- function(Ci,...)Photosyn(Ci=Ci,...)


