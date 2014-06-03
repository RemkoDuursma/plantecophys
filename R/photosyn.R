#' @title Coupled leaf gas exchange model
#' @description A coupled photosynthesis - stomatal conductance model, based on the Farquhar model of photosynthesis, and a Ball-Berry type model of stomatatal conductance. Includes temperature sensitivity of photosynthetic parameters, dark respiration (optionally calculated from leaf temperature), and mesophyll conductance. 
#' @param VPD Vapour pressure deficit (kPa)
#' @param Ca Atmospheric CO2 concentration (ppm)
#' @param PPFD Photosynthetic photon flux density ('PAR') (mu mol m-2 s-1)
#' @param Tleaf Leaf temperature (degrees C)
#' @param Patm Atmospheric pressure (kPa)
#' @param RH Relative humidity (in \%)
#' @param gsmodel One of BBOpti (Medlyn et al. 2011), BBLeuning (Leuning 1995), or BallBerry (Ball et al. 1987)
#' @param g0,g1 Parameters of Ball-Berry type stomatal conductance models.
#' @param gk Optional, exponent of VPD in gs model (Duursma et al. 2013)
#' @param vpdmin Below vpdmin, VPD=vpdmin, to avoid very high gs.
#' @param D0 Parameter for the BBLeuning stomatal conductance model.
#' @param alpha Quantum yield of electron transport (mol mol-1)
#' @param theta Shape of light response curve.
#' @param Jmax Maximum rate of electron transport at 25 degrees C (mu mol m-2 s-1)
#' @param Vcmax Maximum carboxylation rate at 25 degrees C (mu mol m-2 s-1)
#' @param gmeso Mesophyll conductance (mol m-2 s-1). If not NULL (the default), Vcmax and Jmax are chloroplastic rates.
#' @param Rd Dark respiration rate (mu mol m-2 s-1), optional (if not provided, calculated from Tleaf)
#' @param Rd0 Dark respiration rata at reference temperature (\code{TrefR})
#' @param Q10 Temperature sensitivity of Rd.
#' @param TrefR Reference temperature for Rd (Celcius).
#' @param Rdayfrac Ratio of Rd in the light vs. in the dark.
#' @param EaV, EdVC, delsC Vcmax temperature response parameters
#' @param EaJ, EdVJ, delsJ Jmax temperature response parameters
#' @param Ci Optional, intercellular CO2 concentration (ppm). If not provided, calculated via gs model.
#' @param whichA Which assimilation rate does gs respond to? 
#' @seealso FARAO fitaci AciC4
#' @details The coupled photosynthesis - stomatal conductance model finds the intersection between the supply of CO2 by diffusion, and the demand for CO2 by photosynthesis. See Farquhar and Sharkey (1982) for basic description of this type of model. 
#' 
#' The model of Farquhar et al. (1980) is used to estimate the dependence of photosynthesis rate on Ci.
#'  
#'  The temperature response of photosynthetic parameters, including Vcmax, Jmax, Gammastar, and Km follow Medlyn et al. 2002. 
#'  
#'  At the moment, two stomatal conductance models are implemented, both are Ball-Berry type models. The 'BBOpti' model is a slightly more general form of the model of Medlyn et al. 2011 (see Duursma et al. 2013). It is given by (in notation of the parameters and output variables of \code{Photosyn}),
#'  
#'  GS = G0 + 1.6*(1 + G1/D^GK)*ALEAF/CA
#'  
#'  where GK = 0.5 if stomata behave optimally (cf. Medlyn et al. 2011).
#'  
#'  The 'BBLeuning' model is that of Leuning (1995). It is given by,
#'  
#'  GS = G0 + g1*ALEAF/(Ca * (1 + VPD/D0))
#'  
#'  Note that this model also uses the g1 parameter, but it needs to be set to a much higher value to be comparable in magnitude to the BBOpti model.
#'  
#'  The 'BallBerry' model is that of Ball et al. (1987). It is given by,
#'  
#'  GS = G0 + g1*RH*ALEAF/Ca
#'  
#'  Where RH is relative humidity.
#'  
#'  For the full numerical solution to the Cowan-Farquhar optimization, use the \code{\link{FARAO}} function.
#'  
#'  If the mesophyll conductance is provided, it is assumed that Vcmax and Jmax are the chloroplastic rates, and leaf photosynthesis is calculated following Ethier and Livingston (2004).
#'  
#'  If Ci is provided as an input, this function calculates an A-Ci curve. Otherwise, Ci is calculated from the intersection between the 'supply' and 'demand' relationships, using the stomatal conductance model of Medlyn et al. (2011). 
#'  
#'  By default, the \code{Photosyn} function returns the hyperbolic minimum of Vcmax and Jmax-limited photosynthetic rates. This is to avoid the discontinuity at the transition between the two rates. Both Ac and Aj are also returned should they be needed. Note that those rates are output as gross photosynthetic rates!
#' @references 
#' Duursma, R.A., Payton, P., Bange, M.P., Broughton, K.J., Smith, R.A., Medlyn, B.E., Tissue, D. T., 2013,  Near-optimal response of instantaneous transpiration efficiency to vapour pressure deficit, temperature and [CO2] in cotton (Gossypium hirsutum L.). Agricultural and Forest Meteorology 168 : 168 - 176.
#'
#'Ethier, G. and N. Livingston. 2004. On the need to incorporate sensitivity to CO2 transfer conductance into the Farquhar von Caemmerer Berry leaf photosynthesis model. Plant, Cell & Environment. 27:137-153.
#'
#' Farquhar, G.D., S. Caemmerer and J.A. Berry. 1980. A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta. 149:78-90.
#' 
#' Farquhar, G. D., & Sharkey, T. D. (1982). Stomatal conductance and photosynthesis. Annual review of plant physiology, 33(1), 317-345.
#' 
#' Leuning, R. 1995. A critical-appraisal of a combined stomatal-photosynthesis model for C-3 plants. Plant Cell and Environment. 18:339-355.
#'
#' Medlyn, B.E., E. Dreyer, D. Ellsworth, M. Forstreuter, P.C. Harley, M.U.F. Kirschbaum, X. Le Roux, P. Montpied, J. Strassemeyer, A. Walcroft, K. Wang and D. Loustau. 2002. Temperature response of parameters of a biochemically based model of photosynthesis. II. A review of experimental data. Plant Cell and Environment. 25:1167-1179.
#' 
#' Medlyn, B.E., R.A. Duursma, D. Eamus, D.S. Ellsworth, I.C. Prentice, C.V.M. Barton, K.Y. Crous, P. De Angelis, M. Freeman and L. Wingate. 2011. Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology. 17:2134-2144.
#' @aliases Photosyn Aci
#' @return Returns a dataframe.
#' @examples
#' # Run the coupled leaf gas exchange model, set only a couple of parameters
#' Photosyn(VPD=2, g1=4, Ca=500)
#' 
#' # It is easy to set multiple values for inputs (and these can be mixed with single inputs);
#' r <- Photosyn(VPD=seq(0.5, 4, length=25), Vcmax=50, Jmax=100)
#' with(r, plot(VPD, ALEAF, type='l'))
#' 
#' # Set the mesophyll conductance
#' run1 <- Photosyn(PPFD=seq(50,1000,length=25), gmeso=0.15, Vcmax=40, Jmax=85)
#' with(run1, plot(PPFD, GS, type='l'))
#' 
#' # Run A-Ci curve only.
#' arun1 <- Aci(Ci=seq(50, 1200, length=101), Vcmax=40, Jmax=85)
#' arun2 <- Aci(Ci=seq(50, 1200, length=101), Vcmax=30, Jmax=70)
#' with(arun1, plot(Ci, ALEAF, type='l'))
#' with(arun2, points(Ci, ALEAF, type='l', lty=5))
#' @export
#' @rdname Photosyn
Photosyn <- function(VPD=1.5, 
                     Ca=400, 
                     PPFD=1500,
                     Tleaf=25,
                     Patm=101,
                     RH=NULL,
                     
                     gsmodel=c("BBOpti","BBLeuning","BallBerry"),
                     g1=4,
                     g0=0, 
                     gk=0.5,
                     vpdmin=0.5,
                     D0=5,
                      
                     alpha=0.24, 
                     theta=0.85, 
                     Jmax=100, 
                     Vcmax=50, 
                     gmeso=NULL,
                     
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
                     GS = NULL,
                     #AcCi=NULL,
                     Tcorrect=TRUE,  
                     returnParsOnly=FALSE,
                     whichA=c("Ah","Amin","Ac","Aj")){

  
  whichA <- match.arg(whichA)
  gsmodel <- match.arg(gsmodel)
  inputCi <- !is.null(Ci)
  inputGS <- !is.null(GS)
  if(inputCi & inputGS)stop("Cannot provide both Ci and GS")
  
  #---- Constants; hard-wired parameters.
  Rgas <- .Rgas()
  GCtoGW <- 1.57     # conversion from conductance to CO2 to H2O
  
  
  #---- Functions
  # Non-rectangular hyperbola
  Jfun <- function(PPFD, alpha, Jmax, theta){
    (alpha*PPFD + Jmax - 
       sqrt((alpha*PPFD + Jmax)^2 - 4*alpha*theta*PPFD*Jmax))/(2*theta)
  }
  
  
  #---- Do all calculations that can be vectorized
  
  # g1 and g0 are ALWAYS IN UNITS OF H20
  # G0 must be converted to CO2 (but not G1, see below)
  g0 <- g0/GCtoGW
  
  # Leaf respiration
  if(is.null(Rd))
    Rd <- Rdayfrac*Rd0*Q10^((Tleaf-TrefR)/10)
  
  # CO2 compensation point in absence of photorespiration
  GammaStar <- TGammaStar(Tleaf)
  
  # Michaelis-Menten coefficient
  Km <- TKm(Tleaf)
  
  #-- Vcmax, Jmax T responses
  if(Tcorrect){
    Vcmax <- Vcmax * TVcmax(Tleaf,EaV, delsC, EdVC)
    Jmax <- Jmax * TJmax(Tleaf, EaJ, delsJ, EdVJ)
  }
  
  #--- Stop here if only the parameters are required
  if(returnParsOnly){
    return(list(Vcmax=Vcmax, Jmax=Jmax, Km=Km, GammaStar=GammaStar))
  }
  
  # Electron transport rate
  J <- Jfun(PPFD, alpha, Jmax, theta)
  VJ <- J/4
  
  # Medlyn et al. 2011 model gs/A. NOTE: 1.6 not here because we need GCO2!
  if(gsmodel == "BBOpti"){
    vpduse <- VPD
    vpduse[vpduse < vpdmin] <- vpdmin
    GSDIVA <- (1 + g1/(vpduse^(1-gk)))/Ca
  }
  # Leuning 1995 model, without gamma (CO2 compensation point)
  if(gsmodel == "BBLeuning"){
    GSDIVA <- g1 / Ca / (1 + VPD/D0)
    GSDIVA <- GSDIVA / GCtoGW   # convert to conductance to CO2
  }
  if(gsmodel == "BallBerry"){
    if(is.null(RH))RH <- VPDtoRH(VPD, Tleaf)
    RH <- RH / 100
    GSDIVA <- g1 * RH / Ca
    GSDIVA <- GSDIVA / GCtoGW   # convert to conductance to CO2
  }  
    
  # If CI not provided, calculate from intersection between supply and demand
  if(!inputCi){
  
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
      
      CIJ[CIJ < GammaStar] <- GammaStar[CIJ < GammaStar]
      
      CIC <- Ci
    }
    
  
    # Photosynthetic rates, without or with mesophyll limitation
    if(is.null(gmeso)){
      # Get photosynthetic rate  
      Ac <- Vcmax*(CIC - GammaStar)/(CIC + Km)
      Aj <- VJ * (CIJ-GammaStar)/(CIJ + 2*GammaStar)
    
    } else {
    # Ethier and Livingston (2004) (Equation 10).
      A <- -1/gmeso
      BC <- (Vcmax - Rd)/gmeso + CIC + Km
      CC <- Rd*(CIC+Km)-Vcmax*(CIC-GammaStar)
      Ac <- mapply(QUADP, A=A,B=BC,C=CC)
      
      BJ <- (VJ - Rd)/gmeso + CIC + 2.0*GammaStar
      CJ <- Rd*(CIC+2.0*GammaStar) - VJ*(CIC - GammaStar)
      Aj <- mapply(QUADP, A=A,B=BJ,C=CJ)
      
      Ac <- Ac + Rd
      Aj <- Aj + Rd
      
    }
    
  
    # When below light-compensation points, assume Ci=Ca.
    if(!inputCi){
      lesslcp <- vector("logical", length(Aj))
      lesslcp <- Aj-Rd < 0
      
      if(length(Ca) == 1)Ca <- rep(Ca, length(CIJ))
      if(length(GammaStar) == 1)GammaStar <- rep(GammaStar, length(CIJ))
      if(length(VJ) == 1)VJ <- rep(VJ, length(CIJ))
      
      CIJ[lesslcp] <- Ca[lesslcp]
      Aj[lesslcp] <- VJ[lesslcp] * (CIJ[lesslcp] - GammaStar[lesslcp]) / 
        (CIJ[lesslcp] + 2*GammaStar[lesslcp])

      Ci <- ifelse(Aj < Ac, CIJ, CIC)
    }

    # Hyperbolic minimum.
    hmshape <- 0.9999
    Am <- (Ac+Aj - sqrt((Ac+Aj)^2-4*hmshape*Ac*Aj))/(2*hmshape) - Rd
#     
#     # If Ci below some set minimum, A will always be Rubisco limited.
#     # Used only for fitting stubborn A-Ci curves!
#     if(!is.null(AcCi)){
#       Am[Ci < AcCi] <- Ac[Ci < AcCi]
#     }
    
  
    # Calculate conductance to CO2
    if(!inputCi){
      if(whichA == "Ah")GS <- g0 + GSDIVA*Am
      if(whichA == "Aj")GS <- g0 + GSDIVA*(Aj-Rd)
      if(whichA == "Ac")GS <- g0 + GSDIVA*(Ac-Rd)
    } else {
      if(whichA == "Ah")GS <- Am/(Ca - Ci)
      if(whichA == "Aj")GS <- (Aj-Rd)/(Ca - Ci)
      if(whichA == "Ac")GS <- (Ac-Rd)/(Ca - Ci)
    }

    # Output conductance to H2O
    GS <- GS*GCtoGW
    
    # Transpiration rate assuming perfect coupling.
    # Output units are mmol m-2 s-1
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
                      PPFD=PPFD,
                      Patm=Patm)
  
return(df)
}

#' #'@rdname Photosyn
#'@export
Aci <- function(Ci,...)Photosyn(Ci=Ci,...)





