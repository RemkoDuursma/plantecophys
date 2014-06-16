#' Coupled leaf gas exchange model with energy balance
#' @description As \code{\link{Photosyn}}, but calculates the leaf temperature based on the leaf's energy balance. Including sensible and long-wave heat loss, latent heat loss from evaporation, and solar radiation input. 
#' 
#'@details Uses the Penman-Monteith equation to calculate the leaf transpiration rate, and finds Tleaf by solving the leaf energy balance iteratively. In the solution, it is accounted for that stomatal conductance and net radiation depend on Tleaf. There are no simplifications or approximations made to find the solution.
#'@param Wind Wind speed (m s-1)
#'@param VPD The vapour pressure deficit of the air (i.e. not the leaf-to-air VPD)
#'@param Wleaf Leaf width (m)
#'@param StomatalRatio The stomatal ratio (cf. Licor6400 terminology), if it is 1, leaves have stomata only on one side (hypostomatous), 2 for leaves with stomata on both sides (amphistomatous)
#'@param LeafAbs Leaf absorptance of solar radiation
#'@param RH The relative humidity of the air (i.e. not calculated with leaf temperature)
#'@export PhotosynEB
PhotosynEB <- function(Tair=25,
                       VPD=1.5,
                       Wind = 2,   
                       Wleaf = 0.02, 
                       StomatalRatio = 1,   # 2 for amphistomatous
                       LeafAbs = 0.86, 
                       penmon=c("iso","full"),
                       RH=NULL,
                       ...){
  
  
  if(!is.null(RH))VPD <- RHtoVPD(RH,Tair)
  penmon <- match.arg(penmon)
  
  # Non-vectorized function declared here.
  PhotosynEBfun <- function(Tair,
                            VPD,
                            Wind,   
                            Wleaf, 
                            StomatalRatio,   # 2 for amphistomatous
                            LeafAbs,
                            penmon,
                            ...){   # passed to Photosyn
    
    
    m <- match.call(expand.dots=TRUE)
    if("Tleaf" %in% names(m))
      stop("Cannot pass Tleaf to PhotosynEB - it is calculated from energy balance.")
    
    # Wrapper for Photosyn to find gs only  
    gsfun <- function(...)Photosyn(...)$GS
    
    # Find Tleaf. Here, we take into account that Tleaf as solved from
    # energy balance affects gs, so this is the second loop to solve for Tleaf.
    fx <- function(x, Tair, Wind, VPD, Wleaf, StomatalRatio, LeafAbs, penmon, ...){
      newx <- FindTleaf(Tair=Tair, gs=gsfun(Tleaf=x, VPD=VPDairToLeaf(VPD,Tair,x), ...), 
                        Wind=Wind, Wleaf=Wleaf, 
                        StomatalRatio=StomatalRatio, LeafAbs=LeafAbs, penmon=penmon)
      (newx - x)^2
    }

    Tleaf <- optimize(fx, interval=c(Tair-15, Tair+15), Tair=Tair, Wind=Wind, Wleaf=Wleaf, 
                     VPD=VPD, StomatalRatio=StomatalRatio, LeafAbs=LeafAbs, penmon=penmon, ...)$minimum

    # Now run Photosyn
    p <- Photosyn(Tleaf=Tleaf, VPD=VPDairToLeaf(VPD,Tair,Tleaf), ...)
    
    # And energy balance components
    e <- LeafEnergyBalance(Tleaf=Tleaf, Tair=Tair, gs=p$GS, 
                           PPFD=p$PPFD, VPD=p$VPD, Patm=p$Patm, 
                           Wind=Wind, Wleaf=Wleaf, 
                           StomatalRatio=StomatalRatio, LeafAbs=LeafAbs, penmon=penmon,
                           returnwhat="fluxes")
    res <- cbind(p,e)
    
    # Replace ELEAF with energy-balance one.
    res$ELEAF <- res$ELEAFeb
    res$ELEAFeb <- NULL
    res$VPDleaf <- VPDairToLeaf(res$VPD, Tair, res$Tleaf)
    res$Tair <- Tair
    res$Wind <- Wind
    
    return(res)
  }
  
  m <- t(mapply(PhotosynEBfun, Tair=Tair, 
                VPD=VPD, 
                Wind=Wind, Wleaf=Wleaf, StomatalRatio=StomatalRatio, 
                LeafAbs=LeafAbs, penmon=penmon, ..., SIMPLIFY=FALSE))
  
  
  return(as.data.frame(do.call(rbind, m)))
}

# The net leaf energy balance, given that we know Tleaf, gs
LeafEnergyBalance <- function(Tleaf = 21.5, Tair = 20, 
                              gs = 0.15,
                              PPFD = 1500, VPD = 2, Patm = 101,
                              Wind = 2, Wleaf = 0.02, 
                              StomatalRatio = 1,   # 2 for amphistomatous
                              LeafAbs = 0.86,
                              penmon=c("iso","full"),
                              returnwhat=c("balance","fluxes")){   # Leaf absorptance of total solar radiation
                              
  
  returnwhat <- match.arg(returnwhat)
  penmon <- match.arg(penmon)
  
  # Constants
  Boltz <- 5.67 * 10^-8     # w M-2 K-4
  Emissivity <- 0.95        # -
  LatEvap <- 2.54           # MJ kg-1
  CPAIR <- 1010.0           # J kg-1 K-1
  
  H2OLV0 <- 2.501e6         # J kg-1
  H2OMW <- 18e-3            # J kg-1
  AIRMA <- 29.e-3           # mol mass air (kg/mol)
  AIRDENS <- 1.204          # kg m-3
  UMOLPERJ <- 4.57
  DHEAT <- 21.5e-6          # molecular diffusivity for heat

  
  
  # Density of dry air
  AIRDENS <- Patm*1000/(287.058 * Tk(Tair))
  
  # Latent heat of water vapour at air temperature (J mol-1)
  LHV <- (H2OLV0 - 2.365E3 * Tair) * H2OMW
  
  # Const s in Penman-Monteith equation  (Pa K-1)
  SLOPE <- (esat(Tair + 0.1) - esat(Tair)) / 0.1
  
  # Radiation conductance (mol m-2 s-1)
  Gradiation <- 4.*Boltz*Tk(Tair)^3 * Emissivity / (CPAIR * AIRMA)

  # See Leuning et al (1995) PC&E 18:1183-1200 Appendix E
  # Boundary layer conductance for heat - single sided, forced convection
  CMOLAR <- Patm*1000 / (8.314 * Tk(Tair))   # .Rgas() in package...
  Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR
  
  # Free convection
  GRASHOF <- 1.6E8 * abs(Tleaf-Tair) * (Wleaf^3) # Grashof number
  Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25) / Wleaf * CMOLAR
  
  # Total conductance to heat (both leaf sides)
  Gbh <- 2*(Gbhfree + Gbhforced)
  
  # Heat and radiative conductance
  Gbhr <- Gbh + 2*Gradiation
  
  # Boundary layer conductance for water (mol m-2 s-1)
  Gbw <- StomatalRatio * 1.075 * Gbh  # Leuning 1995
  gw <- gs*Gbw/(gs + Gbw)
  
  # Sensible heat flux (W m-2)
  # (positive flux is heat loss from leaf)
  H <- -CPAIR * AIRDENS * (Gbh/CMOLAR) * (Tair - Tleaf)
  
  # Longwave radiation
  # (positive flux is heat loss from leaf)
  Rlongup <- Emissivity*Boltz*Tk(Tleaf)^4
  
  # Rnet
  Rsol <- 2*PPFD/UMOLPERJ   # W m-2
  Rnet <- LeafAbs*Rsol - Rlongup   # full
  Rnetiso <- LeafAbs*Rsol - Emissivity*Boltz*Tk(Tair)^4 # isothermal net radiation
  
  # Transpiration rate
  # Penman-Monteith
  GAMMA <- CPAIR*AIRMA*Patm*1000/LHV

  # Standard Penmon-Monteith (5.23 in Jones 1992)
  if(penmon=="full")ET <- (1/LHV) * (SLOPE * Rnet + 1000*VPD * Gbh * CPAIR * AIRMA) / (SLOPE + GAMMA * Gbh/gw)
  
  # iso-thermal version
  if(penmon=="iso")ET <- (1/LHV) * (SLOPE * Rnetiso + 1000*VPD * Gbh * CPAIR * AIRMA) / (SLOPE + GAMMA * Gbhr/gw)

  # Latent heat loss
  lambdaET <- LHV * ET
  
  # Net energy balance for the leaf
  EnergyBal <- Rnet - lambdaET - H

  if(returnwhat == "balance")return(EnergyBal)
  
  if(returnwhat == "fluxes"){
    
    l <- data.frame(ELEAFeb=1000*ET, Rnet=Rnet, Rlongup=Rlongup, 
                    H=H, lambdaET=lambdaET, gw=gw, Gbh=Gbh)
    return(l)
  }
}
  
  
# Calculate Tleaf from energy balance, given that we know gs
FindTleaf <- function(gs, Tair, penmon,  ...){
  
  Tleaf <- try(uniroot(LeafEnergyBalance, interval=c(Tair-15, Tair+15), 
                   gs=gs, Tair=Tair, penmon=penmon, ...)$root)
   if(inherits(Tleaf, "try-error"))browser()
  
return(Tleaf)
}








