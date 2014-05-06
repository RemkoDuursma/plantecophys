

# The net leaf energy balance, given that we know Tleaf, gs
LeafEnergyBalance <- function(Tleaf = 21.5, Tair = 20, 
                              gs = 0.15,
                              PPFD = 1500, VPD = 2, Patm = 101,
                              Wind = 2, Wleaf = 0.02, 
                              StomatalRatio = 1,   # 2 for amphistomatous
                              LeafAbs = 0.86){   # Leaf absorptance of total solar radiation
                              
  
  # Constants
  Boltz <- 5.67 * 10^-8     # w M-2 K-4
  Emissivity <- 0.95        # -
  LatEvap <- 2.54           # MJ kg-1
  CPAIR <- 1010.0           # 1010
  
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
  
  # Stuff for Penman-Monteith - don't actually need this
  # # Const s in Penman-Monteith equation  (Pa K-1)
  # SLOPE <- (esat(Tair + 0.1) - esat(Tair)) / 0.1
  # 
  # # Radiation conductance (mol m-2 s-1)
  # Gradiation <- 4.*Boltz*(Tair+273.15)^3 * Emissivity / (CPAIR * AIRMA)
  
  # See Leuning et al (1995) PC&E 18:1183-1200 Appendix E
  # Boundary layer conductance for heat - single sided, forced convection
  CMOLAR <- Patm*1000 / (8.314 * Tk(Tair))   # .Rgas() in package...
  Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR
  
  # Free convection
  GRASHOF <- 1.6E8 * abs(Tleaf-Tair) * (Wleaf^3) # Grashof number
  Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25) / Wleaf * CMOLAR
  
  # Total conductance to heat (both leaf sides)
  Gbh <- 2*(Gbhfree + Gbhforced)
  
  # Boundary layer conductance for water
  Gbw <- StomatalRatio * 1.075 * Gbh  # Leuning 1995
  gw <- gs*Gbw/(gs + Gbw)
  
  # Sensible heat flux (W m-2)
  # (positive flux is heat loss from leaf)
  H <- -CPAIR * AIRDENS * (Gbh/CMOLAR) * (Tair - Tleaf)
  
  # Longwave radiation
  # (positive flux is heat loss from leaf)
  Rlongup <- Emissivity*Boltz*((Tleaf+273.15)^4)
  
  # Rnet - full calculation without using Gradiation.
  Rsol <- 2*PPFD/UMOLPERJ   # W m-2
  Rnet <- LeafAbs*Rsol - Rlongup
  
  # Transpiration rate
  # Jones 1992 Eq. 5.17 (this avoids the need for the Penmon-Monteith equation)
  e <- esat(Tleaf) - VPD*1000   # Water vapour partial pressure (Pa)
  ET <- gw*(AIRDENS*0.622/(Patm*1000))*(esat(Tleaf) - e)   # mol m-2 s-1
  
  # Latent heat loss
  lambdaET <- LHV * ET
  
  # Net energy balance for the leaf
  EnergyBal <- Rnet - lambdaET - H

return(EnergyBal)
}
  
  
# Calculate Tleaf from energy balance, given that we know gs
FindTleaf <- function(gs, Tair, ...){
  
  Tleaf <- uniroot(LeafEnergyBalance, interval=c(Tair-15, Tair+15), 
                   gs=gs, Tair=Tair, ...)$root
  
return(Tleaf)
}









  
