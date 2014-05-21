#' FARquhar And Opti
#' 
#' @description The numerical solution of the optimal stomatal conductance model, coupled with the Farquhar model of photosynthesis. The model of Medlyn et al. (2011) is an approximation to this full numeric solution.
#' @param lambda The marginal cost of water (mol mol-1)
#' @param Ca The CO2 concentration. 
#' @param photo Which photosynthesis rate should stomata respond to? Defaults to 'BOTH', i.e. the minimum of Vcmax and Jmax limited rates.
#' @param C4 If TRUE, uses the C4 photosynthesis routine (\code{\link{AciC4}})
#' @param ... All other parameters are passed to \code{\link{Aci}}
#' @author Remko Duursma
#' @details This model finds the Ci that maximizes A - lambda*E (Cowan & Farquhar 1977, see also Medlyn et al. 2011).
#' @references 
#' Cowan, I. and G.D. Farquhar. 1977. Stomatal function in relation to leaf metabolism and environment. Symposia of the Society for Experimental Biology. 31:471-505.
#' 
#' Medlyn, B.E., R.A. Duursma, D. Eamus, D.S. Ellsworth, I.C. Prentice, C.V.M. Barton, K.Y. Crous, P. De Angelis, M. Freeman and L. Wingate. 2011. Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology. 17:2134-2144.
#' @export
FARAO <- function(lambda=0.002, Ca=400, VPD=1, 
                  photo=c("BOTH","VCMAX","JMAX"), C4=FALSE, ...){
  
  photo <- match.arg(photo)
  
  fx <- function(Ca,...)optimize(OPTfun, interval=c(0,Ca), 
                                 maximum=TRUE,Ca=Ca,
                                 ...)$maximum
  optimalcis <- mapply(fx,Ca=Ca,lambda=lambda, photo=photo, C4=C4, VPD=VPD,...)
  
  res <- as.data.frame(OPTfun(Ci=optimalcis, retobjfun=FALSE, 
                              Ca=Ca, photo=photo, C4=C4, VPD=VPD,...))
  
  return(res)
}

# This function returns the 'objective function' A - lambda*E
# This is to be optimized by the next function by varying ci.
OPTfun <- function(Ci,              # mu mol mol-1
                    lambda=0.002,   # mol mol-1
                    Ca=400,         # mu mol mol-1
                    VPD=1.5,        # kPa
					          Patm=101,         # ambient pressure, kPa
					          photo=c("BOTH","VCMAX","JMAX"),
                    energybalance=FALSE,
                    retobjfun=TRUE, # if false, returns A, g and E (otherwise sum(A-lambda*E))
					          C4=FALSE,
					          ...){     

  GCtoGW <- 1.57
	VPDmol <- VPD/Patm
	
	photo <- match.arg(photo)
	
  # Given a Ci, calculate photosynthetic rate
  if(!C4)
		run <- Aci(Ci=Ci, VPD=VPD, ...)   # note that VPD does not do anything, just for consistency in I/O
	else 
		run <- AciC4(Ci, VPD=VPD, ...)
	
  if(photo == "BOTH")A <- run$ALEAF
	if(photo == "VCMAX")A <- run$Ac
	if(photo == "JMAX")A <- run$Aj
  
  # Given Ci and A, calculate gs (diffusion constraint)
  gs <- GCtoGW * A / (Ca - Ci)
	    
  # Transpiration rate
  E <- gs*VPDmol
    
    
  # Objective function to be maximized (Cowan-Farquhar condition)
  objfun <- 10^-6*A - lambda*E

if(retobjfun)return(objfun)

if(!retobjfun)return(list( Ci=Ci, ALEAF=A, GS=gs, ELEAF=E*1000, Ac=run$Ac, Aj=run$Aj,
                           Rd=run$Rd, VPD=VPD, Tleaf=run$Tleaf,  Ca=Ca, PPFD=run$PPFD ))
}         


OPTfunEB <- function(Ci,           # mu mol mol-1
                   lambda=0.002,   # mol mol-1
                   Ca=400,         # mu mol mol-1
                   VPD=1.5,        # AIR VPD! kPa
                   Patm=101,       # ambient pressure, kPa
                   Tair=25,
                   Wind=2,
                   Wleaf=0.02,
                   StomatalRatio=1,
                   LeafAbs=0.86,
                   photo=c("BOTH","VCMAX","JMAX"),
                   energybalance=FALSE,
                   retobjfun=TRUE, # if false, returns A, g and E (otherwise sum(A-lambda*E))
                   C4=FALSE,
                   ...){     
  
  GCtoGW <- 1.57
  VPDmol <- VPD/Patm
  
  photo <- match.arg(photo)
  
  gsfun <- function(Ci, VPD, returnwhat=c("gs","all"), ...){
    
    returnwhat <- match.arg(returnwhat)
    
    # Given a Ci, calculate photosynthetic rate
    if(!C4)
      run <- Aci(Ci, VPD=VPD, ...)   # note that VPD does not do anything, just for consistency in I/O
    else 
      run <- AciC4(Ci, VPD=VPD, ...)
    
    if(photo == "BOTH")A <- run$ALEAF
    if(photo == "VCMAX")A <- run$Ac
    if(photo == "JMAX")A <- run$Aj
    
    # Given Ci and A, calculate gs (diffusion constraint)
    gs <- GCtoGW * A / (Ca - Ci)
    
  if(returnwhat == "gs")return(gs)
  if(returnwhat == "all")return(list(run=run, GS=gs, A=A))
  }
  
  # Find Tleaf. Here, we take into account that Tleaf as solved from
  # energy balance affects gs (because it affects)
  fx <- function(x, Ci, Tair, Wind, VPD, Wleaf, StomatalRatio, LeafAbs, ...){
    newx <- FindTleaf(Tair=Tair, gs=gsfun(Ci=Ci, Tleaf=x, VPD=VPD, ...), 
                      Wind=Wind, Wleaf=Wleaf, 
                      StomatalRatio=StomatalRatio, LeafAbs=LeafAbs)
    newx - x
  }
  Tleaf <- uniroot(fx, interval=c(Tair-15, Tair+15), Ci=Ci, Tair=Tair, Wind=Wind, VPD=VPD, Wleaf=Wleaf, 
                   StomatalRatio=StomatalRatio, LeafAbs=LeafAbs, ...)$root
  z <- gsfun(Ci=Ci, Tleaf=Tleaf, VPD=VPD, returnwhat="all",...)
  GS <- z$GS
  A <- z$A
  
  # And energy balance components
  e <- LeafEnergyBalance(Tleaf=Tleaf, Tair=Tair, gs=GS, 
                         PPFD=z$run$PPFD, VPD=VPD, Patm=z$run$Patm, 
                         Wind=Wind, Wleaf=Wleaf, 
                         StomatalRatio=StomatalRatio, LeafAbs=LeafAbs,
                         returnwhat="fluxes")

  E <- e$ELEAFeb
  
  # Objective function to be maximized (Cowan-Farquhar condition)
  objfun <- 10^-6*A - lambda*E/1000
  
  if(retobjfun)return(objfun)
  
  if(!retobjfun)return(list( Ci=Ci, ALEAF=A, GS=gs, ELEAF=E*1000, Ac=z$run$Ac, Aj=z$run$Aj,
                             Rd=z$run$Rd, VPD=VPD, Tleaf=Tleaf,  Ca=Ca, PPFD=z$run$PPFD ))
}  




