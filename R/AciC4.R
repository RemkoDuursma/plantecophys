#' C4 Photosynthesis
#' 
#' @description An implementation of the A-Ci curve for C4 plants, based on von Caemmerer et al. (1999)
#' @author Rhys Whitley
#' @param Ci
#' @param PPFD
#' @param Tleaf
#' @param VPMAX25
#' @param Vcmax
#' @param Vpr PEP regeneration (mu mol m-2 s-1)
#' @param alpha Fraction of PSII activity in the bundle sheath
#' @param gbs Bundle shead conductance (mol m-2 s-1)
#' @param O2 Mesophyll O2 concentration
#' @param JMAX25
#' @param x Partitioning factor for electron transport
#' @param THETA Shape parameter of the non-rectangular hyperbola
#' @param Q10
#' @param RD0
#' @param RTEMP
#' @param TBELOW
#' @param DAYRESP
#' @param Q10F
#' @param FRM Fraction of dark respiration that is mesophyll respiration (Rm)
#' @export
AciC4 <- function(Ci,
	PPFD=1500, 
	Tleaf = 25,
	VPMAX25=120, 
	Vcmax=60, 
	Vpr=80,         
	alpha=0.0,		  
	gbs=3e-3,		 
	O2=210,			 
	JMAX25=400,		
	x=0.4, 			 
	THETA=0.7,   
	Q10 = 2.3,
	RD0=1, 
  RTEMP=25, 
  TBELOW=0, 
  DAYRESP=1, 
  Q10F=2, 
	FRM=0.5	,...){

	Tk <- Tleaf+273.15
	
	# Temperature effects on Vcmax, Vpmax and Jmax (Massad et al. 2007)
    # This function returns value between 0 and 1.
    Arrhenius <- function(Tk, Ea, Hd, DS){
			R <- 8.3144
            fTk <- exp( Ea*(Tk-298)/(298*R*Tk) ) * 
            		(1+exp( (298*DS-Hd)/(298*R) )) / 
            		(1+exp( (Tk*DS-Hd)/(Tk*R) )) 
    }
	
	# Half the reciprocal for Rubisco specificity (NOT CO2 compensation point)
	low_gammastar <- 1.93e-4
	
	# Michaelis-Menten coefficients for CO2 (Kc, mu mol mol-1) and O (Ko, mmol mol-1) and combined (K)
  Kc <- 650*Q10^((Tleaf-25)/10)
  Kp <- 80*Q10^((Tleaf-25)/10)
  Ko <- 450*Q10^((Tleaf-25)/10)
  K <- Kc*(1+O2/Ko)
          
  # T effects according to Massad et al. (2007)
	Vcmax <- Vcmax*Arrhenius(Tk, 67294, 144568, 472)
	Vpmax <- VPMAX25*Arrhenius(Tk, 70373, 117910, 376)
	Jmax <- JMAX25*Arrhenius(Tk, 77900, 191929, 627)
	
	# # Dark and Mesophyll respiration
	# Rd <- 0.01*Vcmax
	# Rm <- 0.5*Rd
	
	# Day leaf respiration, umol m-2 s-1
  if (Tleaf > TBELOW) {
      Rd <- RD0 * exp(Q10F * (Tleaf - RTEMP)) * DAYRESP
  } else {
      Rd <- 0.0
  }
	Rm <-  FRM*Rd
	
	# PEP carboxylation rate
	Vp <- pmin(Ci*Vpmax/(Ci+Kp),Vpr)
	
	# Quadratic solution for enzyme limited C4 assimilation
	a.c <- 1 - (alpha*Kc)/(0.047*Ko)
	b.c <- -( (Vp-Rm+gbs*Ci) + (Vcmax-Rd) + gbs*K + 
			alpha*low_gammastar/0.047*( low_gammastar*Vcmax+Rd*Kc/Ko ) )
	c.c <- (Vcmax-Rd)*(Vp-Rm+gbs*Ci) - (Vcmax*gbs*low_gammastar*O2 + Rd*gbs*K)
		
	A.enzyme <- (-b.c - sqrt(b.c^2 - 4*a.c*c.c)) / (2*a.c)

	# Non-rectangular hyperbola describing light effect on electron transport rate (J)
    Qp2 <- PPFD*(1-0.15)/2
    J <- (1/(2*THETA))*(Qp2+Jmax - sqrt((Qp2+Jmax)^2-4*THETA*Qp2*Jmax))
	
	# Quadratic solution for light-limited C4 assimilation
	a.j <- 1 - 7*low_gammastar*alpha/(3*0.047)
	b.j <- -( (x*J/2-Rm+gbs*Ci) + ((1-x)*J/3-Rd) + gbs*(7*low_gammastar*O2/3)
				+ alpha*low_gammastar/0.047*((1-x)*J/3+Rd) )
	c.j <- ( (x*J/2-Rm+gbs*Ci)*((1-x)*J/3-Rd) - gbs*low_gammastar*O2*((1-x)*J/3-7*Rd/3) )
		
	A.light <- (-b.j - sqrt(b.j^2 - 4*a.j*c.j)) / (2*a.j)

	# Actual assimilation rate
	An <- pmin(A.enzyme,A.light)
	Ac <- A.enzyme
	Aj <- A.light
	
	# Hyperbolic minimum (Buckley), to avoid discontinuity at transition from Ac to Aj
	shape2 <- 0.999
	Ad <- (Ac+Aj - sqrt((Ac+Aj)^2-4*shape2*Ac*Aj))/(2*shape2) - Rd
	Ac <- Ac - Rd
	Aj <- Aj - Rd
		
	return(list(ALEAF=Ad, An=An, Ac=Ac, Aj=Aj, Vp=Vp, Rd=Rd, Tleaf=Tleaf, PPFD=PPFD))
}
