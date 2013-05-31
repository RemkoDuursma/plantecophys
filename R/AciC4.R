#' C4 Photosynthesis
#' 
#' @description An implementation of the A-Ci curve for C4 plants, based on von Caemmerer et al. (1999)
#' @author Rhys Whitley
#' @param ci
#' @param PPFD
#' @param tleaf
#' @param VPMAX25
#' @param VCMAX25
#' @param Vpr
#' @param alpha
#' @param gbs
#' @param O2
#' @param JMAX25
#' @param x
#' @param THETA
#' @param Q10
#' @param RD0
#' @param RTEMP
#' @param TBELOW
#' @param DAYRESP
#' @param Q10F
#' @param FRM
AciC4 <- function( ci,
	PPFD=1500, 
	tleaf = 25,
	VPMAX25=120, 
	VCMAX25=60, 
	Vpr=80,         # PEP regeneration (mu mol m-2 s-1)
	alpha=0.0,		 # Fraction of PSII activity in the bundle sheath
	gbs=3e-3,		# Bundle shead conductance (mol m-2 s-1)
	O2=210,			# Mesophyll O2 concentration
	JMAX25=400,		
	x=0.4, 			# Partitioning factor for electron transport
	THETA=0.7,      # Shape parameter of the non-rectangular hyperbola
	Q10 = 2.3,
	RD0=1, 
  RTEMP=25, 
  TBELOW=0, 
  DAYRESP=1, 
  Q10F=2, 
	FRM=0.5		# Fraction of dark respiration that is mesophyll respiration (Rm)
	)
{

	Tk <- tleaf+273.15
	
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
    Kc <- 650*Q10^((tleaf-25)/10)
    Kp <- 80*Q10^((tleaf-25)/10)
    Ko <- 450*Q10^((tleaf-25)/10)
    K <- Kc*(1+O2/Ko)
            
    # T effects according to Massad et al. (2007)
	Vcmax <- VCMAX25*Arrhenius(Tk, 67294, 144568, 472)
	Vpmax <- VPMAX25*Arrhenius(Tk, 70373, 117910, 376)
	Jmax <- JMAX25*Arrhenius(Tk, 77900, 191929, 627)
	
	# # Dark and Mesophyll respiration
	# Rd <- 0.01*Vcmax
	# Rm <- 0.5*Rd
	
	# Day leaf respiration, umol m-2 s-1
    if (tleaf > TBELOW) {
        Rd <- RD0 * exp(Q10F * (tleaf - RTEMP)) * DAYRESP
    } else {
        Rd <- 0.0
    }
	Rm <-  FRM*Rd
	
	# PEP carboxylation rate
	Vp <- pmin(ci*Vpmax/(ci+Kp),Vpr)
	
	# Quadratic solution for enzyme limited C4 assimilation
	a.c <- 1 - (alpha*Kc)/(0.047*Ko)
	b.c <- -( (Vp-Rm+gbs*ci) + (Vcmax-Rd) + gbs*K + 
			alpha*low_gammastar/0.047*( low_gammastar*Vcmax+Rd*Kc/Ko ) )
	c.c <- (Vcmax-Rd)*(Vp-Rm+gbs*ci) - (Vcmax*gbs*low_gammastar*O2 + Rd*gbs*K)
		
	A.enzyme <- (-b.c - sqrt(b.c^2 - 4*a.c*c.c)) / (2*a.c)

	# Non-rectangular hyperbola describing light effect on electron transport rate (J)
    Qp2 <- PPFD*(1-0.15)/2
    J <- (1/(2*THETA))*(Qp2+Jmax - sqrt((Qp2+Jmax)^2-4*THETA*Qp2*Jmax))
	
	# Quadratic solution for light-limited C4 assimilation
	a.j <- 1 - 7*low_gammastar*alpha/(3*0.047)
	b.j <- -( (x*J/2-Rm+gbs*ci) + ((1-x)*J/3-Rd) + gbs*(7*low_gammastar*O2/3)
				+ alpha*low_gammastar/0.047*((1-x)*J/3+Rd) )
	c.j <- ( (x*J/2-Rm+gbs*ci)*((1-x)*J/3-Rd) - gbs*low_gammastar*O2*((1-x)*J/3-7*Rd/3) )
		
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
		
	return(list(An=An, Ac=Ac, Aj=Aj, Ad=Ad, Vp=Vp, Rd=Rd))
}
