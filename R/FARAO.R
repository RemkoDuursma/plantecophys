
# FARquhar And Opti
FARAO <- function(Ca=400, ...){
  
  fx <- function(Ca,...)optimize(OPTfun, interval=c(0,Ca), maximum=TRUE,Ca=ca,...)$maximum
  optimalcis <- mapply(fx,Ca=Ca,...)
  
  res <- as.data.frame(OPTfun(Ci=optimalcis, retobjfun=FALSE, Ca=Ca, ...))
  
  return(res)
}

# This function returns the 'objective function' A - lambda*E
# This is to be optimized by the next function by varying ci.
OPTfun <- function(Ci,              # mu mol mol-1
                    lambda=0.002,   # mol mol-1
                    Ca=400,         # mu mol mol-1
                    VPD=1.5,        # kPa
					          Pa=101,         # ambient pressure, kPa
					          gbl=NA,
					          photo=c("BOTH","VCMAX","JMAX"),
                    retobjfun=TRUE, # if false, returns A, g and E (otherwise sum(A-lambda*E))
					          C4=FALSE,
					          ...){     

	a <- 1.6
	VPDmol <- VPD/Pa
	
	photo <- match.arg(photo)
	
  # Given a Ci, calculate photosynthetic rate
  if(!C4){
		if(photo == "BOTH")A <- Aci(Ci=ci,...)$ALEAF  # uses the hyperbolic minimum
		if(photo == "VCMAX")A <- Aci(Ci=ci,...)$Ac    # vcmax limited rate
		if(photo == "JMAX")A <- Aci(Ci=ci,...)$Aj     # jmax limited rate
	}
  if(C4){
		if(photo == "BOTH")A <- AciC4(ci,...)$Ad  # uses the hyperbolic minimum
		if(photo == "VCMAX")A <- AciC4(ci,...)$Ac  # vcmax/vpmax limited
		if(photo == "JMAX")A <- AciC4(ci,...)$Aj  # jmax limited
	}
  
  # Given Ci and A, calculate gs (diffusion constraint)
  gs <- a* A / (ca - ci)
	
  # Boundary layer conductance (optional)
	if(is.na(gbl))
		gtot <- gs
	else
		gtot <- 1/(1/gs + 1/gbl)
	
  # Transpiration rare
  E <- gtot*VPDmol  
  
  # Objective function to be maximizes (Cowan-Farquhar condition)
  objfun <- 10^-6*A - lambda*E

if(retobjfun)return(objfun)

if(!retobjfun)return(list(A=A, gs=gs, E=E, Ci=Ci))
}         