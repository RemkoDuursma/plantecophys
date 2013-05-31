# This function returns the 'objective function' A - lambda*E
# This is to be optimized by the next function by varying ci.
OPTfun <- function(ci,              # mu mol mol-1
                    lambda=0.002,   # mol mol-1
                    ca=380,         # mu mol mol-1
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
	
    if(!C4){
		if(photo == "BOTH")A <- unname(unlist(Aci(CI=ci,...)["ALEAFHM"]))  # uses the hyperbolic minimum
		if(photo == "VCMAX")A <- unname(unlist(Aci(CI=ci,...)["AC"]))  # vcmax limited rate
		if(photo == "JMAX")A <- unname(unlist(Aci(CI=ci,...)["AJ"]))  # jmax limited rate
		
	}
    if(C4){
		if(photo == "BOTH")A <- AciC4(ci,...)$Ad  # uses the hyperbolic minimum
		if(photo == "VCMAX")A <- AciC4(ci,...)$Ac  # vcmax/vpmax limited
		if(photo == "JMAX")A <- AciC4(ci,...)$Aj  # jmax limited
	}
    gs <- a* A / (ca - ci)
	
	if(is.na(gbl))
		gtot <- gs
	else
		gtot <- 1/(1/gs + 1/gbl)
	
    E <- gtot*VPDmol  
                      
    objfun <- 10^-6*A - lambda*E

if(retobjfun)return(objfun)

if(!retobjfun)return(list(A=A, gs=gs, E=E, ci=ci))
}         