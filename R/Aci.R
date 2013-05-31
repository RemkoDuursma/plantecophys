

# OBSOLETE. Only here for variable names.

acifun <- function(
CI,
PAR = 1000,
TLEAF = 20,
JMAX25 = 145.4,
VCMAX25 = 89.5,
RD0 = 0.9,
Q10F = 0.67,
RTEMP = 25,
IECO = 1,
EAVJ = 37259,
EDVJ = 200000,
DELSJ = 640.02,
THETA = 0.4,
AJQ = 0.324,
GMESO = -1.0,
EAVC = 47590,
EDVC = 0.0,
DELSC = 0.0,
TVJUP = -100,
TVJDN = -100,
DAYRESP = 1.0,
TBELOW = -100.0,
HMSHAPE = 0.999,
PATM=101000,
returncols = c('CI','PAR','TLEAF','ALEAF','AJ','AC','RD')
){


# Outputs need to be initialized
ALEAF <- -999
AJ <- -999
AC <- -999
RD <- -999
ALEAFHM <- -999
GAMMASTAR <- -999
VJ <- -999
KM <- -999
VCMAX <-  -999

if(any(is.na(c(TLEAF,PAR,CI))))return(rep(NA,length(returncols)))

f <- .Fortran("aci", as.double(CI),
                     as.double(PAR),
                     as.double(TLEAF),
					 as.double(PATM),
					 as.double(JMAX25),
					 as.integer(IECO),
					 as.double(EAVJ),
					 as.double(EDVJ),
					 as.double(DELSJ),
					 as.double(VCMAX25),
					 as.double(EAVC),
					 as.double(EDVC),
					 as.double(DELSC),
					 as.double(TVJUP),
					 as.double(TVJDN),
					 as.double(THETA),
					 as.double(AJQ),
					 as.double(GMESO),
					 as.double(RD0),
					 as.double(Q10F),
					 as.double(RTEMP),
					 as.double(DAYRESP),
					 as.double(TBELOW),
					 as.double(HMSHAPE),
					 as.double(AC),
					 as.double(AJ),
					 as.double(ALEAF),
					 as.double(ALEAFHM),
					 as.double(RD),
					 as.double(GAMMASTAR),
					 as.double(VJ),
					 as.double(KM),
					 as.double(VCMAX),
					 PACKAGE='GasExchangeR')
		 
vec <- Reduce(c,f)
names(vec) <- c('CI','PAR','TLEAF','PATM','JMAX25','IECO','EAVJ','EDVJ','DELSJ',
 'VCMAX25','EAVC','EDVC','DELSC','TVJUP','TVJDN','THETA','AJQ','GMESO','RD0','Q10F','RTEMP','DAYRESP',
 'TBELOW','HMSHAPE','AC','AJ','ALEAF','ALEAFHM', 'RD', 'GAMMASTAR','VJ','KM','VCMAX')

return(vec[returncols])
}

Aci <- function(returncols=c('CI','PAR','TLEAF','AC','AJ','ALEAF','ALEAFHM','RD','GAMMASTAR','VJ','KM','VCMAX'),...){
	oldform <- formals(acifun)
	formals(acifun)$returncols <- returncols
	dfr <- as.data.frame(t(mapply(acifun, ...)))
	formals(acifun) <- oldform
	return(dfr)
}

