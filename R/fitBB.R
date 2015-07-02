fitBB <- function(df, varnames=list(ALEAF="Photo", GS="Cond", VPD="VpdL", Ca="CO2S",RH="RH"),
                  gsmodel=c("BBOpti","BBLeuning","BallBerry"),
                  fitg0=FALSE){
  
  gsmodel <- match.arg(gsmodel)
  
  gs <- df[,varnames$GS]  
  if(is.null(gs))stop("GS data missing - check varnames.")
  vpd <- df[,varnames$VPD]
  if(is.null(vpd))stop("VPD data missing - check varnames.")
  aleaf <- df[,varnames$ALEAF]  
  if(is.null(aleaf))stop("ALEAF data missing - check varnames.")
  ca <- df[,varnames$Ca]
  if(is.null(ca))stop("Ca data missing - check varnames.")
  if("RH" %in% names(varnames)){
    
    rh <- df[,varnames$RH]
    if(is.null(rh) & gsmodel == "BallBerry"){
      stop("To fit Ball-Berry you must first add RH to the dataset.")
    }
    if(max(rh, na.rm=TRUE) > 1){
      message("RH provided in % converted to relative units.")
      rh <- rh / 100
    }
  }
  
  if(gsmodel == "BBOpti"){
  if(!fitg0){
    fit <- nls(gs ~ 1.6*(1 + g1/VPD)*(aleaf/ca), start=list(g1=2)) 
  } else {
    fit <- nls(gs ~ g0 + 1.6*(1 + g1/VPD)*(aleaf/ca), start=list(g1=2)) 
  }
  }
  if(gsmodel == "BBLeuning"){
  if(!fitg0){
    fit <- nls(gs ~ 1.6*aleaf*g1/Ca/(1 + vpd/D0), start=list(g1=2, D0=1.5))
  } else {
    fit <- nls(gs ~ g0 + 1.6*aleaf*g1/Ca/(1 + vpd/D0), start=list(g1=2, D0=1.5))
  }
  }
  if(gsmodel == "BallBerry"){
    if(!fitg0){
      fit <- nls(gs ~ 1.6*g1*aleaf/(rh*Ca), start=list(g1=2))
    } else {
      fit <- nls(gs ~ g0 + 1.6*g1*aleaf/(rh*Ca), start=list(g1=2))
    }
  }

return(fit)  

}
