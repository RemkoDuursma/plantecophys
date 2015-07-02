#' Fit Ball-Berry type models of stomatal conductance
#' @description Fits one of three versions of the Ball-Berry type stomatal conductance models to 
#' observations of stomatal conductance (gs), photosynthesis (A), atmospheric CO2 concentration (Ca) 
#' and vapour pressure deficit (VPD).
#' @param df Input dataframe, containing all variables needed to fit the model.
#' @param varnames List of names of variables in the input dataframe. Relative humidity (RH) is only 
#' needed when the original Ball-Berry model is to be fit.
#' @param gsmodel One of BBOpti (Medlyn et al. 2011), BBLeuning (Leuning 1995), or BallBerry (Ball and Berry 1987).
#' @param fitg0 If TRUE, also fits the intercept term (g0, the 'residual conductance'). Default is FALSE.
#' @export
#' @references 
#' #' Leuning, R. 1995. A critical-appraisal of a combined stomatal-photosynthesis model for C-3 plants. Plant Cell and Environment. 18:339-355.
#'
#' Medlyn, B.E., R.A. Duursma, D. Eamus, D.S. Ellsworth, I.C. Prentice, C.V.M. Barton, K.Y. Crous, P. De Angelis, M. Freeman and L. Wingate. 2011. Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology. 17:2134-2144.
#' @importFrom stats nls
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
    fit <- try(nls(gs ~ 1.6*(1 + g1/vpd)*(aleaf/ca), start=list(g1=4)) )
  } else {
    fit <- try(nls(gs ~ g0 + 1.6*(1 + g1/vpd)*(aleaf/ca), start=list(g1=4, g0=0.005)) )
  }
  }
  if(gsmodel == "BBLeuning"){
  if(!fitg0){
    fit <- try(nls(gs ~ 1.6*aleaf*g1/Ca/(1 + vpd/D0), start=list(g1=4, D0=1.5)))
  } else {
    fit <- try(nls(gs ~ g0 + 1.6*aleaf*g1/Ca/(1 + vpd/D0), start=list(g1=4, D0=1.5, g0=0.005)))
  }
  }
  if(gsmodel == "BallBerry"){
    if(!fitg0){
      fit <- try(nls(gs ~ 1.6*g1*aleaf/(rh*Ca), start=list(g1=4)))
    } else {
      fit <- try(nls(gs ~ g0 + 1.6*g1*aleaf/(rh*Ca), start=list(g1=4, g0=0.005)))
    }
  }

l <- list()
l$gsmodel <- gsmodel
l$varnames <- varnames
l$fitg0 <- fitg0
l$data <- df
l$success <- !inherits(fit, "try-error")
l$coef <- if(l$success)coef(fit) else NA
l$fit <- fit
l$n <- length(residuals(fit))
class(l) <- "BBfit"

return(l)  
}






