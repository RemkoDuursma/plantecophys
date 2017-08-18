

psil_e <- function(ELEAF, kl, psis){

  psil <- psis - ELEAF / kl
  
  psil[!is.finite(psil)] <- psis

psil
}

fsig_tuzet <- function(psil, sf=3.2, psif=-1.9){
  (1+exp(sf*psif))/(1+exp(sf*(psif-psil)))
}



PhotosynTuzet_f <- function(g1=4,
                          Ca=400,
                          psis=0, 
                          kl=2, 
                          sf=3, 
                          psif=-2,
                          ...){

  vn <- as.list(match.call())[-1]
  if("gsmodel" %in% vn){
    stop("Cannot define gsmodel with PhotosynTuzet.", .call=FALSE)
  }
  
  O <- function(psil, psis, kl, sf, psif, g1, Ca, ...){
    
    p <- Photosyn(g1=g1, Ca=Ca, gsmodel="BBdefine", BBmult=(g1/Ca)*fsig_tuzet(psil, sf, psif), ...)
    psilout <- psil_e(p$ELEAF, kl, psis)
    
    psil - psilout  # objective function: psil in = psil out.
  }
  
  topt <- uniroot(O, c(-20,0), psis=psis, kl=kl, sf=sf, psif=psif, g1=g1, Ca=Ca, ...)
  
  p <- Photosyn(g1=g1, Ca=Ca, gsmodel="BBdefine", BBmult=(g1/Ca)*fsig_tuzet(topt$root, sf, psif), ...)
  p <- cbind(p, data.frame(PSIL=topt$root))

return(p)
}

#' Coupled leaf gas exchange model with Tuzet stomatal conductance
#' @description An implementation of the coupled photosynthesis - stomatal conductance model, using the Tuzet et al. (2003) model of stomatal conductance. Accepts all arguments of \code{\link{Photosyn}} (except \code{gsmodel}, of course).
#' @param g1 The slope parameter. Note that the default value should be much higher than that used in the Medlyn et al. (2011) model to give comparable predictions.
#' @param Ca Atmospheric CO2 concentration.
#' @param psis Soil water potential (MPa). Note that soil-to-root hydraulic conductance is not implemented.
#' @param kl Leaf-specific hydraulic conductance (mmol m-2 s-1 MPa-1)
#' @param sf Shape parameter (-) of sigmoidal function of leaf water potential (see Tuzet et al. 2003)
#' @param psif Leaf water potential at which stomatal conductance is 50\% of maximum (MPa).
#' @param \dots All other arguments in \code{\link{Photosyn}}
#'@export
PhotosynTuzet <- function(g1=8,
                          Ca=400,
                          psis=0, 
                          kl=2, 
                          sf=3, 
                          psif=-2,
                          ...){
  as.data.frame(t(mapply(PhotosynTuzet_f, g1=g1, Ca=Ca, psis=psis, kl=kl, sf=sf, psif=psif, ...)))
}


