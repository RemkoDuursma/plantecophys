#' Fit the Farquhar-Berry-von Caemmerer model of leaf photosynthesis
#' @description Fits the Farquhar-Berry-von Caemmerer model of photosynthesis to measurements of 
#' photosynthesis and intercellular \eqn{CO_2}{CO2} concentration (Ci). Estimates Jmax, Vcmax, Rd 
#' and their standard errors. A simple plotting method is also included, as well as the function 
#' \code{\link{fitacis}} which quickly fits multiple A-Ci curves (see its help page). Temperature 
#' dependencies of the parameters are taken into account following Medlyn et al. (2002), see 
#' \code{\link{Photosyn}} for more details.
#' @param data Dataframe with Ci, Photo, Tleaf, PPFD (the last two are optional). For \code{fitacis}, 
#' also requires a grouping variable.
#' @param varnames List of names of variables in the dataset (see Details).
#' @param Tcorrect If TRUE, Vcmax and Jmax are corrected to 25C. Otherwise, Vcmax and Jmax 
#' are estimated at measurement temperature.
#' @param Patm Atmospheric pressure (kPa)
#' @param citransition If provided, fits the Vcmax and Jmax limited regions separately (see Details).
#' @param quiet If TRUE, no messages are written to the screen.
#' @param startValgrid If TRUE (the default), uses a fine grid of starting values to increase 
#' the chance of finding a solution.
#' @param fitmethod Method to fit the A-Ci curve. Either 'default' (Duursma 2015), 'bilinear' (See Details), or 'onepoint' (De Kauwe et al. 2016).
#' @param algorithm Passed to \code{\link{nls}}, sets the algorithm for finding parameter values.
#' @param fitTPU Logical (default FALSE). Attempt to fit TPU limitation (fitmethod set to 'bilinear' 
#' automatically if used). See Details.
#' @param alphag When estimating TPU limitation (with \code{fitTPU}), an 
#' additional parameter (see Details).
#' @param useRd If Rd provided in data, and useRd=TRUE (default is FALSE), uses measured Rd in 
#' fit. Otherwise it is estimated from the fit to the A-Ci curve.
#' @param PPFD Photosynthetic photon flux density ('PAR') (mu mol m-2 s-1)
#' @param Tleaf Leaf temperature (degrees C)
#' @param alpha Quantum yield of electron transport (mol mol-1)
#' @param theta Shape of light response curve.
#' @param gmeso Mesophyll conductance (mol m-2 s-1 bar-1). If not NULL (the default), Vcmax 
#' and Jmax are chloroplastic rates.
#' @param EaV,EdVC,delsC Vcmax temperature response parameters
#' @param EaJ,EdVJ,delsJ Jmax temperature response parameters
#' @param Km,GammaStar Optionally, provide Michaelis-Menten coefficient for Farquhar model, and 
#' Gammastar. If not provided, they are calculated with a built-in function of leaf temperature.
#' @param id Names of variables (quoted, can be a vector) in the original dataset to be stored
#'  in the result. Most useful when using \code{\link{fitacis}}, see there for examples of its use.
#' @param \dots Further arguments (ignored at the moment).
#'
#'
#' @details 
#' 
#' \subsection{Fitting method}{The default method to fit A-Ci curves (set by 
#' \code{fitmethod="default"}) uses non-linear regression to fit the A-Ci curve. No assumptions 
#' are made on which part of the curve is Vcmax or Jmax limited. Normally, all three parameters 
#' are estimated: Jmax, Vcmax and Rd, unless Rd is provided as measured (when \code{useRd=TRUE}, 
#' and Rd is contained in the data). This is the method as described by Duursma (2015, Plos One).
#' 
#' The 'bilinear' method to fit A-Ci curves (set by \code{fitmethod="bilinear"}) linearizes 
#' the Vcmax and Jmax-limited regions, and applies linear regression twice to estimate first 
#' Vcmax and Rd, and then Jmax (using Rd estimated from the Vcmax-limited region). The transition 
#' point is found as the one which gives the best overall fit to the data (i.e. all possible 
#' transitions are tried out, similar to Gu et al. 2010, PCE). The advantage of this method is that 
#' it \emph{always} returns parameter estimates, so it should be used in cases where the default 
#' method fails. Be aware, though, that the default method fails mostly when the curve contains 
#' bad data (so check your data before believing the fitted parameters).
#' 
#' When \code{citransition} is set, it splits the data into a Vcmax-limited (where Ci < citransition), 
#' and Jmax-limited region (Ci > citransition). Both parameters are then estimated separately for 
#' each region (Rd is estimated only for the Vcmax-limited region). \bold{Note} that the actual 
#' transition point as shown in the standard plot of the fitted A-Ci curve may be quite different 
#' from that provided, since the fitting method simply decides which part of the dataset to use 
#' for which limitation, it does not constrain the actual estimated transition point directly. 
#' See the example below. If \code{fitmethod="default"}, it applies non-linear regression to 
#' both parts of the data, and when fitmethod="bilinear", it uses linear regression on the 
#' linearized photosynthesis rate. Results will differ between the two methods (slightly).}
#' 
#' The 'onepoint' fitting method is a very simple estimation of Vcmax and Jmax for every point 
#' in the dataset, simply by inverting the photosynthesis equation. See De Kauwe et al. (2016) 
#' for details. The output will give the original data with Vcmax and Jmax added (note you can 
#' set \code{Tcorrect} as usual!). For increased reliability, this method only works if 
#' dark respiration (Rd) is included in the data (\code{useRd} is set automatically when 
#' setting \code{fitmethod='one-point'}). This method is not recommended for full A-Ci curves, 
#' but rather for spot gas exchange measurements, when a simple estimate of Vcmax or Jmax 
#' is needed, for example when separating stomatal and non-stomatal drought effects on 
#' photosynthesis (Zhou et al. 2013, AgForMet). The user will have to decide whether the Vcmax 
#' or Jmax rates are used in further analyses. This fitting method can not be used in \code{fitacis}, because Vcmax and Jmax are already estimated for each point in the dataset.
#' 
#' \subsection{TPU limitation}{Optionally, the \code{fitaci} function estimates the triose-phosphate 
#' utilization (TPU) rate. The TPU can act as another limitation on photosynthesis, and can be 
#' recognized by a 'flattening out' of the A-Ci curve at high Ci. When \code{fitTPU=TRUE}, the 
#' fitting method used will always be 'bilinear'. The TPU is estimated by trying out whether the 
#' fit improves when the last n points of the curve are TPU-limited (where n=1,2,...). When TPU is 
#' estimated, it is possible (though rare) that no points are Jmax-limited (in which case estimated 
#' Jmax will be NA). A minimum of two points is always reserved for the estimate of Vcmax and Rd. 
#' An additional parameter (\code{alphag}) can be set that affects the behaviour at high Ci (see 
#' Ellsworth et al. 2015 for details, and also \code{\link{Photosyn}}). See examples.}
#' 
#' \subsection{Temperature correction}{When \code{Tcorrect=TRUE} (the default), Jmax and Vcmax 
#' are re-scaled to 25C, using the temperature response parameters provided (but Rd is always 
#' at measurement temperature). When \code{Tcorrect=FALSE}, estimates of all parameters are at 
#' measurement temperature. If TPU is fit, it is never corrected for temperature. Important 
#' parameters to the fit are GammaStar and Km, both of which are calculated from leaf temperature 
#' using standard formulations. Alternatively, they can be provided as known inputs.}
#' 
#' \subsection{Mesophyll conductance}{It is possible to provide an estimate of the mesophyll 
#' conductance as input (\code{gmeso}), in which case the fitted Vcmax and Jmax are to be interpreted 
#' as chloroplastic rates. When using gmeso, it is recommended to use the 'default' fitting 
#' method (which will use the Ethier&Livingston equations inside \code{Photosyn}). It is also 
#' implemented with the 'bilinear' method but it requires more testing (and seems to give 
#' some strange results). When gmeso is set to a relatively low value, the resulting fit 
#' may be quite strange.}
#' 
#' \subsection{Other parameters}{The A-Ci curve parameters depend on the values of a number 
#' of other parameters. For Jmax, PPFD is needed in order to express it as the asymptote. If 
#' PPFD is not provided in the dataset, it is assumed to equal 1800 mu mol m-2 s-1 (in which 
#' case a warning is printed). It is possible to either provide PPFD as a variable in the 
#' dataset (with the default name 'PARi', which can be changed), or as an argument to 
#' the \code{fitaci} directly.}
#' 
#' \subsection{Plotting and summarizing}{The default \strong{plot} of the fit is constructed 
#' with \code{\link{plot.acifit}}, see Examples below. When plotting the fit, the A-Ci curve 
#' is simulated using the \code{\link{Aci}} function, with leaf temperature (Tleaf) and PPFD 
#' set to the mean value for the dataset. The \strong{coefficients} estimated in the fit 
#' (Vcmax, Jmax, and usually Rd) are extracted with \code{coef}. The summary of the fit is 
#' the same as the 'print' method, that is \code{myfit} will give the same output as 
#' \code{summary(myfit)} (where \code{myfit} is an object returned by \code{fitaci}).
#' 
#' Because fitaci returns the fitted \code{\link{nls}} object, more details on statistics of 
#' the fit can be extracted with standard tools. The Examples below shows the use of the 
#' \pkg{nlstools} to extract many details of the fit at once. The fit also includes the 
#' \strong{root mean squared error} (RMSE), which can be extracted as \code{myfit$RMSE}. This 
#' is a useful metric to compare the different fitting methods.}
#' 
#' \subsection{Atmospheric pressure correction}{Note that atmospheric pressure (Patm) is taken 
#' into account, assuming the original data are in molar units (Ci in mu mol mol-1, or ppm). 
#' During the fit, Ci is converted to mu bar, and Km and Gammastar are recalculated accounting 
#' for the effects of Patm on the partial pressure of oxygen. When plotting the fit, though,
#'  molar units are shown on the X-axis. Thus, you should get (nearly) the same fitted curve when
#'   Patm was set to a value lower than 100kPa, but the fitted Vcmax and Jmax will be higher. 
#'   This is because at low Patm, photosynthetic capacity has to be higher to achieve the same 
#'   measured photosynthesis rate.}
#' 
#' @troubleshooting From time to time, the \code{fitaci} function returns an error, indicating 
#' that the A-Ci curve could not be fit. In the majority of cases, this indicates a bad curve 
#' that probably could not (and should not) be fit in any case. Inspect the raw data to check 
#' if the curve does not include severe outliers, or large deviations from the expected, typical 
#' A-Ci curve. If the curve looks fine, refit using the option \code{fitmethod="bilinear"}, 
#' which will always return estimated parameters. Whatever you do, do not trust fitted curves 
#' without inspecting the raw data, and the fit of the model to the data.
#' 
#' 
#' @references
#' Duursma, R.A., 2015. Plantecophys - An R Package for Analysing and Modelling Leaf Gas 
#' Exchange Data. PLoS ONE 10, e0143346. doi:10.1371/journal.pone.0143346
#' 
#' De Kauwe, M. G. et al. 2016. A test of the 'one-point method' for estimating maximum carboxylation capacity from field-measured, light-saturated photosynthesis. New Phytol 210, 1130-1144.
#' 
#' @return A list of class 'acifit', with the following components:
#' \describe{
#' \item{df}{A dataframe with the original data, including the measured photosynthetic 
#' rate (Ameas), the fitted photosynthetic rate (Amodel), Jmax and Vcmax-limited gross 
#' rates (Aj, Ac), TPU-limited rate (Ap), dark respiration (Rd), leaf temperature (Tleaf), 
#' chloroplastic CO2 (Cc), PPFD, atmospheric pressure (Patm), and 'original Ci, i.e. the 
#' Ci used as input (which is different from the Ci used in fitting if Patm was not set to 100kPa)}
#' \item{pars}{Contains the parameter estimates and their approximate standard errors}
#' \item{nlsfit}{The object returned by \code{\link{nls}}, and contains more detail on 
#' the quality of the fit}
#' \item{Tcorrect}{whether the temperature correction was applied (logical)}
#' \item{Photosyn}{A copy of the \code{\link{Photosyn}} function with the arguments adjusted for 
#' the current fit. That is, Vcmax, Jmax and Rd are set to those estimated in the fit, and Tleaf and 
#' PPFD are set to the mean value in the dataset. All other parameters that were set in fitaci are 
#' also used (e.g. temperature dependency parameters, TPU, etc.).}
#' \item{Ci_transition}{The Ci at which photosynthesis transitions from Vcmax 
#' to Jmax limited photosynthesis.}
#' \item{Rd_measured}{Logical - was Rd provided as measured input?}
#' \item{GammaStar}{The value for GammaStar, either calculated or provided to the fit.}
#' \item{Km}{he value for Km, either calculated or provided to the fit.}
#' \item{kminput}{Was Km provided as input? (If FALSE, it was calculated from Tleaf)}
#' \item{gstarinput}{Was GammaStar provided as input? (If FALSE, it was calculated from Tleaf)}
#' \item{fitmethod}{The fitmethod uses, either default or bilinear}
#' \item{citransition}{The input citransition (NA if it was not provided as input)}
#' \item{gmeso}{The mesophyll conductance used in the fit (NA if it was not set)}
#' \item{fitTPU}{Was TPU fit?}
#' \item{alphag}{The value of alphag used in estimating TPU.}
#' \item{RMSE}{The Root-mean squared error, calculated as \code{sqrt(sum((Ameas-Amodel)^2))}.}
#' \item{runorder}{The data returned in the 'df' slot are ordered by Ci, but in rare cases the 
#' original order of the data contains information; 'runorder' is 
#' the order in which the data were provided.}
#' 
#' }
#' @examples
#' \dontrun{
#' # Fit an A-Ci curve on a dataframe that contains Ci, Photo and optionally Tleaf and PPFD. 
#' # Here, we use the built-in example dataset 'acidata1'.
#' f <- fitaci(acidata1)
#' 
#' # Note that the default behaviour is to correct Vcmax and Jmax for temperature, 
#' # so the estimated values are at 25C. To turn this off:
#' f2 <- fitaci(acidata1, Tcorrect=FALSE)
#' 
#' # To use different T response parameters (see ?Photosyn),
#' f3 <- fitaci(acidata1, Tcorrect=TRUE, EaV=25000)
#' 
#' # Make a standard plot
#' plot(f)
#' 
#' # Look at a summary of the fit
#' summary(f)
#' 
#' # Extract coefficients only
#' coef(f)
#' 
#' # The object 'f' also contains the original data with predictions.
#' # Here, Amodel are the modelled (fitted) values, Ameas are the measured values.
#' with(f$df, plot(Amodel, Ameas))
#' abline(0,1)
#' 
#' # The fitted values can also be extracted with the fitted() function:
#' fitted(f)
#' 
#' # The non-linear regression (nls) fit is stored as well,
#' summary(f$nlsfit)
#' 
#' # Many more details can be extracted with the nlstools package
#' library(nlstools)
#' overview(f$nlsfit)
#'  
#' # The curve generator is stored as f$Photosyn:
#' # Calculate photosynthesis at some value for Ci, using estimated 
#' # parameters and mean Tleaf, PPFD for the dataset.
#' f$Photosyn(Ci=820)
#' 
#' # Photosynthetic rate at the transition point:
#' f$Photosyn(Ci=f$Ci_transition)$ALEAF
#' 
#' # Set the transition point; this will fit Vcmax and Jmax separately. Note that the *actual* 
#' # transition is quite different from that provided, this is perfectly fine : 
#' # in this case Jmax is estimated from the latter 3 points only (Ci>800), but the actual 
#' # transition point is at ca. 400ppm.
#' g <- fitaci(acidata1, citransition=800)
#' plot(g)
#' g$Ci_transition
#' 
#' # Use measured Rd instead of estimating it from the A-Ci curve. 
#' # The Rd measurement must be added to the dataset used in fitting, 
#' # and you must set useRd=TRUE.
#' acidata1$Rd <- 2
#' f2 <- fitaci(acidata1, useRd=TRUE)
#' f2
#' 
#' # Fit TPU limitation
#' ftpu <- fitaci(acidata1, fitTPU=TRUE, PPFD=1800, Tcorrect=TRUE)
#' plot(ftpu)
#' }
#' @export
#' @rdname fitaci
#' @importFrom stats lm
fitaci <- function(data, 
                   varnames=list(ALEAF="Photo", Tleaf="Tleaf", 
                                 Ci="Ci", PPFD="PARi", Rd="Rd"),
                   Tcorrect=TRUE, 
                   Patm=100,
                   citransition=NULL,
                   quiet=FALSE, startValgrid=TRUE, 
                   fitmethod=c("default","bilinear","onepoint"),
                   algorithm="default", 
                   fitTPU=FALSE,
                   alphag=0,
                   
                   useRd=FALSE,
                   PPFD=NULL,
                   Tleaf=NULL,
                   
                   alpha=0.24,
                   theta=0.85,
                   gmeso=NULL,
                   EaV = 82620.87,
                   EdVC = 0,
                   delsC = 645.1013,
                   EaJ = 39676.89,
                   EdVJ = 200000,
                   delsJ = 641.3615,
                   
                   GammaStar = NULL,
                   Km = NULL,
                   id=NULL,
                   ...){
  
  fitmethod <- match.arg(fitmethod)
  if(fitTPU & fitmethod == "default")fitmethod <- "bilinear"
  
  if(fitTPU & fitmethod == "onepoint"){
    Stop("TPU limitation not implemented with one-point method.")
  }
  
  if(fitmethod == "onepoint")useRd <- TRUE
  
  gstarinput <- !is.null(GammaStar)
  kminput <- !is.null(Km)
  if(is.null(gmeso))gmeso <- -999  # cannot pass NULL value to nls
  
  # Check data 
  if(nrow(data) == 0){
    Stop("No rows in data - check observations.")
  }
  
  # Make sure data is a proper dataframe, not some tibble or other nonsense.
  data <- as.data.frame(data)
  
  # Extra parameters not used at the moment; warn when provided.
  chkDots(...)
  
  # Add PPFD and Tleaf to data, if needed (uses default values, or input values)
  data <- set_PPFD(varnames, data, PPFD, quiet)
  data <- set_Tleaf(varnames, data, Tleaf, quiet)
  
  # Set measured Rd if provided (or warn when provided but not used)
  Rd_meas <- set_Rdmeas(varnames, data, useRd, citransition, quiet)
  haveRd <- !is.na(Rd_meas)
  
  # Must have Rd with one-point method.
  if(!haveRd && fitmethod == "onepoint"){
    Stop("Must add Rd to the dataset (and check that varnames is set correctly).")
  }
  
  # Extract Ci and apply pressure correction
  data$Ci_original <- data[,varnames$Ci]
  data$Ci <- data[,varnames$Ci] * Patm/100
  
  # Extract measured net leaf photosynthesis
  data$ALEAF <- data[,varnames$ALEAF]
  
  # Calculate Km and GammaStar, if not input 
  # (Photosyn does this as well, but have to repeat it here 
  # since we cannot have NULL in call to nls)
  Km_v <- if(!kminput) TKm(data$Tleaf, Patm) else rep(Km, nrow(data))
  GammaStar_v <- if(!gstarinput){
    TGammaStar(data$Tleaf, Patm)
  } else {
    rep(GammaStar, nrow(data))
  }
  
  # Citransition not defined:
  # default method = full non-linear model
  # bilinear = optimize transition point in bilinear fit
  if(is.null(citransition)){
    
    if(fitmethod == "default"){
    
      f <- do_fit_method1(data, haveRd, Rd_meas, Patm, startValgrid, 
                          Tcorrect, algorithm,alpha,theta,gmeso,EaV,
                          EdVC,delsC,EaJ,EdVJ,delsJ,GammaStar_v,Km_v)
    } 
    if(fitmethod == "bilinear"){
      
      f <- do_fit_method_bilinear_bestcitrans(data, haveRd, fitTPU, alphag, Rd_meas, 
                                              Patm, Tcorrect, algorithm,
                                              alpha,theta,gmeso,EaV,EdVC,delsC,
                                              EaJ,EdVJ,delsJ,
                                              GammaStar_v, Km_v)
    }
    if(fitmethod == "onepoint"){
      f <- do_fit_method_bilinear(data, haveRd, alphag, Rd_meas, Patm, 
                                  NA,NA,  # don't matter
                                  Tcorrect=Tcorrect, algorithm,
                                  alpha,theta,gmeso,EaV,EdVC,delsC,
                                  EaJ,EdVJ,delsJ,
                                  GammaStar, Km, onepoint=TRUE)
      return(f)
    }
  }
  
  # Ci transition defined. 
  # default - two nonlinear fits (first Vcmax, then Jmax region)
  # bilinear - two linear fits
  if(!is.null(citransition)){
    
    # NOTE-- bug in default method when citransition specified. Switch off for now.
    fitmethod <- "bilinear"
    
    # if(fitmethod == "default"){
    #   f <- do_fit_method2(data, haveRd, Rd_meas, Patm, citransition, Tcorrect, 
    #                       algorithm,alpha,theta,gmeso,EaV,EdVC,
    #                       delsC,EaJ,EdVJ,delsJ,GammaStar_v,Km_v)
    # }
    if(fitmethod == "bilinear"){
      f <- do_fit_method_bilinear(data, haveRd, alphag, Rd_meas, Patm, citransition, 
                                  NULL, Tcorrect, algorithm,
                                alpha,theta,gmeso,EaV,EdVC,delsC,EaJ,EdVJ,delsJ,
                                GammaStar_v, Km_v)
    }
    
  }
  
  
  # TPU. If it is not estimated, put something in because we need it later.
  if(!("TPU" %in% names(f)))f$TPU <- 1000
  
  # Only used to add 'Amodel' to the output
  acirun <- do_acirun(data,f,Patm,Tcorrect,  # Note several parameters are stored in 'f'
                      alpha=alpha,theta=theta,
                      gmeso=gmeso,EaV=EaV,
                      EdVC=EdVC,delsC=delsC,
                      EaJ=EaJ,EdVJ=EdVJ,
                      delsJ=delsJ,GammaStar=GammaStar_v, Km=Km_v)

  # If Ap is never actually limiting, set estimated TPU to 1000
  # This means there is no evidence for TPU-limitation.
  if(!any(acirun$Ap < acirun$Aj))f$TPU <- 1000
  
  # Organize output
  l <- list()  
  runorder <- order(acirun$Ci)
  l$df <- acirun[runorder,]
  if(!is.null(id)){
    if(any(!id %in% names(data))){
      Warning("id ignored - not all variables are in dataset provided")
    } else {
      l$df <- cbind(l$df, data[id])
    }
  }
  l$id <- id
  l$pars <- f$pars
  if(fitTPU){
    val <- ifelse(f$TPU < 1000, f$TPU, NA)
    tpu <- matrix(c(val,NA), ncol=2)
    colnames(tpu) <- colnames(l$pars)
    rownames(tpu) <- "TPU"
    l$pars <- rbind(l$pars, tpu)
  }
  
  l$nlsfit <- f$fit
  l$Tcorrect <- Tcorrect
  
  # Save function itself, the formals contain the parameters 
  # used to fit the A-Ci curve.
  # First save Tleaf, PPFD in the formals (as the mean of the dataset)
  formals(Photosyn)$Tleaf <- mean(data$Tleaf)
  formals(Photosyn)$Patm <- Patm
  formals(Photosyn)$PPFD <- mean(data$PPFD)
  formals(Photosyn)$Vcmax <- l$pars[1]
  formals(Photosyn)$Jmax <- l$pars[2]
  formals(Photosyn)$Rd <- l$pars[3]
  formals(Photosyn)$TPU <- f$TPU
  formals(Photosyn)$alphag <- alphag
  formals(Photosyn)$Tcorrect <- Tcorrect
  formals(Photosyn)$alpha <- alpha
  formals(Photosyn)$theta <- theta
  formals(Photosyn)$gmeso <- gmeso
  formals(Photosyn)$EaV <- EaV
  formals(Photosyn)$EdVC <- EdVC
  formals(Photosyn)$delsC <- delsC
  formals(Photosyn)$EaJ <- EaJ
  formals(Photosyn)$EdVJ <- EdVJ
  formals(Photosyn)$delsJ <- delsJ
  
  if(gstarinput)formals(Photosyn)$GammaStar <- GammaStar
  if(kminput)formals(Photosyn)$Km <- Km
  
  l$Photosyn <- Photosyn
  
  # Store Ci at which photosynthesis transitions from 
  # Jmax to Vcmax limitation
  l$Ci_transition <- findCiTransition(l$Photosyn)
  l$Rd_measured <- haveRd
  
  # Save GammaStar and Km (either evaluated at mean temperature, 
  # or input if provided)
  l$GammaStar <- mean(GammaStar_v)
  l$Km <- mean(Km_v)
  l$kminput <- kminput
  l$gstarinput <- gstarinput
  
  l$fitmethod <- fitmethod
  l$citransition <- ifelse(is.null(citransition), NA, citransition)
  l$gmeso <- ifelse(is.null(gmeso) || gmeso < 0, NA, gmeso)
  l$fitTPU <- fitTPU
  l$alphag <- alphag
  l$RMSE <- rmse_acifit(l)
  l$runorder <- runorder
  
  class(l) <- "acifit"
  
return(l)
}


#-- Component functions of fitaci

rmse_acifit <- function(x)sqrt(sum((x$df$Ameas - x$df$Amodel)^2))

guess_Jmax <- function(data, Rd_guess, Patm, Tcorrect){
  
  # Guess Jmax from max A, T-corrected gammastar (starting value)
  maxCi <- max(data$Ci)
  mi <- which.max(data$Ci)

  maxPhoto <- data$ALEAF[mi]

  Tl <- data$Tleaf[mi]
  gammastar <- TGammaStar(Tl,Patm)
  VJ <- (maxPhoto+Rd_guess) / ((maxCi - gammastar) / (maxCi + 2*gammastar))
  Jmax_guess <- VJ*4
  if(Tcorrect){
    Teffect <- TJmax(Tl,  EaJ=39676.89, delsJ=641.3615, EdVJ=200000)
    Jmax_guess <- Jmax_guess / Teffect
  }
  
  return(Jmax_guess)
}



guess_Vcmax <- function(data, Jmax_guess, Rd_guess, Patm, Tcorrect){
  # Guess Vcmax, from section of curve that is definitely Vcmax-limited
  dato <- data[data$Ci < 150 & data$Ci > 60 & data$ALEAF > 0,]
  if(nrow(dato) > 0){
    Km <- TKm(dato$Tleaf,Patm)
    gammastar <- TGammaStar(dato$Tleaf,Patm)
    vcmax <- with(dato, (ALEAF + Rd_guess) / ((Ci - gammastar)/(Ci + Km)))
    Vcmax_guess <- median(vcmax)
  } else {
    Vcmax_guess <- Jmax_guess/1.8 
  }
  if(Tcorrect){
    Teffect <- TVcmax(mean(data$Tleaf), EaV=82620.87, delsC=645.1013, EdVC=0)
    Vcmax_guess <- Vcmax_guess / Teffect
  }
  return(Vcmax_guess)
}


guess_Rd <- function(haveRd, Rd_meas, default=1.5){
  if(haveRd){
    Rd_guess <- Rd_meas
  } else {
    Rd_guess <- default
  }  
  return(Rd_guess)
}

set_PPFD <- function(varnames, data, PPFD, quiet){
  
  if(!varnames$PPFD %in% names(data) & is.null(PPFD)){
    data$PPFD <- 1800
    if(!quiet)Warning("PARi not in dataset; assumed PARi = 1800.")
  } else {
    if(!is.null(PPFD))
      data$PPFD <- PPFD
    else
      data$PPFD <- data[,varnames$PPFD]
  }
  return(data)
}



set_Tleaf <- function(varnames, data, Tleaf, quiet){
  # Check if Tleaf is provided
  if(!varnames$Tleaf %in% names(data) & is.null(Tleaf)){
    data$Tleaf <- 25
    if(!quiet)Warning("Tleaf not in dataset; assumed Tleaf = 25.")
  } else {
    if(!is.null(Tleaf))
      data$Tleaf <- Tleaf
    else
      data$Tleaf <- data[,varnames$Tleaf]
  }
return(data)
}



# Check if Rd is provided and whether we want to set it as known.
set_Rdmeas <- function(varnames, data, useRd, citransition, quiet){
  
  Rd_meas <- NA
  if(!is.null(varnames$Rd)){ 
    if(varnames$Rd %in% names(data) && useRd){
      
      # Has to be a single unique value for this dataset
      Rd_meas <- data[,varnames$Rd]
      Rd_meas <- unique(Rd_meas)
      if(length(Rd_meas) > 1){
        Stop(paste("If Rd provided as measured, it must be a single",
                   "unique value for an A-Ci curve."))
      }
      
      # Use positive value throughout.
      Rd_meas <- abs(Rd_meas)
      haveRd <- TRUE
      
      if(!is.null(citransition))
        Stop("At the moment cannot provide citransition as well as measured Rd.")
    }
    if(varnames$Rd %in% names(data) && !useRd){
      if(!quiet)
        message(paste("Rd found in dataset but useRd set to FALSE.",
                      "Set to TRUE to use measured Rd."))
    }
  }
  return(Rd_meas)
}





do_fit_method1 <- function(data, haveRd, Rd_meas, Patm, startValgrid, 
                           Tcorrect, algorithm,
                           alpha,theta,gmeso,EaV,EdVC,delsC,EaJ,EdVJ,delsJ,
                           GammaStar,Km){
  
  # Guess Rd (starting value)
  Rd_guess <- guess_Rd(haveRd, Rd_meas)
  Jmax_guess <- guess_Jmax(data, Rd_guess, Patm, Tcorrect)
  Vcmax_guess <- guess_Vcmax(data, Jmax_guess, Rd_guess, Patm, Tcorrect)
  
  # Fine-tune starting values; try grid of values around initial estimates.
  if(startValgrid){
    if(!haveRd){
      aciSS <- function(Vcmax, Jmax, Rd){
        Photo_mod <- acifun_wrap(data$Ci, PPFD=data$PPFD, 
                                 Vcmax=Vcmax, Jmax=Jmax, 
                                 Rd=Rd, Tleaf=data$Tleaf, Patm=Patm,
                                 TcorrectVJ=Tcorrect,
                                 alpha=alpha,theta=theta,
                                 gmeso=gmeso,EaV=EaV,
                                 EdVC=EdVC,delsC=delsC,
                                 EaJ=EaJ,EdVJ=EdVJ,
                                 delsJ=delsJ,Km=Km,GammaStar=GammaStar)
        
        SS <- sum((data$ALEAF - Photo_mod)^2)
        return(SS)
      }
      d <- 0.4
      n <- 25
      gg <- expand.grid(Vcmax=seq(Vcmax_guess*(1-d),Vcmax_guess*(1+d),length=n),
                        Rd=seq(Rd_guess*(1-d),Rd_guess*(1+d),length=n))
      
      m <- with(gg, mapply(aciSS, Vcmax=Vcmax, Jmax=Jmax_guess, Rd=Rd))
      ii <- which.min(m)
      Vcmax_guess <- gg$Vcmax[ii]
      Rd_guess <- gg$Rd[ii]
      
    } else {
      
      aciSS <- function(Vcmax, Jmax, Rd){
        Photo_mod <- acifun_wrap(data$Ci, PPFD=data$PPFD, 
                                 Vcmax=Vcmax, Jmax=Jmax, 
                                 Rd=Rd, Tleaf=data$Tleaf, Patm=Patm,
                                 TcorrectVJ=Tcorrect,
                                 alpha=alpha,theta=theta,
                                 gmeso=gmeso,EaV=EaV,
                                 EdVC=EdVC,delsC=delsC,
                                 EaJ=EaJ,EdVJ=EdVJ,
                                 delsJ=delsJ,Km=Km,GammaStar=GammaStar)
        SS <- sum((data$ALEAF - Photo_mod)^2)
        return(SS)
      }
      d <- 0.3
      n <- 20
      gg <- data.frame(Vcmax=seq(Vcmax_guess*(1-d),Vcmax_guess*(1+d),length=n))
      
      m <- with(gg, mapply(aciSS, Vcmax=Vcmax, Jmax=Jmax_guess, Rd=Rd_meas))
      ii <- which.min(m)
      Vcmax_guess <- gg$Vcmax[ii]
      
    }
    
  }
  
  if(!haveRd){
    # Fit Vcmax, Jmax and Rd
    nlsfit <- try(nls(ALEAF ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, 
                                      Jmax=Jmax, Rd=Rd, Tleaf=Tleaf, Patm=Patm,
                                      TcorrectVJ=Tcorrect,
                                      alpha=alpha,theta=theta,
                                      gmeso=gmeso,EaV=EaV,
                                      EdVC=EdVC,delsC=delsC,
                                      EaJ=EaJ,EdVJ=EdVJ,
                                      delsJ=delsJ,Km=Km,GammaStar=GammaStar),
                  algorithm=algorithm,
                  data=data, control=nls.control(maxiter=500, minFactor=1/10000),
                  start=list(Vcmax=Vcmax_guess, Jmax=Jmax_guess, Rd=Rd_guess)), silent=TRUE)
    
    if(inherits(nlsfit, "try-error")){
      Stop("Could not fit curve - check quality of data or fit using fitmethod='bilinear'.")
    }
    
    p <- coef(nlsfit)
    pars <- summary(nlsfit)$coefficients[,1:2]
  } else {
    
    # Fit Vcmax and Jmax
    nlsfit <- try(nls(ALEAF ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, 
                                      Jmax=Jmax, Rd=Rd_meas, Tleaf=Tleaf, Patm=Patm,
                                      TcorrectVJ=Tcorrect,
                                      alpha=alpha,theta=theta,
                                      gmeso=gmeso,EaV=EaV,
                                      EdVC=EdVC,delsC=delsC,
                                      EaJ=EaJ,EdVJ=EdVJ,
                                      delsJ=delsJ,Km=Km,GammaStar=GammaStar),
                  algorithm=algorithm,
                  data=data, control=nls.control(maxiter=500, minFactor=1/10000),
                  start=list(Vcmax=Vcmax_guess, Jmax=Jmax_guess)), silent=TRUE)
    
    if(inherits(nlsfit, "try-error")){
      Stop(paste("Could not fit curve -",
                 "check quality of data or fit using fitmethod='bilinear'."))
    }
    
    p <- coef(nlsfit)
    p[[3]] <- Rd_meas
    names(p)[3] <- "Rd"
    
    pars <- summary(nlsfit)$coefficients[,1:2]
    pars <- rbind(pars, c(Rd_meas, NA))
    rownames(pars)[3] <- "Rd"
  }
  
  
  return(list(pars=pars, fit=nlsfit))
}



# do_fit_method2 <- function(data, haveRd, Rd_meas, Patm, citransition, 
#                            Tcorrect,algorithm,alpha,theta,gmeso,EaV,EdVC,
#                            delsC,EaJ,EdVJ,delsJ,GammaStar,Km){
#   
#   # Guess Rd (starting value)
#   Rd_guess <- guess_Rd(haveRd, Rd_meas)
#   Jmax_guess <- guess_Jmax(data, Rd_guess, Patm, Tcorrect)
#   Vcmax_guess <- guess_Vcmax(data, Jmax_guess, Rd_guess, Patm, Tcorrect)
#   
#   # If citransition provided, fit twice.
#   dat_vcmax <- data[data$Ci < citransition,]
#   dat_jmax <- data[data$Ci >= citransition,]
#   
#   if(nrow(dat_vcmax) > 0){
#     nlsfit_vcmax <- nls(ALEAF ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=Vcmax, 
#                                             Jmax=10^6, Rd=Rd, Tleaf=mean(Tleaf), 
#                                             Patm=Patm, returnwhat="Ac",
#                                             TcorrectVJ=Tcorrect,
#                                             alpha=alpha,theta=theta,
#                                             gmeso=gmeso,EaV=EaV,
#                                             EdVC=EdVC,delsC=delsC,
#                                             EaJ=EaJ,EdVJ=EdVJ,
#                                             delsJ=delsJ,Km=Km,GammaStar=GammaStar),
#                         algorithm=algorithm,
#                         data=dat_vcmax, 
#                         control=nls.control(maxiter=500, minFactor=1/10000),
#                         start=list(Vcmax=Vcmax_guess, Rd=Rd_guess))
#     p1 <- coef(nlsfit_vcmax)
#   } else {
#     nlsfit_vcmax <- NULL
#     p1 <- c(Vcmax=NA, Rd=NA)
#   }
#   
#   # Don't re-estimate Rd; it is better estimated from the Vcmax-limited region.
#   Rd_vcmaxguess <- if(!is.null(nlsfit_vcmax))coef(nlsfit_vcmax)[["Rd"]] else 1.5
#   
#   if(nrow(dat_jmax) > 0){
#     nlsfit_jmax <- nls(ALEAF ~ acifun_wrap(Ci, PPFD=PPFD, Vcmax=10000, 
#                                            Jmax=Jmax, Rd=Rd_vcmaxguess, 
#                                            Tleaf=mean(Tleaf), Patm=Patm, returnwhat="Aj",
#                                            TcorrectVJ=Tcorrect,
#                                            alpha=alpha,theta=theta,
#                                            gmeso=gmeso,EaV=EaV,
#                                            EdVC=EdVC,delsC=delsC,
#                                            EaJ=EaJ,EdVJ=EdVJ,
#                                            delsJ=delsJ,Km=Km,GammaStar=GammaStar),
#                        algorithm=algorithm,
#                        data=dat_jmax, control=nls.control(maxiter=500, minFactor=1/10000),
#                        start=list(Jmax=Jmax_guess))
#     p2 <- coef(nlsfit_jmax)
#   } else { 
#     nlsfit_jmax <- NULL
#     p2 <- c(Jmax=NA)
#   }
#   
#   p <- c(p1[1],p2,p1[2])
#   pars1 <- if(!is.null(nlsfit_vcmax)){
#     summary(nlsfit_vcmax)$coefficients[,1:2] 
#   } else {
#     matrix(rep(NA,4),ncol=2)
#   }
#   
#   pars2 <- if(!is.null(nlsfit_jmax)){
#     summary(nlsfit_jmax)$coefficients[,1:2] 
#   } else {
#     matrix(rep(NA,2),ncol=2)
#   }
#   
#   pars <- rbind(pars1[1,],pars2,pars1[2,])
#   rownames(pars) <- c("Vcmax","Jmax","Rd")
#   nlsfit <- list(nlsfit_vcmax=nlsfit_vcmax,nlsfit_jmax=nlsfit_jmax)
#   
#   return(list(pars=pars, fit=nlsfit))      
# }



do_fit_method_bilinear <- function(data, haveRd, alphag, Rd_meas, 
                                   Patm, citransition, citransition2=NULL,
                                   Tcorrect, algorithm,
                                   alpha,theta,gmeso,EaV,EdVC,delsC,
                                   EaJ,EdVJ,delsJ,
                                   GammaStar, Km, onepoint=FALSE){
  
  # Calculate T-dependent parameters
  ppar <- Photosyn(Tleaf=data$Tleaf, Patm=Patm, Tcorrect=Tcorrect,
                   alpha=alpha,theta=theta,
                   gmeso=gmeso,EaV=EaV,
                   EdVC=EdVC,delsC=delsC,
                   EaJ=EaJ,EdVJ=EdVJ,
                   delsJ=delsJ,
                   returnParsOnly=TRUE)
  if(!is.null(GammaStar))ppar$GammaStar <- GammaStar
  if(!is.null(Km))ppar$Km <- Km
  if(is.null(citransition2))citransition2 <- 10^6
  
  # Calculate Cc if gmeso included
  if(!is.null(gmeso) && gmeso > 0){
    Conc <- data$Ci - data$ALEAF/gmeso
  } else {
    Conc <- data$Ci
  }
  
  # Linearize
  data$vcmax_pred <- (Conc - ppar$GammaStar)/(Conc + ppar$Km)
  data$Jmax_pred <- (Conc - ppar$GammaStar)/(Conc + 2*ppar$GammaStar)
  data$TPU_part <- (Conc - ppar$GammaStar)/(Conc - (1+3*alphag)*ppar$GammaStar)
    
  if(onepoint){
    
    orig_dat <- data[, 1:(match("Ci_original", names(data))-1)]
    
    orig_dat$Vcmax <- (data$ALEAF + Rd_meas)  / data$vcmax_pred
    orig_dat$Jmax <- (data$ALEAF + Rd_meas)  / data$Jmax_pred
    
    if(Tcorrect){
      orig_dat$Jmax <- orig_dat$Jmax / TJmax(mean(data$Tleaf), EaJ, delsJ, EdVJ)
      orig_dat$Vcmax <- orig_dat$Vcmax / TVcmax(mean(data$Tleaf),EaV, delsC, EdVC)
    }
    
    return(orig_dat)
  }
  
  # Fit Vcmax and Rd from linearized portion
  datv <- data[data$Ci < citransition & data$Ci < citransition2,]
  if(nrow(datv) == 0){
    return(list(pars=NA, fit=NA, TPU=NA, success=FALSE))
  }
  
  if(!haveRd){
    fitv <- lm(ALEAF ~ vcmax_pred, data=datv)
    Rd_fit <- coef(fitv)[[1]]
    Vcmax_fit <- coef(fitv)[[2]]
  } else {
    # If using measured Rd, add to Anet, and remove intercept from fit.
    datv$ALEAFg <- datv$ALEAF + Rd_meas
    fitv <- lm(ALEAFg ~ vcmax_pred-1, data=datv)
    Rd_fit <- -Rd_meas
    Vcmax_fit <- coef(fitv)[[1]]
  }
  
  # Fit Jmax from linearized portion
  datj <- data[data$Ci >= citransition & data$Ci < citransition2,]
  datp <- data[data$Ci >= citransition2,]

  # Manual fix: if only one point for TPU, and none for Jmax, abandon fit.
  # In this case it would be more defensible to use the single point for Jmax.
  if(nrow(datp) == 1 && nrow(datj) == 0){
    return(list(pars=NA, fit=NA, TPU=NA, success=FALSE))
  }
  
  
  if(nrow(datj) > 0){
    
    # Fit gross photo using fitted Rd
    # (Rd_fit is negative)
    datj$Agross <- datj$ALEAF - Rd_fit
    
    # One point, calculate directly
    if(nrow(datj) == 1){
      J_fit <- with(datj, 4 * Agross / Jmax_pred)
    } else {
      fitj <- lm(Agross ~ Jmax_pred-1, data=datj)
      J_fit <- 4 * coef(fitj)[[1]]
    }
    
    # And solve for Jmax from inverse non-rect. hyperbola
    Jmax_fit <- inverseJfun(mean(data$PPFD), alpha, J_fit, theta)
    
  } else {
    Jmax_fit <- 10^6  # not elegant but will do for now 
                      # (avoids trouble elsewhere)
  }
  
  # TPU
  if(nrow(datp) > 0){
    datp$Agross <- datp$ALEAF - Rd_fit
    
    tpu_vals <- (1/3) * datp$Agross / datp$TPU_part
    
    TPU <- mean(tpu_vals)
  } else {
    TPU <- 1000  # same as default in Photosyn
  }
  
  # The above estimates are at the measured Tleaf.
  # Express at 25C?
  if(Tcorrect){
    Jmax_fit <- Jmax_fit / TJmax(mean(data$Tleaf), EaJ, delsJ, EdVJ)
    Vcmax_fit <- Vcmax_fit / TVcmax(mean(data$Tleaf),EaV, delsC, EdVC)
  }
  
  ses <- summary(fitv)$coefficients[,2]
  if(!haveRd){
    pars <- matrix(c(Vcmax_fit, Jmax_fit, -Rd_fit,
                   ses[2],NA,ses[1]), ncol=2)
  } else {
    pars <- matrix(c(Vcmax_fit, Jmax_fit, -Rd_fit,
                     ses[1],NA,NA), ncol=2)
  }
  
  
  rownames(pars) <- c("Vcmax","Jmax","Rd")
  colnames(pars) <- c("Estimate","Std. Error")
  
  
  return(list(pars=pars, fit=fitv, TPU=TPU, success=TRUE))
}

do_fit_method_bilinear_bestcitrans <- function(data, haveRd, fitTPU, alphag, 
                                               Rd_meas, Patm, Tcorrect,
                                               algorithm,alpha,theta,gmeso,EaV,
                                               EdVC,delsC,EaJ,EdVJ,
                                               delsJ,GammaStar, Km){
  
  # Possible Ci transitions
  ci <- data$Ci
  nci <- length(ci)
  citransitions <- diff(ci)/2 + ci[-nci]
  
  # at least two Ci values to estimate Vcmax and Rd, so delete first
  citransitions1 <- citransitions[-1]
  
  # Possible transitions to TPU
  if(!fitTPU){
    citransitions2 <- max(ci) + 1  # outside range, on purpose
  } else {
    # start at top, all the way down, leave lowest 2 points alone
    citransitions2 <- c(max(ci) + 1, rev(citransitions1))
  }
  citransdf <- expand.grid(ci1=citransitions1, ci2=citransitions2)
  citransdf <- citransdf[citransdf$ci1 <= citransdf$ci2,]
  SS <- c()
  
  # Note that Tcorrect is set to FALSE inside the loop. If Tcorrect is needed, it is done
  # after the loop finishes (avoids a bug).
  for(i in seq_len(nrow(citransdf))){
    
    fit <- do_fit_method_bilinear(data, haveRd, alphag, Rd_meas, Patm, 
                                  citransdf$ci1[i], citransdf$ci2[i], 
                                  Tcorrect=FALSE, algorithm,
                                  alpha,theta,gmeso,EaV,EdVC,delsC,
                                  EaJ,EdVJ,delsJ,
                                  GammaStar, Km)
    
    if(fit$success && !any(is.na(fit$pars[,"Estimate"]))){
        run <- do_acirun(data,fit,Patm,Tcorrect=FALSE,
                         alpha=alpha,theta=theta,
                         gmeso=gmeso,EaV=EaV,
                         EdVC=EdVC,delsC=delsC,
                         EaJ=EaJ,EdVJ=EdVJ,
                         delsJ=delsJ,GammaStar=GammaStar,Km=Km)
        
        SS[i] <- sum((run$Ameas - run$Amodel)^2)  
    } else {
      SS[i] <- 10^6
    }
  }
  
  # Best Ci transitions
  bestcis <- citransdf[which.min(SS),]
  
  f <- do_fit_method_bilinear(data, haveRd, alphag, Rd_meas, Patm, 
                              bestcis$ci1, bestcis$ci2, Tcorrect, algorithm,
                              alpha,theta,gmeso,EaV,EdVC,delsC,EaJ,EdVJ,delsJ,
                              GammaStar, Km)
  
  if(f$pars["Jmax","Estimate"] < 0){
    Stop("Cannot invert light response curve to estimate Jmax - increase alpha or theta.")
  }

  
  return(f)  
}

# Wrapper around Photosyn; this wrapper will be sent to nls. 
acifun_wrap <- function(Ci,..., TcorrectVJ, returnwhat="ALEAF"){
  r <- Photosyn(Ci=Ci,Tcorrect=TcorrectVJ,...)
  if(returnwhat == "ALEAF")return(r$ALEAF)
  if(returnwhat == "Ac")return(r$Ac - r$Rd)
  if(returnwhat == "Aj")return(r$Aj - r$Rd)
}


# Using fitted coefficients, get predictions from model.
do_acirun <- function(data,f,Patm,Tcorrect,...){
  
  acirun <- Photosyn(Ci=data$Ci, 
                     Vcmax=f$pars[1,1], Jmax=f$pars[2,1], 
                     Rd=f$pars[3,1], 
                     TPU=f$TPU,
                     PPFD=data$PPFD, 
                     Tleaf=data$Tleaf,
                     Patm=Patm,
                     Tcorrect=Tcorrect,...)
  
  acirun$Ameas <- data$ALEAF
  acirun$ELEAF <- NULL
  acirun$GS <- NULL
  acirun$Ca <- NULL
  acirun$Ci_original <- data$Ci_original
  names(acirun)[names(acirun) == "ALEAF"] <- "Amodel"
  
  # shuffle
  avars <- match(c("Ci","Ameas","Amodel"),names(acirun))
  acirun <- acirun[,c(avars, setdiff(1:ncol(acirun), avars))]
  
  return(acirun)
}


