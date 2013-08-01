\name{fitaci}
\alias{fitaci}
\alias{fitacis}
\alias{plot.acifit}
\alias{plot.acifits}
\title{Fit the Farquhar Berry von Caemmerer model of photosynthesis}
\usage{
  fitaci(dat,
    varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", Ci = "Ci", PPFD = "PARi"),
    nlsmethod = "default", quiet = FALSE, ...)

  \method{plot}{acifit} (x, what = c("data", "model"),
    xlim = NULL, ylim = NULL,
    whichA = c("Ac", "Aj", "Amin"), add = FALSE, pch = 19,
    addzeroline = TRUE, addlegend = !add,
    transitionpoint = TRUE, ...)

  fitacis(data, group, progressbar = TRUE...)

  \method{plot}{acifits} (x,
    how = c("manyplots", "oneplot"), ...)
}
\arguments{
  \item{dat}{Dataset with Ci, Photo, Tleaf, PPFD (the last
  two are optional).}

  \item{varnames}{List of names of variables (see
  Details).}

  \item{nlsmethod}{Method passed to nls2 ('algorithm' in
  nls2).}

  \item{quiet}{If TRUE, no messages are written to the
  screen.}

  \item{group}{For batch analysis using \code{fitacis}, the
  name of the grouping variable in the dataframe.}

  \item{what}{The default is to plot both the data and the
  model fit, or specify 'data' or 'model' to plot one of
  them.}

  \item{add}{If TRUE, adds to the current plot}

  \item{pch}{The plotting symbol for the data}

  \item{xlim}{Limits for the X axis, if left blank
  estimated from data}

  \item{ylim}{Limits for the Y axis, if left blank
  estimated from data}

  \item{whichA}{By default all three photosynthetic rates
  are plotted (Aj=Jmax-limited (blue), Ac=Vcmax-limited
  (red), Hyperbolic minimum (black)). Or, specify one or
  two of them.}

  \item{addzeroline}{If TRUE, the default, adds a dashed
  line at y=0}

  \item{addlegend}{If TRUE, adds a legend (by default does
  not add a legend if add=TRUE)}

  \item{how}{If 'manyplots', produces a single plot for
  each A-Ci curve. If 'oneplot' overlays all of them.}
}
\value{
  A list of class 'acifit', with five components:
  \describe{ \item{df}{A dataframe with the original data,
  the fitted photosynthetic rate (Amodel), Jmax and
  Vcmax-limited gross rates (Aj, Ac)} \item{pars}{Contains
  the parameter estimates and their approximate standard
  errors} \item{nlsfit}{The object returned by
  \code{\link{nls2}}, and contains more detail on the
  quality of the fit} \item{Photosyn}{A copy of the
  \code{\link{Photosyn}} function with the arguments
  adjusted for the current fit. That is, Vcmax, Jmax and Rd
  are set to those estimated in the fit, and Tleaf and PPFD
  are set to the mean value in the dataset.}
  \item{Ci_transition}{The Ci at which photosynthesis
  transitions from Vcmax to Jmax limited photosynthesis.} }
}
\description{
  Fits the Farquhar model of photosynthesis to measurements
  of photosynthesis and intercellular \eqn{CO_2}{CO2}
  concentration (Ci). Estimates Jmax, Vcmax, Rd and their
  (approximate) standard errors. Temperature dependencies
  are taken into account, see \code{\link{Photosyn}}.
}
\details{
  Uses non-linear regression to fit an A-Ci curve. No
  assumptions are made on which part of the curve is Vcmax
  or Jmax limited. Three parameters are estimated, Jmax,
  Vcmax (both at 25deg C) and Rd (at the measurement
  temperature).

  When plotting the fit, the A-Ci curve is simulated using
  the \code{\link{Aci}} function, with leaf temperature
  (Tleaf) and PPFD set to the mean value for the dataset.
  If PPFD is not provided in the dataset, it is assumed to
  equal 1800 mu mol m-2 s-1.
}
\examples{
# Fit an A-Ci curve on a dataframe that contains Ci, Photo and optionally Tleaf and PPFD. Here, we use the built-in example dataset 'acidata1'.
f <- fitaci(acidata1)

# Make a standard plot
plot(f)

# Look at a summary of the fit
summary(f)

# Extract coefficients only
coef(f)

# The object 'f' also contains the original data with predictions.
# Here, Amodel are the modelled (fitted) values, Ameas are the measured values.
with(f$df, plot(Amodel, Ameas))
abline(0,1)

# The fitted values can also be extracted with the fitted() function:
fitted(f)

# The non-linear regression (nls) fit is stored as well,
summary(f$nlsfit)

# The curve generator is stored as f$Photosyn:
# Calculate photosynthesis at some value for Ci, using estimated parameters and mean Tleaf, PPFD for the dataset
f$Photosyn(Ci=820)

# Photosynthetic rate at the transition point:
f$Photosyn(Ci=f$Ci_transition)$ALEAF
}

