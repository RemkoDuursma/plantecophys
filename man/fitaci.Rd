\name{fitaci}
\alias{fitaci}
\alias{fitacis}
\alias{plot.acifit}
\title{Fit the Farquhar Berry von Caemmerer model of photosynthesis}
\usage{
  fitaci(dat,
    varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", Ci = "Ci", PPFD = "PPFD"),
    nlsmethod = "default", quiet = FALSE, ...)

  \method{plot}{acifit} (x, ...)

  fitacis(data, group, ...)
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
}
\value{
  A list of class 'acifit', with three components: 'df' is
  a dataframe with the original data, and the fitted
  photosynthetic rate (Amodel). 'pars' contains the
  parameter estimates and their approximate standard
  errors. 'nlsfit' is the object returned by
  \code{\link{nls2}}, and contains more detail on the
  quality of the fit.
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
}

