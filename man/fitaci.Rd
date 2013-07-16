\name{fitaci}
\alias{fitaci}
\alias{fitacis}
\alias{plot.acifit}
\alias{print.acifit}
\alias{summary.acifit}
\title{Fit A-Ci curve}
\usage{
  fitaci(dat,
    varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", Ci = "Ci", PPFD = "PPFD"),
    nlsmethod = "default", quiet = FALSE)

  \method{plot}{acifit} (x, ...)

  \method{print}{acifit} (x, ...)

  summary.acifit(x, ...)

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
\description{
  Fits the Farquhar model of photosynthesis to measurements
  of photosynthesis and intercellular \eqn{CO_2}{CO2}
  concentration (Ci). Estimates Jmax, Vcmax, Rd and their
  (approximate) standard errors.
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

