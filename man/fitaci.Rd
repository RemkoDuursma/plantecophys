\name{fitaci}
\alias{fitaci}
\alias{plot.acifit}
\alias{print.acifit}
\alias{summary.acifit}
\title{Fit A-Ci curve}
\usage{
  fitaci(dat,
    varnames = list(ALEAF = "Photo", Tleaf = "Tleaf", Ci = "Ci", PPFD = "PPFD"),
    nlsmethod = "default")

  plot.acifit(x, ...)

  print.acifit(x, ...)

  summary.acifit(x, ...)
}
\arguments{
  \item{dat}{Dataset with Ci, Photo, Tleaf, PPFD (the last
  two are optional).}

  \item{varnames}{List of names of variables (see
  Details).}

  \item{nlsmethod}{Method passed to nls2 (!!not implemented
  yet)}
}
\description{
  Fits the Farquhar model of photosynthesis to measurements
  of photosynthesis and intercellular CO2 concentration
  (Ci). Estimates Jmax, Vcmax, Rd and their (approximate)
  standard errors.
}
\details{
  Fits the entire model at once (no splicing)... details to
  follow. nls2.
}

