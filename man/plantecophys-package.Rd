\name{plantecophys-package}
\alias{plantecophys-package}
\alias{plantecophys}
\docType{package}
\title{
\packageTitle{plantecophys}
}
\description{
\packageDescription{plantecophys}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{plantecophys}
\packageIndices{plantecophys}

The following functions are the main tools in \pkg{plantecophys}:

\enumerate{
  \item \code{\link{fitaci}} (and \code{\link{fitaci}}) fits A-Ci curves to data.
  \item \code{\link{Photosyn}} can be used to simulate A-Ci curves (or \code{\link{Aci}}), and simulate from a coupled leaf gas exchange model.
  \item \code{\link{fitBB}} fits Ball-Berry-type stomatal conductance models to data.
  \item \code{\link{FARAO}} is an implementation of a numeric solution to Cowan-Farquhar optimization of stomatal conductance.
  \item \code{\link{RHtoVPD}} converts relative humidity to vapour pressure deficit (and more similar functions on that help page).
}
 
The package also includes the following example datasets:

\enumerate{
  \item \code{\link{acidata1}} A dataset with a single A-Ci curve.
  \item \code{\link{manyacidat}} A dataset with many A-Ci curves.
}

}
\author{
Remko Duursma

Maintainer: Remko Duursma
}
\references{
Duursma, R.A., 2015. Plantecophys - An R Package for Analysing and Modelling Leaf Gas Exchange Data. PLoS ONE 10, e0143346. doi:10.1371/journal.pone.0143346
}
\keyword{ package }
