\name{Aci}
\alias{Aci}
\alias{Photosyn,}
\title{Coupled leaf gas exchange model}
\usage{
  Aci(Ci, ...)
}
\arguments{
  \item{VPD}{Vapour pressure deficit (kPa)}

  \item{Ca}{Atmospheric CO2 concentration (ppm)}

  \item{PPFD}{Photosynthetic photon flux density ('PAR')
  (mu mol m-2 s-1)}

  \item{Tleaf}{Leaf temperature (degrees C)}

  \item{g0,g1}{Parameters of Medlyn et al. (2011) model of
  stomatal conductance (gs).}

  \item{gk}{Optional, exponent of VPD in gs model (Duursma
  et al. 2013)}

  \item{vpdmin}{Below vpdmin, VPD=vpdmin, to avoid very
  high gs.}

  \item{alpha}{Quantum yield of electron transport (mol
  mol-1)}

  \item{theta}{Shape of light response curve.}

  \item{Jmax}{Maximum rate of electron transport at 25
  degrees C (mu mol m-2 s-1)}

  \item{Vcmax}{Maximum carboxylation rate at 25 degrees C
  (mu mol m-2 s-1)}

  \item{Rd0}{Dark respiration rata at reference temperature
  (\code{TrefR})}

  \item{Q10}{Temperature sensitivity of Rd.}

  \item{TrefR}{Reference temperature for Rd.}

  \item{Rdayfrac}{Ratio of Rd in the light vs. in the
  dark.}

  \item{EaV,}{EdVC, delsC Vcmax temperature response
  parameters}

  \item{EaJ,}{EdVJ, delsJ Jmax temperature response
  parameters}

  \item{Ci}{Optional, intercellular CO2 concentration}

  \item{whichA}{Which assimilation rate does gs respond
  to?}
}
\value{
  A dataframe.
}
\description{
  Farquhar-Ball-Berry (Medlyn 2011 as default), with T
  response.
}
\details{
  If Ci is provided as an input, this function calculates
  an A-Ci curve. Otherwise, Ci is calculated from the
  intersection between the 'supply' and 'demand'
  relationships, using the stomatal conductance model of
  Medlyn et al. (2011).
}

