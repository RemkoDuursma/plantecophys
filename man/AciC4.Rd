\name{AciC4}
\alias{AciC4}
\title{C4 Photosynthesis}
\usage{
AciC4(Ci, PPFD = 1500, Tleaf = 25, VPMAX25 = 120, VCMAX25 = 60,
  Vpr = 80, alpha = 0, gbs = 0.003, O2 = 210, JMAX25 = 400, x = 0.4,
  THETA = 0.7, Q10 = 2.3, RD0 = 1, RTEMP = 25, TBELOW = 0,
  DAYRESP = 1, Q10F = 2, FRM = 0.5)
}
\arguments{
  \item{Ci}{}

  \item{PPFD}{}

  \item{Tleaf}{}

  \item{VPMAX25}{}

  \item{VCMAX25}{}

  \item{Vpr}{PEP regeneration (mu mol m-2 s-1)}

  \item{alpha}{Fraction of PSII activity in the bundle
  sheath}

  \item{gbs}{Bundle shead conductance (mol m-2 s-1)}

  \item{O2}{Mesophyll O2 concentration}

  \item{JMAX25}{}

  \item{x}{Partitioning factor for electron transport}

  \item{THETA}{Shape parameter of the non-rectangular
  hyperbola}

  \item{Q10}{}

  \item{RD0}{}

  \item{RTEMP}{}

  \item{TBELOW}{}

  \item{DAYRESP}{}

  \item{Q10F}{}

  \item{FRM}{Fraction of dark respiration that is mesophyll
  respiration (Rm)}
}
\description{
An implementation of the A-Ci curve for C4 plants, based on
von Caemmerer et al. (1999)
}
\author{
Rhys Whitley
}

