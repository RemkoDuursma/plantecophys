\name{FARAO}
\alias{FARAO}
\title{FARquhar And Opti}
\usage{
FARAO(lambda = 0.002, Ca = 400, VPD = 1, gbl = NA, photo = c("BOTH",
  "VCMAX", "JMAX"), C4 = FALSE, ...)
}
\arguments{
  \item{lambda}{The marginal cost of water (mol mol-1)}

  \item{Ca}{The CO2 concentration.}

  \item{gbl}{Optional, the boundary layer conductance (mol
  m-2 s-1). Somewhat experimental!}

  \item{photo}{Which photosynthesis rate should stomata
  respond to? Defaults to 'BOTH', i.e. the minimum of Vcmax
  and Jmax limited rates.}

  \item{C4}{If TRUE, uses the C4 photosynthesis routine
  (\code{\link{AciC4}})}

  \item{...}{All other parameters are passed to
  \code{\link{Aci}}}
}
\description{
The numerical solution of the optimal stomatal conductance
model, coupled with the Farquhar model of photosynthesis.
The model of Medlyn et al. (2011) is an approximation to
this full numeric solution.
}
\details{
This model finds the Ci that maximizes A - lambda*E (Cowan
& Farquhar 1977, see also Medlyn et al. 2011).
}
\author{
Remko Duursma
}
\references{
Cowan, I. and G.D. Farquhar. 1977. Stomatal function in
relation to leaf metabolism and environment. Symposia of
the Society for Experimental Biology. 31:471-505.

Medlyn, B.E., R.A. Duursma, D. Eamus, D.S. Ellsworth, I.C.
Prentice, C.V.M. Barton, K.Y. Crous, P. De Angelis, M.
Freeman and L. Wingate. 2011. Reconciling the optimal and
empirical approaches to modelling stomatal conductance.
Global Change Biology. 17:2134-2144.
}

