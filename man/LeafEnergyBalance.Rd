\name{LeafEnergyBalance}
\alias{LeafEnergyBalance}
\title{Coupled leaf gas exchange model with energy balance}
\usage{
LeafEnergyBalance(Tleaf = 21.5, Tair = 20, gs = 0.15, PPFD = 1500,
  VPD = 2, Patm = 101, Wind = 2, Wleaf = 0.02, StomatalRatio = 1,
  LeafAbs = 0.86, returnwhat = c("balance", "fluxes"))
}
\arguments{
  \item{Wind}{Wind speed (m s-1)}

  \item{Wleaf}{Leaf width (m)}

  \item{StomatalRatio}{The stomatal ratio (cf. Licor6400
  terminology), if it is 1, leaves have stomata only on one
  side (hypostomatous), 2 for leaves with stomata on both
  sides (amphistomatous)}

  \item{LeafAbs}{Leaf absorptance of solar radiation}
}
\description{
As \code{\link{Photosyn}}, but calculates the leaf
temperature based on the leaf's energy balance. Including
sensible and long-wave heat loss, latent heat loss from
evaporation, and solar radiation input.
}
\details{
Uses the Penman-Monteith equation to calculate the leaf
transpiration rate, and finds Tleaf by solving the leaf
energy balance iteratively. In the solution, it is
accounted for that stomatal conductance and net radiation
depend on Tleaf. There are no simplifications or
approximations made to find the solution.
}

