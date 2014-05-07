\name{RHtoVPD}
\alias{DewtoVPD}
\alias{RHairToLeaf}
\alias{RHleafToAir}
\alias{RHtoVPD}
\alias{VPDairToLeaf}
\alias{VPDleafToAir}
\alias{VPDtoRH}
\alias{esat}
\title{Conversions between relative humidity, vapour pressure deficit and dewpoint}
\usage{
RHtoVPD(RH, TdegC, Pa = 101)

VPDtoRH(VPD, TdegC, Pa = 101)

esat(TdegC, Pa = 101)

DewtoVPD(Tdew, TdegC, Pa = 101)

VPDleafToAir(VPD, Tleaf, Tair)

VPDairToLeaf(VPD, Tair, Tleaf)

RHleafToAir(RH, Tleaf, Tair)

RHairToLeaf(RH, Tair, Tleaf)
}
\arguments{
  \item{RH}{Relative humidity (\%)}

  \item{TdegC}{Temperature (degrees C) (either leaf or
  air)}

  \item{Tair}{Air temperature (degrees C)}

  \item{Tleaf}{Leaf temperature (degrees C)}

  \item{VPD}{Vapour pressure deficit (kPa)}

  \item{Pa}{Atmospheric pressure}

  \item{Tdew}{Dewpoint temperature}
}
\description{
Converts relative humidity (\%) to vapour pressure deficit
(kPa), or vice versa. To convert from relative humidity,
use the \code{RHtoVPD} function, use \code{VPDtoRH} for the
other way around. The water vapor saturation pressure is
calculated with \code{esat} (using Jones 1992). Also
provided is \code{DewtoVPD} to convert from dewpoint
temperature to VPD. The functions \code{VPDleafToAir} and
\code{VPDairToLeaf} convert VPD from a leaf temperature to
an air-temperature basis and vice versa.
}
\author{
Remko Duursma
}
\references{
Jones, H.G. 1992. Plants and microclimate: a quantitative
approach to environmental plant physiology. 2nd Edition.,
2nd Edn. Cambridge University Press, Cambridge. 428 p.
}

