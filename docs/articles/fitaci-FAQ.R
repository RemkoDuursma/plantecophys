## ----eval=FALSE----------------------------------------------------------
#  citation("plantecophys")

## ------------------------------------------------------------------------
packageVersion("plantecophys")

## ------------------------------------------------------------------------
library(plantecophys)

# Assume a mesophyll conductance of 0.2 mol m-2 s-1 bar-1
f <- fitaci(acidata1, gmeso=0.2)

## ------------------------------------------------------------------------
# Assume a mesophyll conductance of 0.2 mol m-2 s-1 bar-1
acidata1$Cc <- with(acidata1, Ci - Photo/0.2)

# Fit normally, but make sure to use Cc!
f <- fitaci(acidata1, varnames=list(ALEAF="Photo", Ci="Cc", Tleaf="Tleaf", PPFD="PARi"))

