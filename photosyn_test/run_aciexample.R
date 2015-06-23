
library(plantecophys)

dfr <- read.csv("photosyn_test/aciexample.csv")

c3 <- subset(dfr, Identity == 3)
f <- fitaci(c3)


fits <- fitacis(dfr, "Identity")



c3$Rd <- 2
f2 <- fitaci(c3, useRd=T)

g <- fitaci(c3, varnames=list(ALEAF="Photo", Tleaf="Tleaf", Ci="Ci", PPFD="PARi"))
