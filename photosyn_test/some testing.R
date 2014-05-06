
aci2 <- read.csv("photosyn_test/aci_court2.csv")
aci2$ID <- with(aci2, paste(plot,"-",pot))

# names(aci2)[names(aci2) == "Ps"] <- "Photo"

fits <- fitacis(aci2, group="ID",
                varnames=list(ALEAF="Ps", Tleaf="Temp", Ci="Ci", PPFD="PPFD"))
