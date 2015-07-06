

df <- read.csv("http://files.figshare.com/1886204/WUEdatabase_merged_Lin_et_al_2015_NCC.csv")

.varnames <- list(ALEAF="Photo", GS="Cond",VPD="VPD", Ca="CO2S")

dr <- subset(df, Datacontrib == "John Drake")
fitBB(dr, varnames=.varnames)

pa <- subset(df, Species == "Phillyrea angustifolia")
fitBB(pa, varnames=.varnames, fitg0=T)

dfs <- split(df, df$Species)
l <- lapply(dfs, fitBB, varnames=.varnames)


fitBB(dr, varnames=.varnames)
f <- fitBB(dr, varnames=.varnames, fitg0=TRUE)


dr$RH <- VPDtoRH(dr$VPD, dr$Tleaf)
.varnames$RH <- "RH"
f <- fitBB(dr, varnames=.varnames, fitg0=TRUE, gsmodel="BallBerry")
