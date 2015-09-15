
library(HIEv)
met <- downloadCSV(c("OUTMET","20131231"))
met <- subset(met, PAR > 100 & VPD > 0.8 & TARef < 30)
met$VPD <- with(met, RHtoVPD(RHref_al, TARef))
met <- subset(met, !is.na(VPD) & !is.na(PAR) & !is.na(TARef))
              

Lambda <- 0.002
n <- 50
ii <- sample(1:nrow(met),n)
runcon <- FARAO2(PPFD=met$PAR[ii], VPD=met$VPD[ii], Tleaf=met$TARef[ii],
                 energybalance=FALSE, lambda=Lambda)

bb <- function(dfr)with(dfr, ALEAF/(sqrt(VPD)*Ca))

run1 <- FARAO2(PPFD=met$PAR[ii], VPD=met$VPD[ii], Tair=met$TARef[ii],
              energybalance=TRUE, lambda=Lambda, Wind=0.35)

plot(bb(runcon), runcon$GS, ylim=c(0,0.25), xlim=c(0,0.04))
points(bb(run1), run1$GS, col="red")
abline(lm(runcon$GS ~ bb(runcon)),col="black")
abline(lm(run1$GS ~ bb(run1)),col="red")




