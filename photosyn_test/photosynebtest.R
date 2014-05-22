

dfr <- expand.grid(VPD=seq(1,4,length=25), Wind=c(0.1,0.5,5,50))
x <- PhotosynEB(Wind=dfr$Wind, VPD=dfr$VPD)
y <- Photosyn(VPD=dfr$VPD)


windows(7,7)
par(mar=c(5,5,2,2), cex.lab=1.2)
palette(colorRampPalette(c("red","blue"))(length(unique(x$Wind))))
with(x, plot(VPD, ALEAF/ELEAF, pch=19, col=as.factor(Wind),
             xlab="VPD (Air) (kPa)", ylab=expression(ITE~(A/E)~(mu*mol~mmol^-1))))
legend("topright",  levels(as.factor(x$Wind)), title="Wind speed", col=palette(),pch=19)
with(y[order(y$VPD),], points(VPD, ALEAF/ELEAF, type='l'))
legend("top", c("Photosyn","PhotosynEB"), lty=c(1,-1), pch=c(-1,19))

windows(7,7)
par(mar=c(5,5,2,2), cex.lab=1.2)
with(x, plot(VPDleaf, ALEAF/ELEAF, pch=19, col=as.factor(Wind),
             xlab="VPD (Leaf) (kPa)", ylab=expression(ITE~(A/E)~(mu*mol~mmol^-1))))
legend("topright",  levels(as.factor(x$Wind)), title="Wind speed", col=palette(),pch=19)
with(y[order(y$VPD),], points(VPD, ALEAF/ELEAF, type='l'))
legend("top", c("Photosyn","PhotosynEB"), lty=c(1,-1), pch=c(-1,19))



# for Bonan et al. 2014
vcmaxs <- seq(50, 150, length=25)
jmaxs <- 2 * vcmaxs

run <- Photosyn(Vcmax=vcmaxs, Jmax=jmaxs, VPD=2, Tleaf=25, PPFD=1800, g1=3)
run2 <- FARAO(Vcmax=vcmaxs, Jmax=jmaxs, VPD=2, Tair=25, PPFD=1800, Wind=0.7, Wleaf=0.1, energybalance=TRUE)

with(run2, plot(GS, ALEAF, type='l', xlim=c(0,0.6), ylim=c(0,40)))
with(run, points(GS, ALEAF, type='l', col="blue"))






