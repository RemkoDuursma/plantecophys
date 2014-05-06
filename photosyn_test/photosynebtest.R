

dfr <- expand.grid(VPD=seq(1,4,length=25), Wind=c(0.1,0.5,5,50))
x <- PhotosynEB(Wind=dfr$Wind, VPD=dfr$VPD)
x$Wind <- dfr$Wind
y <- Photosyn(VPD=dfr$VPD)


windows(7,7)
par(mar=c(5,5,2,2), cex.lab=1.2)
palette(colorRampPalette(c("red","blue"))(length(unique(x$Wind))))
with(x, plot(VPD, ALEAF/ELEAF, pch=19, col=as.factor(Wind),
             xlab="VPD (Air) (kPa)", ylab=expression(ITE~(A/E)~(mu*mol~mmol^-1))))
legend("topright",  levels(as.factor(x$Wind)), title="Wind speed", col=palette(),pch=19)
with(y[order(y$VPD),], points(VPD, ALEAF/ELEAF, type='l'))
legend("top", c("Photosyn","PhotosynEB"), lty=c(1,-1), pch=c(-1,19))

