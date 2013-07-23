# make figure



ci <- seq(50, 800, length=101)

acirun <- Aci(Ci=ci)


jpeg("man/figures/supplydemand.jpg", width=500, height=250, quality=100, pointsize=8)
par(xaxs="i", yaxs="i", cex.lab=2, mar=c(6,6,2,1), cex.axis=1.6)
with(acirun, plot(Ci, ALEAF, type='l',
                  ylim=c(0,20), xlim=c(0,800),
                  col="blue",lwd=2,
                  xlab=expression(italic(C)[i]~~(ppm)),
                  ylab=expression(italic(A)[leaf]~~(mu*mol~m^-2~s^-1))))

p <- Photosyn(VPD=2.5)

gc <- p$GS / 1.6
abline(gc*p$Ca, -gc, col="forestgreen", lwd=2)

points(p$Ci, p$ALEAF, pch=19, col="red", cex=3)

legend("bottomright", c(expression("Demand curve: "~A==f(italic(C)[i])),
                                   "Supply curve: "~A==g[CO2]*(italic(C)[a]-italic(C)[i]),
                                  "Operating point"),
       lwd=c(2,2,-1),pch=c(-1,-1,19), pt.cex=c(-1,-1,3),
       col=c("blue","forestgreen","red"), inset=0.01, cex=1.54, bty='n')
dev.off()

