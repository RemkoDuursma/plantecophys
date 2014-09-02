

windlabel <- expression(Wind~speed~~(m~s^-1))


library(magicaxis)
gslow <- 0.02
gshigh <- 0.5

Winds <- exp(seq(log(0.1), log(10), length=25))

rlow <- sapply(Winds, function(x)plantecophys:::FindTleaf(gs=gslow, Wind=x, Tair=25, penmon="iso"))
rhigh <- sapply(Winds, function(x)plantecophys:::FindTleaf(gs=gshigh, Wind=x, Tair=25, penmon="iso"))
rlow2 <- sapply(Winds, function(x)plantecophys:::FindTleaf(gs=gslow, Wind=x, Tair=25, penmon="full"))
rhigh2 <- sapply(Winds, function(x)plantecophys:::FindTleaf(gs=gshigh, Wind=x, Tair=25, penmon="full"))


plot(log10(Winds), rlow, ylim=c(20,30), type='l', col="red", axes=FALSE,
     xlab=windlabel, ylab=expression(T[leaf]~~(degree*C)))
points(log10(Winds), rhigh, type='l', col="blue")
# points(log10(Winds), rlow2, type='l', col="red", lty=5)
# points(log10(Winds), rhigh2, type='l', col="blue", lty=5)

abline(h=25)
magaxis(side=1, unlog=1)
axis(2)
box()
legend("topright", legend=sapply(c(bquote(g[s] == .(gslow)),
                          bquote(g[s] == .(gshigh)),
                          bquote(T[air])
                          ),as.expression), lty=1, col=c("red","blue","black"),
       cex=0.9)







f <- function(w,gs,penmon="full",...){
  
  tleaf <- plantecophys:::FindTleaf(gs=gs, Wind=w, Tair=25, penmon=penmon, VPD=1)
  flux <- plantecophys:::LeafEnergyBalance(Tleaf=tleaf, Wind=w, Tair=25, gs=gs, penmon=penmon, VPD=1,
                                           returnwhat="fluxes")
  flux$Tleaf <- tleaf

return(flux)
}


# Wind speed and E/gs
gss <- seq(0.02, 0.5, length=25)
windlow <- 0.1
windhigh <- 10

wlow <- do.call(rbind, lapply(gss, function(x)f(windlow,x)))
whigh <- do.call(rbind, lapply(gss, function(x)f(windhigh,x)))

plot(gss, wlow$ELEAFeb, type='l', col="red", ylim=c(0,5),xlim=c(0,0.6),
     xlab=expression(g[s]~~(mol~m^-2~s^-1)),
     ylab=expression(E[leaf]~~(mmol~m^-2~s^-1)))

points(gss, whigh$ELEAFeb, type='l', col="blue")

legend("topleft", legend=sapply(c(bquote(wind == .(windlow)),
                                   bquote(wind == .(windhigh))),as.expression), 
       lty=1, col=c("red","blue"), cex=0.9)




rlow <- do.call(rbind,lapply(Winds, function(x)f(x, gs=gslow, penmon="iso")))
rhigh <- do.call(rbind,lapply(Winds, function(x)f(x, gs=gshigh, penmon="iso")))
rlow2 <- do.call(rbind,lapply(Winds, function(x)f(x, gs=gslow, penmon="full")))
rhigh2 <- do.call(rbind,lapply(Winds, function(x)f(x, gs=gshigh, penmon="full")))


plot(log10(Winds), rlow$ELEAFeb/gslow/1000,  type='l', col="red", ylim=c(0,0.02),
     axes=FALSE,  xlab=windlabel, ylab=expression(E[leaf]/g[s]~~(mol~mol^-1)))
points(log10(Winds), rhigh$ELEAFeb/gshigh/1000, type='l', col="blue")
# points(log10(Winds), rlow2$ELEAFeb/gslow/1000, type='l', col="red", lty=5)
# points(log10(Winds), rhigh2$ELEAFeb/gshigh/1000, type='l', col="blue", lty=5)
abline(h=0.01)

magaxis(side=1, unlog=1)
axis(2)
box()
legend("topright", legend=sapply(c(bquote(g[s] == .(gslow)),
                                   bquote(g[s] == .(gshigh)),
                                   bquote(VPD)), as.expression),
       lty=1, col=c("red","blue","black"))


f <- lapply(Winds, function(x)FARAO2(Wind=x, energybalance=TRUE))
f <- do.call(rbind,f)












