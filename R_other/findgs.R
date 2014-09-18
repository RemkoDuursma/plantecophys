

findgs <- function(Tleaftarget, Tair, ...){
  
  
  O <- function(gs, Tleaftarget,Tair,...){
    
    Tl <- FindTleaf(gs, Tair,...)
    (Tl - Tleaftarget)^2
    
  }
  o <- optimize(O, c(0,80), Tleaftarget=Tleaftarget, Tair=Tair,...) 
  
return(o$minimum)
}



# T response at constant VPD
tairs <- seq(28, 34, length=20)

# Target is 3C below Tair
r1 <- sapply(tairs, function(x)findgs(x-3, x, PPFD=500))

# Target is 28
r2 <- sapply(tairs, function(x)findgs(28, x, PPFD=500))


windows(6,6)
par(mar=c(5,5,2,2), cex.lab=1.3)
plot(tairs, r1, type='l', ylim=c(0,2), col="forestgreen",
     xlab=expression(T[air]~~(degree*C)),
     ylab=expression(Stomatal~conductance~~(mol~m^-2~s^-1)))
text(29.5, 0.95, expression(Target == T[air] - 3~degree*C), cex=1.3, col="forestgreen")
points(tairs, r2, type='l', lty=5, col="darkorange")
text(32.8, 1.8, expression(Target == 28~degree*C), cex=1.3, col="darkorange")
dev.copy2pdf(file="cotton_tairgs.pdf")


# VPD response at constant T
vpds <- seq(1.5, 4, length=20)

# Tair is 30, target is 29
r1 <- sapply(vpds, function(x)findgs(29, 30, PPFD=500, VPD=x))
r2 <- sapply(vpds, function(x)findgs(29, 31, PPFD=500, VPD=x))


windows(6,6)
par(mar=c(5,5,2,2), cex.lab=1.3)
plot(vpds, r1, type='l', ylim=c(0,1), col="forestgreen",
     xlab=expression(VPD~~(kPa)),
     ylab=expression(Stomatal~conductance~~(mol~m^-2~s^-1)))
text(2, 0.12, expression(T[air] == 30~degree*C), cex=1.3, col="forestgreen")
points(vpds, r2, type='l', lty=5, col="darkorange")
text(2.5, 0.5, expression(T[air] == 31~degree*C), cex=1.3, col="darkorange")
legend("topright", expression(Target==29~degree*C), bty='n', cex=1.3)

dev.copy2pdf(file="cotton_vpdgs.pdf")


# PPFD response at constant VPD and Tair

ppfds <- seq(5, 800, length=25)

# Tair is 30, target is 29
r1 <- sapply(ppfds, function(x)findgs(29, 30, PPFD=x, VPD=2))
r2 <- sapply(ppfds, function(x)findgs(28, 30, PPFD=x, VPD=2))


windows(6,6)
par(mar=c(5,5,2,2), cex.lab=1.3)
plot(ppfds, r1, type='l', ylim=c(0,1), col="forestgreen",
     xlab=expression(PPFD~~(mu*mol~m^-2~s^-1)),
     ylab=expression(Stomatal~conductance~~(mol~m^-2~s^-1)))
text(600, 0.12, expression(Target == 29~degree*C), cex=1.3, col="forestgreen")
points(ppfds, r2, type='l', lty=5, col="darkorange")
text(600, 0.65, expression(Target == 28~degree*C), cex=1.3, col="darkorange")
legend("topright", expression(T[air]==30~degree*C), bty='n', cex=1.3)
dev.copy2pdf(file="cotton_ppfdgs.pdf")




# Illustration effect high boundary layer

windlabel <- expression(Wind~speed~~(m~s^-1))

gslow <- 0.05
gshigh <- 0.5

Winds <- exp(seq(log(0.1), log(10), length=25))

# Calculate ELEAF given known gs; calculate Tleaf from energy balance.
f <- function(w,gs,...){
  
  tleaf <- plantecophys:::FindTleaf(gs=gs, Wind=w, Tair=25,VPD=1, ...)
  flux <- plantecophys:::LeafEnergyBalance(Tleaf=tleaf, Wind=w, Tair=25, gs=gs, 
                                           VPD=1,returnwhat="fluxes",...)
  flux$Tleaf <- tleaf
  
  return(flux)
}

# Wind speed and E/gs
gss <- seq(0.02, 0.5, length=25)
windlow <- 0.1
windhigh <- 10

wlow <- do.call(rbind, lapply(gss, function(x)f(windlow,x, Wleaf=0.15)))
whigh <- do.call(rbind, lapply(gss, function(x)f(windhigh,x, Wleaf=0.05)))


windows(6,6)
par(mar=c(5,5,2,2), cex.lab=1.3)
plot(gss, wlow$ELEAFeb, type='l', col="red", ylim=c(0,5),xlim=c(0,0.6),
     xlab=expression(Stomatal~conductance~~(mol~m^-2~s^-1)),lwd=2,
     ylab=expression(Leaf~transpiration~~(mmol~m^-2~s^-1)))

points(gss, whigh$ELEAFeb, type='l', col="black", lwd=2)

legend("topleft", legend=c("Low (field)","High (cuvette)"),
       title="Boundary layer conductance",
       lty=1, col=c("red","black"), lwd=2, bty='n', cex=1.15)
dev.copy2pdf(file="gsEleaf_withgbl.pdf")



# Leaf temperature effects.
# Calculate Tleaf given known gs, from energy balance

library(magicaxis)
Winds <- exp(seq(log(0.1), log(10), length=25))

palette(colorRampPalette(c("blue","red"))(4))

r1 <- sapply(Winds, function(x)FindTleaf(gs=0.75, Wind=x, Tair=25))
r2 <- sapply(Winds, function(x)FindTleaf(gs=0.5, Wind=x, Tair=25))
r3 <- sapply(Winds, function(x)FindTleaf(gs=0.25, Wind=x, Tair=25))
r4 <- sapply(Winds, function(x)FindTleaf(gs=0.125, Wind=x, Tair=25))

windows(6,6)
par(mar=c(5,5,2,2), cex.lab=1.3)
plot(log10(Winds), r1, ylim=c(20,30), type='l', col=palette()[1], axes=FALSE,lwd=2,
     xlab=windlabel, ylab=expression(T[leaf]~~(degree*C)))
points(log10(Winds), r2, type='l', col=palette()[2],lwd=2)
points(log10(Winds), r3, type='l', col=palette()[3],lwd=2)
points(log10(Winds), r4, type='l', col=palette()[4],lwd=2)

abline(h=25)
magaxis(side=1, unlog=1)
axis(2)
box()
dev.copy2pdf(file="cotton_tleafwindgs.pdf")



gss <- seq(0.025, 1.5, length=101)
r <- sapply(gss, function(x)FindTleaf(gs=x, Wind=0.4, Wleaf=0.2, Tair=31, PPFD=500))

target <- 29
gsa <- findgs(target, 31,Wind=0.4, Wleaf=0.2, PPFD=500)

windows(6,6)
par(mar=c(5,5,2,2), cex.lab=1.3)
plot(gss, r, type='l', lwd=2, 
     xlab=expression(Stomatal~conductance~~(mol~m^-2~s^-1)),
     ylab=expression(Leaf~temperature~~(degree*C)),
     xlim=c(0,1.5), ylim=c(27,33))
abline(h=31, lty=5)
segments(x0=par()$usr[1], x1=gsa, y0=target, y1=target, col='blue', lwd=2)
arrows(x0=gsa, x1=gsa, y0=target, y1=par()$usr[3], length=0.2, col='blue', lwd=2)
text(1, 31+0.2, "Air temperature")
dev.copy2pdf(file="cotton_targetsolve.pdf")



