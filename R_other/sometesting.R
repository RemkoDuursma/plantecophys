


# dA/dE
# I'd like to plot dA/dE vs. gs (cf. Buckley 2014), but OPTfun does not do this
# -> Photosyn has to take GS as input, solve for Ci and A (quadratic easily set up)


getdAdE <- function(Ci,...){
  delta <- 1e-03
  r1 <- Photosyn(Ci=Ci, ...)
  r2 <- Photosyn(Ci=Ci+delta, ...)
  dA <- r2$ALEAF - r1$ALEAF
  dE <- r2$ELEAF - r1$ELEAF
return(dA/dE)
}

Cis <- seq(100,379,length=101)
z <- sapply(Cis, getdAdE, VPD=1.5, Ca=380)

lam <- 0.002
f <- FARAO(lambda=lam, Ca=380, VPD=1.5)

# This shows that 
# a) dA2/d2E < 0 (as required by Cowan-Farquhar)
# b) FARAO works OK
plot(Cis,z,ylim=c(0,10), lwd=2, type='l',
     xlab="Ci",ylab=expression(dA/dE))
points(f$Ci, 1000*lam, pch=19, col="red")


# now with energy balance
getdAdE2 <- function(Ci,...,returnwhat=c("dAdE","both")){
  returnwhat <- match.arg(returnwhat)
  delta <- 1e-03
  r1 <- PhotosynEB(Ci=Ci, ...)
  r2 <- PhotosynEB(Ci=Ci+delta, ...)
  dA <- r2$ALEAF - r1$ALEAF
  dE <- r2$ELEAF - r1$ELEAF
  if(returnwhat == "dAdE")
    return(dA/dE)
  else
    return(c(dA=dA, dE=dE))
}

run0 <- sapply(Cis, getdAdE2, Wind=50)
run1 <- sapply(Cis, getdAdE2, Wind=1)
run2 <- sapply(Cis, getdAdE2, Wind=0.2)

plot(Cis, run0, type='l', ylim=c(-1,10), xlab="Ci", ylab="dA/dE")
points(Cis, run1, type='l', col="blue")
points(Cis, run2, type='l', col="red")


r <- as.data.frame(do.call(rbind, lapply(Cis, getdAdE2, Wind=0.2, returnwhat="both")))

p <- PhotosynEB(Ci=seq(250,350,length=101))


# finding discontinuities

f <- plantecophys:::OPTfunEB
fv <- Vectorize(f)

curve(fv(Ci=x), from=50, to=380)
curve(fv(Ci=x,Wind=0.01), from=50, to=380, col="red")
curve(fv(Ci=x,Wind=0.011), from=50, to=380, col="red")
curve(fv(Ci=x,Wind=0.1), from=50, to=380, col="red")
curve(fv(Ci=x,Wind=0.5), from=50, to=380, col="blue", add=T)


Cis <- seq(100, 370, length=101)
run1 <- as.data.frame(do.call(rbind,lapply(Cis, 
                      function(x)f(x, Wind=0.01, retobjfun=FALSE))))
run2 <- as.data.frame(do.call(rbind,lapply(Cis, 
                      function(x)f(x, Wind=0.5, retobjfun=FALSE))))


with(run1, plot(GS, Tleaf))

gv <- Vectorize(plantecophys:::FindTleaf)
curve(gv(x,25), from=0.01, to=0.8)




FFindTleaf <- function(Ci=100, Tair=25, Wind=2, VPD=1.5, Wleaf=0.02, 
                       StomatalRatio=1,
                       LeafAbs=0.86,Ca=380,...){

  
  gsfun <- function(Ci, VPD, returnwhat=c("gs","all"), ...){
    
    returnwhat <- match.arg(returnwhat)
  
    run <- Aci(Ci, VPD=VPD, ...)   # note that VPD does not do anything, just for consistency in I/O
    
    A <- run$ALEAF
    
    # Given Ci and A, calculate gs (diffusion constraint)
    gs <- 1.57 * A / (Ca - Ci)
    
    if(returnwhat == "gs")return(gs)
    if(returnwhat == "all")return(list(run=run, GS=gs, A=A))
  }
  
  fx <- function(x, Ci, Tair, Wind, VPD, Wleaf, StomatalRatio, LeafAbs, ...){
    newx <- plantecophys:::FindTleaf(Tair=Tair, gs=gsfun(Ci=Ci, Tleaf=x, VPD=VPD, ...), 
                      Wind=Wind, Wleaf=Wleaf, 
                      StomatalRatio=StomatalRatio, LeafAbs=LeafAbs)
    newx - x
  }

  fx2 <- function(x, Ci, Tair=25, Wind=0.01, VPD=1.5, Wleaf=0.02, StomatalRatio=1, LeafAbs=0.86, ...){
    newx <- plantecophys:::FindTleaf(Tair=Tair, gs=gsfun(Ci=Ci, Tleaf=x, VPD=VPD, ...), 
                                     Wind=Wind, Wleaf=Wleaf, 
                                     StomatalRatio=StomatalRatio, LeafAbs=LeafAbs)
    newx - x
  }
  browser()
  
  Tleaf <- uniroot(fx, interval=c(Tair-15, Tair+15), Ci=Ci, Tair=Tair, Wind=Wind, VPD=VPD, Wleaf=Wleaf, 
                   StomatalRatio=StomatalRatio, LeafAbs=LeafAbs, ...)$root

return(Tleaf)
}


yv <- Vectorize(FFindTleaf)
curve(yv(x, Wind=0.1), from=200, to=379, type='o')

curve(yv(x, Wind=0.1), from=292, to=295, type='o')




