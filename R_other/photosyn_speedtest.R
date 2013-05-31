
setwd("G:/Work/Projects/WUE Project/plantecophys package")
library(reshape)

# Another speed test


  library(GasExchangeR)  # photosyn
  source("./R_other/photosyn2.R")  # photosyn2
  source("./R/photosyn.R")  # Photosyn


  makedf <- function(n){
    
    PPFD <- runif(n, 40, 1500)
    Tleaf <- runif(n, 5,40)
    Ca <- runif(n, 380,550)
    VPD <- runif(n, 1,3)
    
  data.frame(PPFD=PPFD, Tleaf=Tleaf, Ca=Ca, VPD=VPD)  
  }
  
  v <- list()
  ns <- 10^seq(3,4.5,length=15)
  for(i in seq_along(ns)){
  
    df <- makedf(ns[i])
    
    t1 <- system.time(photosyn(PAR=df$PPFD, VPD=df$VPD, TLEAF=df$Tleaf))
    t2 <- system.time(photosyn2(PPFD=df$PPFD, VPD=df$VPD, Tleaf=df$Tleaf))
    t3 <- system.time(Photosyn(PPFD=df$PPFD, VPD=df$VPD, Tleaf=df$Tleaf))
    v[[i]] <- c(t1[[3]], t2[[3]], t3[[3]])
  }
  
  v <- as.data.frame(do.call(rbind, v))
  names(v) <- c("t1","t2","t3")
  
  v$n <- ns
  vm <- melt(v, "n")
  names(vm) <- c("n","model","time")
  
  windows()
  palette(c("red","forestgreen","blue"))
  with(vm, plot(log10(n), log10(time), pch=19, col=model))
  legend("topleft", c("GasEchangeR::photosyn","photosyn2","plantecophys::Photosyn"), 
         fill=palette())

  ## Scaling of simulation time
  p <- coef(lm(log10(time) ~ model + model:log10(n)-1, data=vm))  
  matrix(p, ncol=2, dimnames=list(NULL, c("constant","exponent")))
  
  # So n=100,000 takes approx this many seconds:
  100000 * 10^p[[1]]
  100000 * 10^p[[2]]
  100000 * 10^p[[3]]
  
  
  
