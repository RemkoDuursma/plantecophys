plot.acifit <- function(x, ...){
  
  
  Ci <- seq(min(x$data$Ci), max(x$data$Ci), length=101)
  pred <- acifun(Ci=Ci, Vcmax=x$pars[[1]], Jmax=x$pars[[2]], Rd=x$pars[[3]],
                   Tleaf=mean(x$data$Tleaf), PAR=mean(x$data$PAR))
  
  with(x$data, plot(Ci, Photo, pch=19,
               ylim=c(min(Photo), 1.1*max(Photo)),
               xlim=c(0, max(Ci)),
               xlab=expression(italic(C)[i]~~(ppm)),
               ylab=expression(italic(A)[net]~~(mu*mol~m^-2~s^-1))                 
              ))
  
  with(pred, points(Ci, Aj, type='l', col="blue"))
  with(pred, points(Ci, Ac, type='l', col="red"))
  with(pred, points(Ci, Am, type='l', col="black", lwd=2))
  abline(h=0, lty=3)
  legend("topleft", c(expression(italic(A)[c]),
                      expression(italic(A)[j]),
                                 "Limiting rate"), lty=1, lwd=c(1,1,2), col=c("red","blue","black"))
  
}