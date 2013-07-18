#' @S3method plot acifit
#' @method plot acifit
#' @rdname fitaci
plot.acifit <- function(x, ...){
  
  
  Ci <- with(x$df, seq(min(Ci), max(Ci), length=101))
  
  # Exact model used to fit the A-Ci curve was saved in the object.
  pred <- x$Photosyn(Ci=Ci, Vcmax=x$pars[1], Jmax=x$pars[2], Rd=x$pars[3],
                   Tleaf=mean(x$df$Tleaf), PPFD=mean(x$df$PPFD))
  
  with(x$df, plot(Ci, Ameas, pch=19,
               ylim=c(min(Ameas), 1.1*max(Ameas)),
               xlim=c(0, max(Ci)),
               xlab=expression(italic(C)[i]~~(ppm)),
               ylab=expression(italic(A)[net]~~(mu*mol~m^-2~s^-1)),
               ...
              ))
  
  with(pred, points(Ci, Aj-Rd, type='l', col="blue"))
  with(pred, points(Ci, Ac-Rd, type='l', col="red"))
  with(pred, points(Ci, ALEAF, type='l', col="black", lwd=2))
  abline(h=0, lty=3)
  legend("bottomright", c(expression(italic(A)[c]),
                      expression(italic(A)[j]),
                                 "Limiting rate"), lty=1, lwd=c(1,1,2), col=c("red","blue","black"))
  
}