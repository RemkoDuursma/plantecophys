
library(plantecophys)
f <- fitacis(manyacidat, "Curve", fitmethod="bilinear")

get_coefs <- function(x){
  
  data.frame(citransition=x$Ci_transition,
             A1000=x$Photosyn(Ci=1000)$ALEAF,
             A380=x$Photosyn(Ci=380)$ALEAF)
  
}
coefs <- do.call(rbind, lapply(f, get_coefs))

coefs <- cbind(coef(f), coefs)
