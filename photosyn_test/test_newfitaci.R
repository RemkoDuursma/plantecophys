
library(plantecophys)
acidata1$Tleaf <- 33

x <- fitaci(acidata1, PPFD=1000, fitmethod="transit", Tcorrect=T)
coef(x)
y <- fitaci(acidata1, PPFD=1000, fitmethod="de", Tcorrect=T)
coef(y)


plot(x)
plot(y)



f <- fitaci(acidata1)


get_rmse <- function(x)sqrt(mean((x$df$Ameas - x$df$Amodel)^2))
cis <- seq(200, 700, length=101)
fits <- lapply(cis, function(x)get_rmse(fitaci(acidata1, fitmethod="transit",
                                               PPFD=1800, citransition=x)))
g <- fitaci(acidata1, fitmethod="transit", PPFD=1800)

windows()
par(mfrow=c(2,1), mar=c(3,3,0,0))
plot(cis, fits, xlim=c(0,1200))
plot(g, xlim=c(0,1200))



# Court
dat <- read.csv("photosyn_test/ACi_court2.csv")
dat$CurveID <- with(dat, as.factor(paste(Volume, plot, pot)))
vn <- list(ALEAF="Ps", Tleaf="Temp", Ci="Ci", PPFD="PPFD")
f1_1 <- fitacis(dat, "CurveID", fitmethod="transit", varnames=vn)
f1_2 <- fitacis(dat, "CurveID", fitmethod="def", varnames=vn)

plot(coef(f1_1)$Vcmax, coef(f1_2)$Vcmax)
abline(0,1)

plot(coef(f1_1)$Vcmax_SE, coef(f1_2)$Vcmax_SE)
abline(0,1)

plot(coef(f1_1)$Jmax, coef(f1_2)$Jmax)
abline(0,1)

plot(sapply(f1_1,get_rmse), sapply(f1_2, get_rmse))
abline(0,1)


#
glash <- read.csv("photosyn_test/Aci_data_glasshouse.csv")

f2_1 <- fitacis(glash, "Plant", fitmethod="transit")
f2_2 <- fitacis(glash, "Plant", fitmethod="def")

cf1 <- coef(f2_1)
cf1$RMSE <- sapply(f2_1, get_rmse)

cf2 <- coef(f2_2)
cf2$RMSE <- sapply(f2_2, get_rmse)

p <- merge(cf1, cf2, by="Plant", all=TRUE)

plot(p$Vcmax.x, p$Vcmax.y)
abline(0,1)

plot(p$Vcmax_SE.x, p$Vcmax_SE.y)
abline(0,1)

plot(p$Jmax.x, p$Jmax.y)
abline(0,1)

plot(p$RMSE.x, p$RMSE.y)
abline(0,1)

# problem: because of discontinuity?
plot(f2_1[[6]])
plot(f2_2[[6]])




