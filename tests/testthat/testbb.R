
library(plantecophys)
context("Fit BB")

set.seed(1)
n <- 500
dfr <- data.frame(PPFD=runif(n, 100, 1000), 
                  VPD=runif(n,1,4),
                  Tleaf=runif(n,20,28),
                  ID=rep(letters[1:(n/50)], each=50))
p <- Photosyn(PPFD=dfr$PPFD, VPD=dfr$VPD, Tleaf=dfr$Tleaf)
dfr$gs <- p$GS + rnorm(n, 0, 0.03)
dfr$aleaf <- p$ALEAF
dfr$Ca <- 400
dfr$RH <- VPDtoRH(dfr$VPD, dfr$Tleaf)/100
dfr$RHperc <- VPDtoRH(dfr$VPD, dfr$Tleaf)

fit1 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", VPD="VPD", Ca="Ca"))
fit2 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", VPD="VPD", Ca="Ca"), fitg0=TRUE)
fit3 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", RH="RH", Ca="Ca", VPD="VPD"), 
              gsmodel="BallBerry",
              fitg0=TRUE)

fit3.2 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", RH="RHperc", Ca="Ca", VPD="VPD"), 
              gsmodel="BallBerry",
              fitg0=TRUE)

fit4 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", Ca="Ca", VPD="VPD"), 
              gsmodel="BBLeuning",
              fitg0=TRUE)
fit5 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", Ca="Ca", VPD="VPD"), 
              gsmodel="BBOptiFull",
              fitg0=TRUE)

fit2_2 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", VPD="VPD", Ca="Ca"), fitg0=FALSE)
fit3_2 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", RH="RH", Ca="Ca", VPD="VPD"), 
              gsmodel="BallBerry",
              fitg0=FALSE)

fit3.2_2 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", RH="RHperc", Ca="Ca", VPD="VPD"), 
                gsmodel="BallBerry",
                fitg0=FALSE)

fit4_2 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", Ca="Ca", VPD="VPD"), 
              gsmodel="BBLeuning",
              fitg0=FALSE)
fit5_2 <- fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", Ca="Ca", VPD="VPD"), 
              gsmodel="BBOptiFull",
              fitg0=FALSE)

fits1 <- fitBBs(dfr, "ID",varnames=list(ALEAF="aleaf", GS="gs", Ca="Ca", VPD="VPD"))
fits2 <- fitBBs(dfr, "ID",varnames=list(ALEAF="aleaf", GS="gs", Ca="Ca", VPD="VPD"), fitg0=TRUE)
print(fits1)
print(fits2)
print(coef(fits1))

print(fit2)
print(fit2_2)
print(fit5_2)
print(fits1)

test_that("Fit BB output", {
  expect_named(fit1)
  expect_named(fit2)
  expect_length(coef(fit1),2)
  expect_length(coef(fit2),2)
  expect_equal(coef(fit1)[[1]], 0)
  expect_gt(coef(fit2)[[1]],0)
  expect_equal(coef(fit3), coef(fit3.2))
})


# Missing data
var_names <- list(ALEAF="aleaf", GS="gs", VPD="VPD", Ca="Ca")
dfr_gs <- dfr[, -match("gs", names(dfr))]
dfr_vpd <- dfr[, -match("VPD", names(dfr))]
dfr_aleaf <- dfr[, -match("aleaf", names(dfr))]
dfr_ca <- dfr[, -match("Ca", names(dfr))]
dfr_rh <- dfr[, -match("RH", names(dfr))]

test_that("fitBB exceptions", {
  expect_error(fitBB(dfr, varnames=list(ALEAF="xx", GS="gs", VPD="VPD", Ca="Ca")))
  expect_error(fitBB(dfr, varnames=list(ALEAF="aleaf", GS="gs", VPD="xxx", Ca="Ca")))
  
  expect_error(fitBB(dfr_gs, varnames=var_names))
  expect_error(fitBB(dfr_vpd, varnames=var_names))
  expect_error(fitBB(dfr_aleaf, varnames=var_names))
  expect_error(fitBB(dfr_ca, varnames=var_names))
  expect_error(fitBB(dfr_rh, varnames=var_names, gsmodel="BallBerry"))
})





