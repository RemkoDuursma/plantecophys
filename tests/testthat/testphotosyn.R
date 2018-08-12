
library(plantecophys)
context("Photosyn test")

p0 <- Photosyn()  # default settings
p0_2 <- Photosyn(Tcorrect=FALSE)  # should do nothing

vpds <- seq(0,5,length=20)
p1 <- Photosyn(VPD=vpds)

ppfds <- seq(0,2000,length=20)
p2 <- Photosyn(PPFD=ppfds)

p3_1 <- Photosyn(GammaStar=50)
p3_2 <- Photosyn(GammaStar=60)

p4_1 <- Photosyn(GS=0.3)
p4_2 <- Photosyn(GS=0.2)

p5_1 <- Photosyn(gmeso=0.4)
p5_2 <- Photosyn(gmeso=0.2)

p6_1 <- Photosyn(gsmodel="BBLeuning", TPU = NULL)
p6_2 <- Photosyn(gsmodel="BallBerry", whichA = "Aj")
p6_3 <- Photosyn(gsmodel="BallBerry", RH=VPDtoRH(1.5,25), VPD = NULL, whichA = "Ac")
p6_4 <- Photosyn(gsmodel="BallBerry", RH=VPDtoRH(1.5,25)-2)
p6_5 <- Photosyn(GS=-0.01)

p7 <- Photosyn(Ci=200, Tleaf=c(20,25))  # Used to be 'Jena bug'

p8_1 <- Photosyn(Ci=800, Ca=1200, TPU=8)
p8_2 <- Photosyn(Ci=800, Ca=1200, TPU=5)

p9 <- Photosyn(EdVC = 1)  # for a subclause test

test_that("Photosyn outputs", {
  expect_lte(max(diff(p1$GS)), 0)  # GS should always decrease with VPD
  expect_equal(0, p2$GS[1])
  expect_equal(p0$ALEAF, p0_2$ALEAF)
  expect_equal(length(ppfds), nrow(p2))
  expect_lt(p3_2$ALEAF, p3_1$ALEAF)  # increase GammaStar, decreases ALEAF
  expect_lt(p4_2$ALEAF, p4_1$ALEAF)  # decrease GS, decrease ALEAF
  expect_lt(p5_1$ALEAF, p0$ALEAF)    # include gmeso, decrease ALEAF
  expect_lt(p5_2$ALEAF, p5_1$ALEAF)  # decrease gmeso, decrease ALEAF
  expect_lt(p6_1$GS, p0$GS)          # default Leuning will have lower gs
  expect_lt(p6_2$GS, p0$GS)          # default BallBerry will have much lower gs
  expect_equal(p6_3$GS, p6_2$GS)     # RH calculated from default VPD, should give same answer.
  expect_lt(p6_4$GS, p6_3$GS)        # decrease RH, lower GS
  expect_equal(nrow(p7), 2)
  expect_lt(p8_2$ALEAF, p8_1$ALEAF)
})

test_that("Photosyn names", {
  expect_named(Photosyn(returnParsOnly=TRUE), c("Vcmax","Jmax","Km","GammaStar","VJ"))
  expect_named(p0)
})

test_that("Photosyn exceptions", {
  expect_error(Photosyn(gsmodel="BBmult"))
  expect_error(Photosyn(VPD=NULL, RH=NULL))
  expect_error(Photosyn(GS=0.2, Ci=200))
  expect_error(Photosyn(gsmodel = "BBdefine"))
})




