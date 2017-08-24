library(plantecophys)
context("Leaf energy balance")


p0 <- PhotosynEB()

p0_1 <- PhotosynEB(RH=VPDtoRH(1.5, 2.5))

p1 <- PhotosynEB(Wind=0.2)
p2 <- PhotosynEB(VPD=4)


test_that("Leaf energy balance", {
  expect_gt(p1$Tleaf, p0$Tleaf)
  expect_gt(p2$Tleaf, p0$Tleaf)
})
