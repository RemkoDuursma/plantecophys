library(plantecophys)
context("Unit conversions")


c1 <- VPDtoDew(2, 20)
c2 <- DewtoVPD(10,20)
c3 <- VPDleafToAir(2, 20, 22)
c4 <- VPDairToLeaf(2, 20, 22)
c5 <- RHleafToAir(30, 20, 22)
c6 <- RHairToLeaf(30, 20, 22)

# Should have some tests (esp. converting one way and back)