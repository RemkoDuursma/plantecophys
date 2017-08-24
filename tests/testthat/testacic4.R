library(plantecophys)
context("ACi C4")


a1 <- AciC4(Ci=seq(100,1000,length=25))
a2 <- AciC4(Ci=seq(100,1000,length=25), Tleaf=-5)


test_that("Aci C4 output", {
  expect_named(a1)
  expect_lt(max(a1$ALEAF) - min(a1$ALEAF), 5)
  expect_lt(max(a2$ALEAF), max(a1$ALEAF))
})
