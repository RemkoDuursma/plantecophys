library(plantecophys)
context("FARAO")



f1 <- FARAO()
f2 <- FARAO(energybalance=TRUE)
f3 <- FARAO2()
f4 <- FARAO2(energybalance=TRUE)
f5 <- FARAO(C4=TRUE)

test_that("FARAO output", {
  expect_named(f1)
  expect_named(f2)
  expect_named(f3)
  expect_named(f4)
})
