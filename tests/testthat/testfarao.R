library(plantecophys)
context("FARAO")



f1 <- FARAO()
f2 <- FARAO2()
f3 <- FARAO2(energybalance=TRUE)
f4 <- FARAO(C4=TRUE)

test_that("FARAO output", {
  expect_named(f1)
  expect_named(f2)
  expect_named(f3)
  expect_named(f4)
})
