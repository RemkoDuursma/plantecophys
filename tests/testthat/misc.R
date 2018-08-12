
library(plantecophys)
context("Misc tests")


plantecophys

test_that("QUADP", {
 
  expect_equal(plantecophys:::QUADP(0,0,0), 0)
  expect_equal(plantecophys:::QUADM(0,0,0), 0)
  expect_equal(plantecophys:::QUADM(0,1,1), -1)

  expect_warning(plantecophys:::QUADP(1,2,3))
  expect_warning(plantecophys:::QUADM(1,2,3))
})
