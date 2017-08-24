library(plantecophys)
context("Tuzet")


t1 <- PhotosynTuzet()


test_that("Tuzet model output", {
  expect_named(t1)
  expect_lt(t1$PSIL, 0)
  expect_gt(t1$GS, 0)
})
