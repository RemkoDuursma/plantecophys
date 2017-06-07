
library(plantecophys)
context("Fit ACi curves")

# Load example dataset
acidat <- droplevels(manyacidat[manyacidat$Curve %in% levels(manyacidat$Curve)[1:10],])

# Fit one A-Ci curve
acidatone <- droplevels(subset(acidat, Curve == "1000_2_3"))
fit1 <- fitaci(acidatone,  fitmethod="default")
fit2 <- fitaci(acidatone,  fitmethod="bilinear")

# Tcorrect, or not.
fit3.1 <- fitaci(acidatone,  fitmethod="default", Tcorrect=TRUE, Tleaf=15)
fit3.2 <- fitaci(acidatone,  fitmethod="default", Tcorrect=FALSE, Tleaf=15)
fit3.3 <- fitaci(acidatone,  fitmethod="bilinear", Tcorrect=TRUE, Tleaf=15)
fit3.4 <- fitaci(acidatone,  fitmethod="bilinear", Tcorrect=FALSE, Tleaf=15)

# With mesophyll conductance
fit5 <- fitaci(acidatone,  fitmethod="default", gmeso=0.9)
fit6 <- fitaci(acidatone,  fitmethod="default", gmeso=0.3)

# Specify transition point
fit7.1 <- fitaci(acidatone,  fitmethod="default", citransition=400)
fit7.2 <- fitaci(acidatone,  fitmethod="bilinear", citransition=400)
fit7.3 <- fitaci(acidatone,  fitmethod="bilinear", citransition=fit2$Ci_transition)

# Use measured Rd
acidatone$Rd <- 2
fit8.1 <- fitaci(acidatone, useRd=TRUE, fitmethod="default")
fit8.2 <- fitaci(acidatone, useRd=TRUE, fitmethod="bilinear")
fit8.3 <- fitaci(acidatone, useRd=TRUE, fitmethod="default", Tcorrect=FALSE)
fit8.4 <- fitaci(acidatone, useRd=TRUE, fitmethod="bilinear", Tcorrect=FALSE)

# Fit many A-Ci curves
fits1 <- fitacis(acidat, "Curve",  fitmethod="default", progressbar=FALSE)
fits2 <- fitacis(acidat, "Curve",  fitmethod="bilinear", progressbar=FALSE)
fits3 <- fitacis(acidat, "Curve",  fitTPU=TRUE, progressbar=FALSE)

test_that("Aci curve fit output format", {
  expect_equal(names(coef(fit1)), c("Vcmax","Jmax","Rd"))
  expect_equal(names(coef(fit2)), c("Vcmax","Jmax","Rd"))
  expect_equal(names(coef(fits1)), c("Curve","Vcmax","Jmax","Rd","Vcmax_SE","Jmax_SE","Rd_SE"))
  expect_equal(names(coef(fits2)), c("Curve","Vcmax","Jmax","Rd","Vcmax_SE","Jmax_SE","Rd_SE"))
  expect_equal(names(coef(fits3)), c("Curve","Vcmax","Jmax","Rd","TPU","Vcmax_SE","Jmax_SE","Rd_SE","TPU_SE"))
  expect_equal(nrow(coef(fits1)), nlevels(acidat$Curve))
  expect_equal(nrow(coef(fits2)), nlevels(acidat$Curve))
  expect_s3_class(fit1, "acifit")
  expect_s3_class(fits1, "acifits")
  expect_gt(ncol(coef(fits3)), ncol(coef(fits2)))
})

test_that("Aci curve fitted coefficients",{
  expect_lt(coef(fit1)[1], coef(fit1)[2], "Vcmax","Jmax")
  expect_lt(coef(fit2)[1], coef(fit2)[2], "Vcmax","Jmax")
  expect_gt(coef(fit3.1)[1], coef(fit3.2)[1], "Vcmax Tcorrect", "Vcmax no Tcorrect")
  expect_gt(coef(fit3.3)[1], coef(fit3.4)[1], "Vcmax Tcorrect", "Vcmax no Tcorrect")
  expect_gt(min(coef(fit1)), 0)
  expect_gt(min(coef(fit2)), 0)
  expect_gt(coef(fit6)[1], coef(fit5)[1], "Vcmax gmeso = 0.3", "Vcmax gmeso = 0.9")
  expect_equal(fit7.3$Ci_transition, fit7.2$Ci_transition)
})

test_that("Aci curve RMSE", {
  expect_lt( max(sapply(fits1, "[[", "RMSE")), 5)
  expect_lt( max(sapply(fits2, "[[", "RMSE")), 5)
  expect_lt( max(sapply(fits3, "[[", "RMSE")), 5)
  expect_lt(fit1$RMSE, 5)
  expect_lt(fit2$RMSE, 5)
  expect_lt(fit3.1$RMSE, 5)
  expect_lt(fit3.2$RMSE, 5)
  expect_lt(fit3.3$RMSE, 5)
  expect_lt(fit3.4$RMSE, 5)
  expect_lt(fit5$RMSE, 5)
  expect_lt(fit6$RMSE, 7)
  expect_lt(fit5$RMSE, 5)
  expect_lt(fit7.1$RMSE, 5)
  expect_lt(fit7.2$RMSE, 5)
  expect_lt(fit7.3$RMSE, 5)
  expect_lt(fit8.1$RMSE, 5)
  expect_lt(fit8.2$RMSE, 5)
  expect_lt(fit8.3$RMSE, 5)
  expect_lt(fit8.4$RMSE, 5)
})
