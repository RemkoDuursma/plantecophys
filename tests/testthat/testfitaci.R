
library(plantecophys)
context("Fit ACi curves")

# Load example dataset
acidat <- droplevels(manyacidat[manyacidat$Curve %in% levels(manyacidat$Curve)[1:10],])
acidatone <- droplevels(subset(acidat, Curve == "1000_2_3"))


# Fit one A-Ci curve
fit1 <- fitaci(acidatone,  fitmethod="default")
fit2 <- fitaci(acidatone,  fitmethod="bilinear")

# Tcorrect, or not.
fit3.1 <- fitaci(acidatone,  fitmethod="default", Tcorrect=TRUE, Tleaf=15)
fit3.2 <- fitaci(acidatone,  fitmethod="default", Tcorrect=FALSE, Tleaf=15)
fit3.3 <- fitaci(acidatone,  fitmethod="bilinear", Tcorrect=TRUE, Tleaf=15)
fit3.4 <- fitaci(acidatone,  fitmethod="bilinear", Tcorrect=FALSE, Tleaf=15)
fit3.5 <- fitaci(acidatone,  fitmethod="bilinear", Tcorrect=FALSE, fitTPU=TRUE)

# With mesophyll conductance
fit5 <- fitaci(acidatone,  fitmethod="default", gmeso=0.9)
fit6 <- fitaci(acidatone,  fitmethod="default", gmeso=0.4)

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
fits2_pb <- capture.output(fitacis(acidat, "Curve",  fitmethod="bilinear", progressbar=TRUE))
fits3 <- fitacis(acidat, "Curve",  fitTPU=TRUE, progressbar=FALSE)
fits4 <- fitacis(acidat, "Curve",  fitmethod="bilinear", progressbar=FALSE, id="treatment")

# one-point
fit10.1 <- fitaci(acidatone, fitmethod = "onepoint")
fit10.2 <- fitaci(acidatone, fitmethod = "onepoint", Tcorrect=FALSE)

# Make some not fit
acidat$Photo[c(2,50,100)] <- 50
fits1_re <- fitacis(acidat, "Curve",  fitmethod="default", progressbar=FALSE, quiet=FALSE)

# Pass GammaStar, Km
fit11 <- fitaci(acidatone, fitmethod="bilinear", GammaStar=60, Km=800)

# Fit without PPFD, Tleaf
# Set Rd but don't actually use
acidatone$Tleaf <- acidatone$PARi <- NULL
acidatone$Rd <- 2


# Test print methods. 
print(fits4)
print(fits3)
print(fits2)
print(fits1)
print(fit5)

print(fit1)
print(fit2)
print(fit3.1)
print(fit3.2)
print(fit8.2)
print(fit7.1)
print(fit11)

# Same as print.
summary(fit1)

# 'Test' plot methods.
plot(fit1)
plot(fit1, xlim=c(100,200), ylim=c(0,30), linecols="blue")
plot(fit3.5)

plot(fits1)
plot(fits1, how="oneplot")
plot(fits2, highlight="1000_1_5", how="oneplot")
plot(fits3)
plot(fits4, lwd=2)
plot(fits4, how="oneplot", colour_by_id = TRUE, id_legend=TRUE)
plot(fits4, how="oneplot", colour_by_id = FALSE, id_legend=TRUE)


# 
test_that("Aci curve fit output format", {
  expect_error(fitaci(data.frame()))
  expect_length(fitted(fit1), nrow(acidatone))
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
  expect_warning(fitacis(acidat, "Curve",  fitmethod="bilinear", progressbar=FALSE, id="XXXX"))
  expect_true("treatment" %in% names(coef(fits4)))
  expect_warning(fit12 <- fitaci(acidatone, fitmethod="bilinear"))
  expect_warning(fitaci(acidatone, citransition=400, fitmethod="default"))
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
  expect_gt(fit1$Ci(ALEAF = 10), 150)
  expect_equal(fit1$Ci(ALEAF = 10^6), NA)
})


test_that("Aci curve RMSE", {
  expect_lt( max(vapply(fits1, "[[", numeric(1), "RMSE")), 5)
  expect_lt( max(vapply(fits2, "[[", numeric(1), "RMSE")), 5)
  expect_lt( max(vapply(fits3, "[[", numeric(1), "RMSE")), 5)
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


# Create empty level
acidat_err <- acidat
levels(acidat_err$Curve) <- letters[1:11]

palette(c("blue","red"))

fits4_bad <- fits4
fits4_bad[[1]]$pars[] <- NA

# bad Rd
acidat_err2 <- acidat_err
acidat_err2$Rd <- 2
acidat_err2$Rd[1] <- 2.1


test_that("fitacis exceptions", {
  expect_error(fitacis(acidat, "XXXX"))
  expect_error(fitacis(acidat_err, "Curve"))
  expect_error(plot(fits3, how="oneplot", colour_by_id = TRUE, id_legend=TRUE))
  expect_warning(plot(fits4, how="oneplot", colour_by_id = TRUE, id_legend=TRUE))
  expect_warning(plot(fits4, how="manyplots", add=TRUE))
  expect_error(plot(fits4, highlight="XXX", how="oneplot"))
  expect_true(is.na(coef(fits4_bad)[1,3]))
  expect_error(fitaci(acidatone, fitTPU=TRUE, fitmethod = "onepoint"))
  expect_warning(fitaci(acidatone, fitTPU=FALSE, fitmethod = "onepoint"))
  expect_warning(fitaci(acidatone,  fitmethod="bilinear"))
  expect_warning(fitaci(acidatone,  fitmethod="bilinear", PPFD = 1800))
  expect_error(fitaci(acidat_err2,  fitmethod="bilinear", useRd = TRUE))
  
  
})




