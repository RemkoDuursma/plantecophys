

######################################################################################
######################################################################################
# This section will install the packages require to fit the A:Ci curves.

#install devtools, which is require to install from bitbucket. You only need to run this line ONCE.
install.packages("devtools")

#load the devtools library, so it can be used. You only need to run this line once per R session.
library(devtools)

#install Remko's package, includes the functions fitaci() and fitacis(). Just run ONCE.
install_bitbucket(repo="plantecophys",username="remkoduursma")

#load Remko's package. Run once per R session.
library(plantecophys)
######################################################################################
######################################################################################



######################################################################################
######################################################################################
#read in the data and fit the curves
######################################################################################

#read in the example dataset of 81 A:Ci curves. 
#You will need to edit this code to match where the data are stored on your computer.
dat <- read.csv("./data/Aci_data_glasshouse.csv")

#look at the top bit of the data. The object "dat" is a dataframe. The first column
# is a factor that identifies each curve, and the subsequent columns contain the
# data required to fit A:Ci curves, plus a bit more.
head(dat)

#look at the structure of the data. "Plant" is a factor, the rest of the columns are numeric
# The function will look for columns titled "Photo", "Ci", "Tleaf", and "PPFD". If the
# function can't find these columns, it will use default parameters that may or may not be appropriate for your data.
str(dat)

######################################################################################
#extract just one curve to fit individually, as an example
tKR9 <- subset(dat,Plant=="tKR9")

#fit this one curve
tKR9.fit <- fitaci(tKR9)
#plot the curve
plot(tKR9.fit)
#extract the fitted parameters (scaled to 25 deg C)
tKR9.fit$pars

######################################################################################
#fit all the curves. Pass the name of the column that identifies each curve to the function in quotes.
# In this case, the column "Plant" contains the unique identifier for each curve.
fits <- fitacis(dat,"Plant")

#fits is a complex object- it is a list of 81 aci-curve fits. 
fitted.pars <- coef(fits)

#Done! 81 A:Ci curves fit simply, quickly, and reproducibly.

#now you can analyze the results however you'd like. The output looks believeable, as
#        the ratio of Jmax/Vcmax is similar to the expected ratio of ~1.7
with(fitted.pars,plot(Jmax~Vcmax))
abline(0,1,lty=2)
my.lm <-lm(Jmax~Vcmax,data=fitted.pars)
abline(my.lm)
legend("topleft",legend=paste("slope = ",round(my.lm$coefficients[2],2),", r2 = ",round(summary(my.lm)$r.squared,2)))


######################################################################################
######################################################################################