The plantecophys R package
==========================


Welcome to the home of *plantecophys*, an R package that bundles a number of tools to analyze and visualize plant ecophysiological data, particularly leaf gas exchange data.

The following functions are the core of the package:

**Leaf gas exchange**

* `Photosyn` is an implementation of a coupled leaf gas exchange model, including the Farquhar-Berry-von Caemmerer model of photosynthesis, and Ball-Berry-type stomatal conductance models.
* `Aci` calculates the dependence of photosynthesis (A) on the intercellular CO2 concentration (Ci), with the Farquhar model. 
* `fitaci` and `fitacis` fit the Farquhar model of photosynthesis to observations of A and Ci as measured with a Licor 6400 or similar instrument.
* `FARAO` is an implementation of the stomatal optimization model as suggested by Cowan and Farquhar (1977). It is a full numerical solution to the model, and uses the FARquhar model for photosynthesis.

**Unit conversion**

* `RHtoVPD` converts from relative humidity to VPD. There are others as well (see `?RHtoVPD`).


**visualization**

* `gamplot` is a simple wrapper to produce smoothed regressions with confidence intervals, using generalized additive model fits. This is a useful function to visualize many types of quantitative relationships for which no model can be easily specified.


### Installation

The package is not hosted on CRAN yet, so use the following command to install it. Windows users must have Rtools installed (http://cran.r-project.org/bin/windows/Rtools/)

```
library(devtools)
install_bitbucket("plantecophys","remkoduursma")
```


For questions, comments, and suggestions, please email remkoduursma@gmail.com