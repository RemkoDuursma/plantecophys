
library(plantecophys)

dfr <- read.csv("photosyn_test/aciexample.csv")

c3 <- subset(dfr, Identity == 3)
f <- fitaci(c3)


fits <- fitacis(dfr, "Identity")
