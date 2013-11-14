\name{predict_nls}
\alias{plot.nls}
\alias{predict_nls}
\title{Predict from a non-linear regression object}
\usage{
  predict_nls(object, from = NULL, to = NULL, x = NULL,
    interval = c("none", "confidence"), level = 0.95,
    n = 101, nboot = 999, add = TRUE, ...)

  plot.nls(x, add = FALSE, lwd = c(1, 1), lty = c(1, 5),
    col = c("black", "red"), xlim = NULL, ylim = NULL,
    xlab = NULL, ylab = NULL, ...)
}
\arguments{
  \item{object}{An object of class 'nls'}

  \item{from}{}
}
\description{
  Predicts from a fitted non-linear regression object, as
  returned by \code{\link{nls}}. Also calculates a
  confidence interval using a case bootstrap resampling
  method (using \code{\link{bootCase}}) from the \code{car}
  package.
}
\examples{
\dontrun{
f <- function(t, k=0.5)1 - exp(-k*t)
x <- seq(0,12,length=101)
y <- f(x) + rnorm(101,sd=0.1)
dfr <- data.frame(x=x,y=y)
nls1 <- nls(y ~ b0 + A*(1-exp(-k*x)), start=list(b0=0,A=1,k=1), data=dfr)

p <- predict_nls(nls1, from=min(x), to=max(x), interval="confidence")

plot(x,y)
with(p, {
  lines(x, pred, lty=5)
  lines(x, lwr, lty=4, col="red")
  lines(x, upr, lty=4, col="red")
})

# Compare to predictNLS, which uses a 2nd order Taylor expansion.
library(propagate)
P <- predictNLS(nls1, newdata = data.frame(x=p$x))
lines(p$x, P$summary[,5], col="blue", lty=5)
lines(p$x, P$summary[,6], col="blue", lty=5)
}

