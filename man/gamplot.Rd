\name{gamplot}
\alias{gamplot}
\title{Plot a simple generalized additive model fit}
\usage{
  gamplot(xvarname, yvarname, groupname = "", dfr,
    add = FALSE, CI = TRUE, CIlevel = 0.95,
    polycolor = "grey", line = TRUE, linecolor = "black",
    lwd = 2, optline = FALSE, ...)
}
\arguments{
  \item{xvarname}{Name of the X variable (quoted)}

  \item{yvarname}{Name of the Y variable (quoted)}

  \item{groupname}{Name of the grouping variable (quoted,
  optional)}

  \item{dfr}{Dataframe that contains the above three
  variables}

  \item{add}{Logical, if TRUE adds to current plot}

  \item{CI}{Logical, if TRUE (the default), adds a shaded
  confidence region}

  \item{CIlevel}{The coverage of the confidence region
  (defaults to 0.95).}

  \item{polycolor}{Colour of the polygon}

  \item{line}{Logical, if TRUE adds a line for the mean
  prediction}

  \item{linecolor}{Color of the line (can be a vector with
  colors for each group)}

  \item{lwd}{Width of the line}

  \item{optline}{Logical, if TRUE adds a vertical line at
  the maximum value of Y.}
}
\description{
  A wrapper for \code{gam} from the \code{mgcv} package, to
  produce a nice smoother with a 95% confidence region.
}

