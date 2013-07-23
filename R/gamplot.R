#' Plot a simple generalized additive model fit
#'@description A wrapper for \code{gam} from the \code{mgcv} package, to produce a  smoother with a confidence region.
#'@param xvarname Name of the X variable (quoted)
#'@param yvarname Name of the Y variable (quoted)
#'@param groupname Name of the grouping variable (quoted, optional)
#'@param dfr Dataframe that contains the above three variables
#'@param add Logical, if TRUE adds to current plot
#'@param CI Logical, if TRUE (the default), adds a shaded confidence region
#'@param CIlevel The coverage of the confidence region (defaults to 0.95).
#'@param polycolor Colour of the polygon
#'@param line Logical, if TRUE adds a line for the mean prediction
#'@param linecolor Color of the line (can be a vector with colors for each group)
#'@param lwd Width of the line
#'@param optline Logical, if TRUE adds a vertical line at the maximum value of Y.
#'@param gamfit Optional, an object returned by \code{gam}, which can yield predictions for the dfr provided as input. 
#'
#'@details If a fitted gam object is not provided, this function fits a vanilla smoother function when there is no grouping variable,
#'\preformatted{
#'gamfit <- gam(Y ~ s(X), data=dat)
#'}
#'
#'When a grouping variable is provided, it fits :
#'\preformatted{
#'gamfit <- gam(Y ~ group + s(X, by=group, id=1), data=dat)
#'}
#'
#'See the documentation for \code{\link{mgcv::gam}} in the \code{mgcv} package.
#' @return The fitted gam object, invisibly.
#' @export
gamplot <- function(xvarname, yvarname, groupname="", dfr, add=FALSE,
                    CI=TRUE, CIlevel=0.95, polycolor="grey", line=TRUE, 
                    linecolor="black", lwd=2,
                    optline=FALSE, gamfit=NULL, ...){
  
  
  dfr$Y <- dfr[,yvarname]
  dfr$X <- as.numeric(dfr[,xvarname])
  
  if(groupname != "")
    dfr$group <- as.factor(dfr[,groupname])
  else
    dfr$group <- as.factor("a")
  
  dat <- dfr[complete.cases(dfr),]
  
  if(is.null(gamfit)){
    if(nlevels(dat$group) > 1)
      gamfit <- gam(Y ~ group + s(X, by=group, id=1), data=dat)
    else
      gamfit <- gam(Y ~ s(X), data=dat)
  }
  
  p <- predict(gamfit, dat, se.fit=TRUE)
  
  dat$Ypred <- p[[1]]
  dat$YpredSE <- p[[2]]
  
  if(!add){
    with(dfr, plot(X, Y, type='n', ...))
    box()
  }
  
  
  addgampoly <- function(dataset,...){
    alpha <- 1-CIlevel
    qn <- qnorm(1-alpha/2)
    x <- dataset[order(dataset$X),]
    with(x,
         polygon(x=c(X, rev(X)), y=c(Ypred + qn*YpredSE, rev(Ypred - qn*YpredSE)), 
                 col=polycolor, border=NA))
  }
  
  # gam polygon (CI)
  if(CI)lapply(split(dat,dat$group),function(x)addgampoly(x))
  
  addgamline <- function(dataset,...){
    x <- dataset[order(dataset$X),]
    points(x$X, x$Ypred, type='l',  ...)
  }
  
  # line (mean)
  nl <- nlevels(dat$group)
  nc <- length(linecolor)
  if(nc < nl)linecolor <- rep(linecolor, ceiling(nl/nc))
  
  for(i in 1:nl){
    addgamline(dat[dat$group==levels(dat$group)[i],],
               lwd=lwd, col=linecolor[i],...)
    
  }
  
  addOptline <- function(x,...){
    xmax <- x$X[which.max(x$Ypred)]
    ymax <- max(x$Ypred, na.rm=TRUE)
    segments(x0=xmax, x1=xmax, y0=0, y1=ymax, ... ) 
  }
  
  if(optline){
    for(i in 1:nl){
      Ymaxline(dat[dat$group==levels(dat$group)[i],],
               lwd=lwd, col=linecolor[i],...)
    }
  }
  
return(invisible(gamfit))
}