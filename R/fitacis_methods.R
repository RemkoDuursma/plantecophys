#' @method plot acifits
#' @export
#' @rdname fitacis
#' @importFrom grDevices palette rainbow
plot.acifits <- function(x, how=c("manyplots","oneplot"),
                         highlight=NULL, 
                         ylim=NULL, 
                         xlim=NULL,
                         add=FALSE, 
                         what=c("model","data","none"),
                         colour_by_id = FALSE,
                         id_legend=TRUE,
                         linecol_highlight = "black",
                         ...){
  
  how <- match.arg(how)
  what <- match.arg(what)
  
  if(colour_by_id){

    if(is.null(x[[1]]$id))
      Stop("To colour curves by id, fit with id argument (see ?fitacis).")

    id_fac <- sapply(x, function(fit)unique(fit$df[,fit$id]))
    if(nlevels(id_fac) > length(palette())){
      pal <- rainbow(nlevels(id_fac))
      Warning("Not enough colours in palette, using rainbow().",
              "\nSet your colours with palette() first")
      line_cols <- pal[id_fac]
    } else {
      pal <- palette()
      line_cols <- pal[id_fac]
    }
    
  } else {
    line_cols <- rep("grey", length(x))
  }
  
  if(is.null(ylim)){
    amax <- max(sapply(x, function(x)max(x$df$Amodel)))
    amin <- max(sapply(x, function(x)min(x$df$Amodel)))
    ylim <- c(amin,amax)
  }
  if(is.null(xlim)){
    cimax <- max(sapply(x, function(x)max(x$df$Ci)))
    cimin <- min(sapply(x, function(x)min(x$df$Ci)))
    xlim <- c(cimin,cimax)
  }
  
  if(how == "manyplots"){
    if(add)Warning("Argument 'add' ignored when making multiple plots.")  
    
    for(i in seq_along(x)){
      plot.acifit(x[[i]],main=names(x)[i],xlim=xlim,ylim=ylim,...)
    }
  }
  
  if(how == "oneplot"){
    
    if(!is.null(highlight)){
      if(!highlight %in% names(x))
        stop("Curve ID not found.")
      
      hi <- which(names(x) == highlight)
      
      if(!add){
        plot.acifit(x[[1]], what="none", 
                    ylim=ylim, xlim=xlim, 
                    whichA="Amin", 
                    ...)
      }
      
      for(i in seq_along(x)){
        plot.acifit(x[[i]], what=what, whichA="Amin", add=TRUE,
                    linecols = line_cols[i], ...)  
      }
      plot.acifit(x[[hi]], what=what, whichA="Amin", add=TRUE,
                  linecols = linecol_highlight, ...)  
      
    } else {
      if(!add)
        plot.acifit(x[[1]], what="none", ylim=ylim, xlim=xlim, 
                    addlegend=FALSE,
                    whichA="Amin", ...)
      
      for(i in seq_along(x))
        plot.acifit(x[[i]], what=what, whichA="Amin", add=TRUE, 
                    linecols=line_cols[i], ...)  
      
    }
    
    if(colour_by_id && id_legend){
      legend("topleft", levels(id_fac), lty=1, col=pal, cex=0.8, lwd=2)
    }
    
  }
}


#' @method coef acifits
#' @export
coef.acifits <- function(object,...){
  
  get_pars <- function(object){
    if(all(is.na(object))) NA else as.vector(object$pars)
  }
  
  f <- lapply(object, get_pars)
  
  # Find objects without result (could not be fitted, even with bilinear),
  # and replace with contents of another fit, but all set to NA.
  # (This way, names and structure of coefficients is the same).
  ok <- sapply(f, function(x)!all(is.na(x)))
  if(any(!ok)){
    f[[which(!ok)]] <- f[[which(ok)[1]]]
    f[[which(!ok)]][] <- NA
  }
  
  pars <- as.data.frame(do.call(rbind,f))
  rn <- rownames(object[[which(ok)[1]]]$pars)
  nm <- c(rn, paste0(rn,"_SE"))
  names(pars) <- nm
  
  d <- data.frame(group=names(object))
  names(d) <- attr(object,"group")
  pars <- cbind(d,pars)
  rownames(pars) <- NULL
  
  if(!is.null(object[[1]]$id)){
    get_id <- function(x){
      res <- x$df[x$id]
      res1 <- res[1,,drop=FALSE]
      as.data.frame(lapply(res1, as.character))
    }
    ids <- do.call(rbind,lapply(object, get_id))
    pars <- cbind(pars, ids)
  }
  
  return(pars)
}










#'@export
#'@method print acifits
print.acifits <- function(x, ...){
  
  cat("Result of fitacis.\n\n")
  p <- coef(x)
  
  cat("Fitted", nrow(p), "curves by", attr(x, "groupname"), "grouping variable.")
  
  cat("\nRange in estimated Vcmax:", round(min(p$Vcmax, na.rm=TRUE),2), "-", round(max(p$Vcmax),2))
  cat("\nRange in estimated Jmax:", round(min(p$Jmax, na.rm=TRUE),2), "-", round(max(p$Jmax),2))
  cat("\nUse coef() on the object to see all fitted coefficients.")
  
}

