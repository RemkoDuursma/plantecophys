#'@export
Photosyn_hack <- function(Tleaf=25, ...){
  
  topt <- find_topt(...)
  Aopt <- Photosyn(Tleaf=topt, ...)$ALEAF
  
  l <- list()
  for(i in seq_along(Tleaf)){
    
    if(Tleaf[i] > topt){
      l[[i]] <- Photosyn(Tleaf=Tleaf[i], whichA="Aval", Aval=Aopt, ...)  
    } else {
      l[[i]] <- Photosyn(Tleaf=Tleaf[i], ...)  
    }
    
    
  }
  
  dfr <- do.call(rbind,l)
  dfr$Aopt <- Aopt
  dfr$Topt <- topt
  
  return(dfr)
}


find_topt <- function(...){
  
  O <- function(tleaf)Photosyn(Tleaf=tleaf, ...)$ALEAF
  optimize(O, c(15,40), maximum=TRUE)$maximum
  
}


