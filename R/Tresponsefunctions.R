TGammaStar <- function(Tleaf, 
                       Egamma=37830.0, 
                       value25=42.75){  

  value25*arrh(Tleaf,Egamma)
  
}

TKm <- function(Tleaf,
                Oi = 210,      # O2 concentration
                Ec = 79430.0,  # activation energy for Kc
                Eo = 36380.0,  # activation energy for Ko
                Kc25 = 404.9,  # Kc at 25C
                Ko25 = 278.4  # Ko at 25C
                ){
  
  Ko <- Ko25*arrh(Tleaf, Eo)
  Kc <- Kc25*arrh(Tleaf, Ec)
  Km <- Kc * (1.0 + Oi / Ko)
  
return(Km)
}


