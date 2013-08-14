.Rgas <- function()8.314

# Arrhenius
arrh <- function(Tleaf, Ea){
  Tk <- Tleaf + 273.15
  exp((Ea * (Tk - 298.15)) / (298.15 * .Rgas() * Tk)) 
}

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

# Hard-wired parameters.
TVcmax <- function(Tleaf, EaV, delsC, EdVC){
  
  if(EdVC > 0){
    V1 <- (1+exp((delsC*(25 + 273.15)-EdVC)/(.Rgas()*(25 + 273.15))))
    V2 <- (1+exp((delsC*(Tleaf+273.15)-EdVC)/(.Rgas()*(Tleaf+273.15))))
    f <- V1/V2
  } else f <- 1
  
  exp((Tleaf-25)*EaV/(.Rgas()*(Tleaf+273.15)*(25 + 273.15))) * f
}

# Hard-wired parameters.
TJmax <- function(Tleaf, EaJ, delsJ, EdVJ){
  J1 <- 1+exp((298.15*delsJ-EdVJ)/.Rgas()/298.15)
  J2 <- 1+exp(((Tleaf+273.15)*delsJ-EdVJ)/.Rgas()/(Tleaf+273.15))
  exp(EaJ/.Rgas()*(1/298.15 - 1/(Tleaf+273.15)))*J1/J2
}





