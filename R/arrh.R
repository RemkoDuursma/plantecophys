# Arrhenius
arrh <- function(Tleaf, Ea){
  Tk <- Tleaf + 273.15
  exp((Ea * (Tk - 298.15)) / (298.15 * Rgas * Tk)) 
}