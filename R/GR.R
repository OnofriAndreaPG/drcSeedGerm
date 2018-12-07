GR <- function(mod, Psi, respLev, type="absolute", vcov.){
  #mod <- obj; x <- Psi; respLev <- g
  GR <- ED(mod, respLev=respLev, type=type, Psi=Psi, vcov. = vcovCL, display=F)
  data.frame(Estimate=GR[,1], SE=GR[,2], row.names = paste("GR", respLev, sep=":"))
  }

GTime <- function(mod, Psi, respLev, type="absolute", vcov.){
  #mod <- obj; x <- Psi; respLev <- g
  GT <- ED(mod, respLev=respLev, type=type, Psi=Psi, vcov. = vcovCL, display=F)
  GTval <- 1/GT[,1]
  GTes <- (GT[,2] * 1/GT[,1]^2)
  data.frame(Estimate=GTval, SE=GTes, row.names = paste("T", respLev, sep=":"))
}
