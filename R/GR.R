GR <- function(mod, Psi, respLev, type="absolute", vcov.){
  #mod <- obj; x <- Psi; respLev <- g
  GR <- ED(mod, respLev=respLev, type=type, Psi=Psi, vcov. = vcovCL, display=F)
  row.names(GR) <- gsub("e:", "GR:", row.names(GR), fixed=T)
  GR
  #data.frame(Estimate=GR[,1], SE=GR[,2], row.names = paste("GR", respLev, sep=":"))
  }

GTime <- function(mod, Psi, respLev, type="absolute", vcov.){
  #mod <- obj; x <- Psi; respLev <- g
  GT <- ED(mod, respLev=respLev, type=type, Psi=Psi, vcov. = vcovCL, display=F)
  GTval <- 1/GT[,1]
  GTes <- (GT[,2] * 1/GT[,1]^2)
  returnDF <- data.frame(Estimate=GTval, SE=GTes)
  row.names(returnDF) <- gsub("e:", "T:", row.names(returnDF), fixed=T)
  returnDF
  #row.names = paste("T", respLev, sep=":")
  }
