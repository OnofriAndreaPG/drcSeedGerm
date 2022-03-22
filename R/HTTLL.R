HTTLL.fun <- function(time, Psi, Temp, thetaHT, Tb, Psib50, Kt, delta, sigmaPsib){
  # Hydro-Thermal-Time model for seed germination
  # From Mohsen et al., 2017
  # Rivedere le condizioni
  # Deve essere Psi > Psib50 + Kt*.cond1 e .germ2>0
  .cond1 <- ifelse(Temp > Tb, Temp - Tb, 0)
  .germ1 <- thetaHT/(.cond1 * time)
  #.cond2 <- Psi > Psib50 + Kt*.cond1
  .germ2 <- Psi - Kt*.cond1 - .germ1 + delta
  #.germ2 <- ifelse(.germ2>0, .germ2, 0.0000001)
  .germ3 <- .germ2/(Psib50 + delta)
  .germ3 <- ifelse(.germ3>0, .germ3, 0.0000001)
  .germ4 <- log(.germ3)
  germ <- 1/(1 + exp(-(.germ4/sigmaPsib)))
  germ
  #ifelse(Temp > Tb, germ, 0)
}

HTTLL2.fun <- function(time, Psi, Temp, thetaHT, Tb, To, Psib50, Kt, delta, sigmaPsib){
  # Hydro-Thermal-Time model for seed germination
  # From Mohsen et al., 2017
  # Rivedere le condizioni
  # Deve essere Psi > Psib50 + Kt*.cond1 e .germ2 > 0
  .cond1 <- ifelse(Temp > Tb, Temp - Tb, 0)
  .germ1 <- thetaHT/(.cond1 * time)
  .cond2 <- max(Temp, To)
  .germ2 <- Psi - .germ1 + delta
  .germ2 <- ifelse(.germ2 < 0, 0.000001, .germ2)
  .germ2 <- log(.germ2)
  .germ3 <- Psib50 + delta + Kt * .cond2
  .germ3 <- ifelse(.germ3 < 0, 0.000001, .germ3)
  .germ3 <- log(.germ3)
  germ <- 1/(1 + exp(-(.germ2 - .germ3)/sigmaPsib))
  germ
}


"HTTLL2" <- function(){
fct <- function(x, parm){
  time <- x[,1]; Psi <- x[,2]; Temp <- x[,3]
  thetaHT <- parm[,1]; Tb <- parm[,2]; To <- parm[,3]; Psib50 <- parm[,4]
  Kt <- parm[,5]; delta <- parm[,6]; sigmaPsib <- parm[,7]
  HTTLL2.fun(time, Psi, Temp, thetaHT, Tb, To, Psib50, Kt, delta, sigmaPsib)
}
names <- c("thetaHT", "Tb", "To", "Psib50", "Kt", "delta", "sigmaPsib")
ss <- function(data){
    thetaHT = 850; Tb=1; To = 30; Psib50=-2.5; Kt=0.07; delta=4; sigmaPsib=0.05
    # print(c(thetaHT, Tb, Psib50, Kt, delta, sigmaPsib))
    #
  return(c(thetaHT, Tb, To, Psib50, Kt, delta, sigmaPsib))
}
text <- "Hydrothermal-time model (Mesgaran et al., 2017)"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}
