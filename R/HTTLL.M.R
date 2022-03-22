HTTLL.M.fun <- Vectorize(function(time, Psi, Temp, thetaHT, Tb, Psib50,
                                  Kt, delta, sigmaPsib){
  # Hydro-Thermal-Time model for seed germination
  # From Mohsen et al., 2017
  # log-logistic distribution of base water potential
  # Psib50 increases with temperature for any T > Tb
  # Rivedere le condizioni
  # Deve essere Psi > Psib50 + Kt*.cond1 e .germ2>0
  .cond1 <- ifelse(Temp > Tb, Temp - Tb, 0)
  .germ1 <- thetaHT/(.cond1 * time)
  .germ2 <- Psi - .germ1 + delta
  .germ2 <- ifelse(.germ2 < 0, 0.000001, .germ2)
  .germ2 <- log(.germ2)
  .germ3 <- Psib50 + delta + Kt * .cond1
  .germ3 <- ifelse(.germ3 < 0, 0.000001, .germ3)
  .germ3 <- log(.germ3)
  germ <- 1/(1 + exp(-(.germ2 - .germ3)/sigmaPsib))
  germ})

"HTTLL.M" <- function(){
fct <- function(x, parm){
  time <- x[,1]; Psi <- x[,2]; Temp <- x[,3]
  thetaHT <- parm[,1]; Tb <- parm[,2]; Psib50 <- parm[,3]
  Kt <- parm[,4]; delta <- parm[,5]; sigmaPsib <- parm[,6]
  HTTLL.M.fun(time, Psi, Temp, thetaHT, Tb, Psib50, Kt, delta, sigmaPsib)
}
names <- c("thetaHT", "Tb", "Psib50", "Kt", "delta", "sigmaPsib")
ss <- function(data){
    # thetaHT = 850; Tb=1; To = 30; Psib50=-2.5; Kt=0.07; delta=4; sigmaPsib=0.05
    # print(c(thetaHT, Tb, Psib50, Kt, delta, sigmaPsib))
    #
  # return(c(thetaHT, Tb, To, Psib50, Kt, delta, sigmaPsib))
}
deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d1.2 <- HTTLL.M.fun(x[,1], x[,2], x[,3],(parm[,1] + 10e-6), parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d2.2 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], (parm[,2] + 10e-6), parm[,3], parm[,4], parm[,5], parm[,6])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d3.2 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], (parm[,3] + 10e-6), parm[,4], parm[,5], parm[,6])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d4.2 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], (parm[,4] + 10e-6), parm[,5], parm[,6])
  d4 <- (d4.2 - d4.1)/10e-6

  d5.1 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d5.2 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], (parm[,5] + 10e-6), parm[,6])
  d5 <- (d5.2 - d5.1)/10e-6

  d6.1 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d6.2 <- HTTLL.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], (parm[,6] + 10e-6))
  d6 <- (d6.2 - d6.1)/10e-6

  cbind(d1, d2, d3, d4, d5, d6)
}
text <- "Hydrothermal-time model (Mesgaran et al., 2017)"
returnList <- list(fct=fct, names=names, text=text, deriv1 = deriv1) # ssfct=ss,
class(returnList) <- "drcMean"
invisible(returnList)
}
