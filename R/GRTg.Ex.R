# Exponential with switch-off (From Masin et al., 2017 - modified)
GRTg.Exb.fun <- function(Temp, g, Tb, ThetaT50, k, Tc50, b1, b2) {
  ThetaTg <- ThetaT50 * ( ( (1 - g)/g ) ^ (-1/b1) )
  Tcg <- Tc50 * ( ( (1 - g)/g ) ^ (-1/b2) )
  GR <- ((Temp - Tb)/ThetaTg) * ((1 - exp(k * (Temp - Tcg)))/(1 - exp(k * (Tb - Tcg))))
  return(ifelse(GR < 0 , 0, GR)) }

"GRTg.Exb" <- function(){
fct <- function(x, parm) {
  GR50 <- GRTg.Exb.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  return(GR50) }
names <- c("Tb", "ThetaT50", "k", "Tc50", "b1", "b2")
ss <- function(data){
  Tb <- -1
  ThetaT50 <- 51
  k <- 1.28
  Tc50 <- 31
  b1 <- -5.4
  b2 <- 13
return(c(Tb, ThetaT50, k, Tc50, b1, b2))}

deriv1 <- function(x, parm){
  #Approximation by using finite differences
  # d1.1 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3],
  #                  parm[,4])
  # d1.2 <- GRT.Exb.fun(x, (parm[,1] + 10e-6), parm[,2], parm[,3],
  #                  parm[,4])
  # d1 <- (d1.2 - d1.1)/10e-6
  #
  # d2.1 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3],
  #                  parm[,4])
  # d2.2 <- GRT.Exb.fun(x, parm[,1], (parm[,2] + 10e-6), parm[,3],
  #                  parm[,4])
  # d2 <- (d2.2 - d2.1)/10e-6
  #
  # d3.1 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3],
  #                  parm[,4])
  # d3.2 <- GRT.Exb.fun(x, parm[,1], parm[,2], (parm[,3] + 10e-6),
  #                  parm[,4])
  # d3 <- (d3.2 - d3.1)/10e-6
  #
  # d4.1 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3],
  #                  parm[,4])
  # d4.2 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3],
  #                  (parm[,4] + 10e-6))
  # d4 <- (d4.2 - d4.1)/10e-6
  #
  # cbind(d1, d2, d3, d4)
}

text <- "Exponential effect of temperature on GR50 (Type II - Masin et al., 2017)"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text, deriv1 = deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}

