#Hydrotime models ######################################
HTL.fun <- function(time, Psi, ThetaH, Psib50, sigma){
  plogis((Psi - (ThetaH/time) - Psib50)/sigma) }

"HTL" <- function(){
fct <- function(x, parm){
  HTL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3]) }
names <- c("ThetaH", "Psib50", "sigma")
text <- "Hydrotime model with logistic distribution of Psib (Mesgaran et al., 2013)"
ss <- function(data){
  x1 <- data[, 1]
  x2 <- data[, 2]
  y <- data[, 3]
  pseudoY <- qnorm((y+10e-6)*0.99)
  mod <- lm(pseudoY ~ I(1/x1) + x2)
  sigmaPsib <- 1/coef(mod)[3]
  Psib50 <- -coef(mod)[1]*sigmaPsib
  ThetaH <- -coef(mod)[2]*sigmaPsib
  return(c(ThetaH, Psib50, sigmaPsib))
}
GR <- function(parms, respl, reference="control", type="relative", Psi){
  HTL.gra <- function(thetaH, Psib50, sigma, Psi, g) {
    GR <- - (sigma * qlogis(g) - Psi + Psib50 )/ thetaH
    GR <- ifelse(GR > 0, GR, 0)
  }
  thetaH <- as.numeric(parms[1])
  Psib50 <- as.numeric(parms[2])
  sigma <- as.numeric(parms[3])
  g <- respl/100
  if(type=="absolute"){

    EDp <- HTL.gra(thetaH, Psib50, sigma, Psi, g)

    #Approximation of derivatives(finite differences)
    d1.1 <- HTL.gra(thetaH, Psib50, sigma, Psi, g)
    d1.2 <- HTL.gra(thetaH + 10e-6, Psib50, sigma, Psi, g)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HTL.gra(thetaH, Psib50, sigma, Psi, g)
    d2.2 <- HTL.gra(thetaH, Psib50  + 10e-6, sigma, Psi, g)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HTL.gra(thetaH, Psib50, sigma, Psi, g)
    d3.2 <- HTL.gra(thetaH, Psib50, sigma + 10e-6, Psi, g)
    d3 <- (d3.2 - d3.1)/10e-6
    EDder <- c(d1, d2, d3)
  } else{ if(type=="relative") {
    .Pmax <- plogis((Psi - Psib50 )/sigma)
    grel <- .Pmax*g
    EDp <- HTL.gra(thetaH, Psib50, sigma, Psi, grel)

    #Approximation of derivatives(finite differences)
    d1.1 <- HTL.gra(thetaH, Psib50, sigma, Psi, grel)
    d1.2 <- HTL.gra(thetaH + 10e-6, Psib50, sigma, Psi, grel)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HTL.gra(thetaH, Psib50, sigma, Psi, grel)
    d2.2 <- HTL.gra(thetaH, Psib50  + 10e-6, sigma, Psi, grel)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HTL.gra(thetaH, Psib50, sigma, Psi, grel)
    d3.2 <- HTL.gra(thetaH, Psib50, sigma + 10e-6, Psi, grel)
    d3 <- (d3.2 - d3.1)/10e-6
    EDder <- c(d1, d2, d3)
  } }
  return(list(EDp, EDder))
}

deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HTL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3])
  d1.2 <- HTL.fun(x[,1], x[,2], (parm[,1] + 10e-6), parm[,2], parm[,3])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3])
  d2.2 <- HTL.fun(x[,1], x[,2], parm[,1], (parm[,2] + 10e-6), parm[,3])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3])
  d3.2 <- HTL.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6))
  d3 <- (d3.2 - d3.1)/10e-6

  cbind(d1, d2, d3)
}

returnList <- list(fct=fct, ssfct=ss, names=names, text=text, edfct=GR, deriv1=deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
