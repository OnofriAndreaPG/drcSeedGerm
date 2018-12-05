##From: Mesgaran et al., 2013
HTLL.fun <- function(time, Psi, thetaH, delta, Psib50, sigma){
  .germ2 <- Psi - thetaH/time + delta
  .germ2 <- ifelse(.germ2 < 0, 0.000001, .germ2)
  .germ3 <- .germ2/(Psib50 + delta)
  germ <- 1/(1 + exp(-(log(.germ3)/sigma)))
  germ
}
"HTLL" <- function(){
fct <- function(x, parm){
  time <- x[,1]; Psi <- x[,2]
  thetaH <- parm[,1]; delta <- parm[,2]
  Psib50 <- parm[,3]; sigma <- parm[,4]
  HTLL.fun(time, Psi, thetaH, delta, Psib50, sigma)
}
text <- "Hydrotime model with log-logistic distribution of Psib (Mesgaran et al., 2013)"
names <- c("thetaH", "delta", "Psib50", "sigma")
ss <- function(data){
  x1 <- data[, 1]
  x2 <- data[, 2]
  y <- data[, 3]
  delta <- - (min(x2) - 0.05)
  pseudoY <- qnorm((y+10e-6)*0.99)
  mod <- lm(pseudoY ~ I(1/x1) + x2)
  sigma <- 1/coef(mod)[3]
  Psib50 <- -coef(mod)[1]*sigma
  thetaH <- -coef(mod)[2]*sigma
  return(c(thetaH, delta, Psib50, sigma))
}

GR <- function(parms, respl, reference="control", type="relative", Psi){
  HTLL.gra <- function(thetaH, delta, Psib50, sigma, Psi, g) {
    .temp1 <- sigma*(-log((1 - g)/g) ) + log(Psib50 + delta)
    .temp2 <- Psi + delta - exp(.temp1)
    GR <- .temp2 / thetaH
    GR <- ifelse(GR > 0, GR, 0)
  }
  thetaH <- as.numeric(parms[1])
  delta <- as.numeric(parms[2])
  Psib50 <- as.numeric(parms[3])
  sigma <- as.numeric(parms[4])
  g <- respl/100
  if(type=="absolute"){

    EDp <- HTLL.gra(thetaH, delta, Psib50, sigma, Psi, g)

    #Approximation of derivatives(finite differences)
    d1.1 <- HTLL.gra(thetaH, delta, Psib50, sigma, Psi, g)
    d1.2 <- HTLL.gra(thetaH + 10e-6, delta, Psib50, sigma, Psi, g)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HTLL.gra(thetaH, delta, Psib50, sigma, Psi, g)
    d2.2 <- HTLL.gra(thetaH, delta  + 10e-6, Psib50, sigma, Psi, g)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HTLL.gra(thetaH, delta, Psib50, sigma, Psi, g)
    d3.2 <- HTLL.gra(thetaH, delta, Psib50  + 10e-6, sigma, Psi, g)
    d3<- (d3.2 - d3.1)/10e-6

    d4.1 <- HTLL.gra(thetaH, delta, Psib50, sigma, Psi, g)
    d4.2 <- HTLL.gra(thetaH, delta, Psib50, sigma + 10e-6, Psi, g)
    d4 <- (d4.2 - d4.1)/10e-6

    EDder <- c(d1, d2, d3, d4)
  } else{ if(type=="relative") {
    .Pmax <- HTLL.fun(Inf, Psi, thetaH, delta, Psib50, sigma)
    grel <- .Pmax*g
    EDp <- HTLL.gra(thetaH, delta, Psib50, sigma, Psi, grel)

    #Approximation of derivatives(finite differences)
    d1.1 <- HTLL.gra(thetaH, delta, Psib50, sigma, Psi, grel)
    d1.2 <- HTLL.gra(thetaH + 10e-6, delta, Psib50, sigma, Psi, grel)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HTLL.gra(thetaH, delta, Psib50, sigma, Psi, grel)
    d2.2 <- HTLL.gra(thetaH, delta  + 10e-6, Psib50, sigma, Psi, grel)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HTLL.gra(thetaH, delta, Psib50, sigma, Psi, grel)
    d3.2 <- HTLL.gra(thetaH, delta, Psib50  + 10e-6, sigma, Psi, grel)
    d3<- (d3.2 - d3.1)/10e-6

    d4.1 <- HTLL.gra(thetaH, delta, Psib50, sigma, Psi, grel)
    d4.2 <- HTLL.gra(thetaH, delta, Psib50, sigma + 10e-6, Psi, grel)
    d4 <- (d4.2 - d4.1)/10e-6

    EDder <- c(d1, d2, d3, d4)
  } }
  return(list(EDp, EDder))
}

deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HTLL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4])
  d1.2 <- HTLL.fun(x[,1], x[,2], (parm[,1] + 10e-6), parm[,2], parm[,3], parm[4])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTLL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4])
  d2.2 <- HTLL.fun(x[,1], x[,2], parm[,1], (parm[,2] + 10e-6), parm[,3], parm[4])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTLL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4])
  d3.2 <- HTLL.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6), parm[4])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- HTLL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4])
  d4.2 <- HTLL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4] + 10e-6)
  d4 <- (d4.2 - d4.1)/10e-6

  cbind(d1, d2, d3, d4)
}

returnList <- list(fct=fct, ssfct=ss, names=names, text=text, edfct=GR, deriv1=deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
