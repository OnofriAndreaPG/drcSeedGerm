##From: Mesgaran et al., 2013
HTW1.fun <- function(time, Psi, thetaH, delta, mu, sigma){
  .germ2 <- Psi - thetaH/time + delta
  .germ2 <- ifelse(.germ2 < 0, 0.000001, .germ2)
  .germ3 <- .germ2/(mu + delta)
  germ <- exp(-exp(-(log(.germ3)/sigma)))
  germ
}
"HTW1" <- function(){
fct <- function(x, parm){
  time <- x[,1]; Psi <- x[,2]
  thetaH <- parm[,1]; delta <- parm[,2]
  mu <- parm[,3]; sigmaPsib <- parm[,4]
  HTW1.fun(time, Psi, thetaH, delta, mu, sigmaPsib)
}
names <- c("thetaH", "delta", "mu", "sigma")
ss <- function(data){
  x1 <- data[, 1]
  x2 <- data[, 2]
  y <- data[, 3]
  delta <- - (min(x2)-0.05)
  pseudoY <- qnorm((y+10e-6)*0.99)
  mod <- lm(pseudoY ~ I(1/x1) + x2)
  sigmaPsib <- 1/coef(mod)[3]
  mu <- -coef(mod)[1]*sigmaPsib
  thetaH <- -coef(mod)[2]*sigmaPsib
  return(c(thetaH, delta, mu, sigmaPsib))
}
text <- "Hydrotime model with Weibull (Type I) distribution of Psib (Mesgaran et al., 2013)"
GR <- function(parms, respl, reference="control", type="relative", Psi){
  HTW1.gra <- function(thetaH, delta, mu, sigma, Psi, g) {
    .temp1 <- sigma*(-log(-log(g)) ) + log(mu + delta)
    .temp2 <- Psi + delta - exp(.temp1)
    GR <- .temp2 / thetaH
    GR <- ifelse(GR > 0, GR, 0)
  }
  thetaH <- as.numeric(parms[1])
  delta <- as.numeric(parms[2])
  mu <- as.numeric(parms[3])
  sigma <- as.numeric(parms[4])
  g <- respl/100
  if(type=="absolute"){

    EDp <- HTW1.gra(thetaH, delta, mu, sigma, Psi, g)

    #Approximation of derivatives(finite differences)
    d1.1 <- HTW1.gra(thetaH, delta, mu, sigma, Psi, g)
    d1.2 <- HTW1.gra(thetaH + 10e-6, delta, mu, sigma, Psi, g)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HTW1.gra(thetaH, delta, mu, sigma, Psi, g)
    d2.2 <- HTW1.gra(thetaH, delta  + 10e-6, mu, sigma, Psi, g)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HTW1.gra(thetaH, delta, mu, sigma, Psi, g)
    d3.2 <- HTW1.gra(thetaH, delta, mu  + 10e-6, sigma, Psi, g)
    d3<- (d3.2 - d3.1)/10e-6

    d4.1 <- HTW1.gra(thetaH, delta, mu, sigma, Psi, g)
    d4.2 <- HTW1.gra(thetaH, delta, mu, sigma + 10e-6, Psi, g)
    d4 <- (d4.2 - d4.1)/10e-6

    EDder <- c(d1, d2, d3, d4)
  } else{ if(type=="relative") {
    .Pmax <- HTW1.fun(Inf, Psi, thetaH, delta, mu, sigma)
    grel <- .Pmax*g
    EDp <- HTW1.gra(thetaH, delta, mu, sigma, Psi, grel)

    #Approximation of derivatives(finite differences)
    d1.1 <- HTW1.gra(thetaH, delta, mu, sigma, Psi, grel)
    d1.2 <- HTW1.gra(thetaH + 10e-6, delta, mu, sigma, Psi, grel)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HTW1.gra(thetaH, delta, mu, sigma, Psi, grel)
    d2.2 <- HTW1.gra(thetaH, delta  + 10e-6, mu, sigma, Psi, grel)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HTW1.gra(thetaH, delta, mu, sigma, Psi, grel)
    d3.2 <- HTW1.gra(thetaH, delta, mu  + 10e-6, sigma, Psi, grel)
    d3<- (d3.2 - d3.1)/10e-6

    d4.1 <- HTW1.gra(thetaH, delta, mu, sigma, Psi, grel)
    d4.2 <- HTW1.gra(thetaH, delta, mu, sigma + 10e-6, Psi, grel)
    d4 <- (d4.2 - d4.1)/10e-6

    EDder <- c(d1, d2, d3, d4)
  } }
  return(list(EDp, EDder))
}

deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HTW1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4])
  d1.2 <- HTW1.fun(x[,1], x[,2], (parm[,1] + 10e-6), parm[,2], parm[,3], parm[4])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTW1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4])
  d2.2 <- HTW1.fun(x[,1], x[,2], parm[,1], (parm[,2] + 10e-6), parm[,3], parm[4])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTW1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4])
  d3.2 <- HTW1.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6), parm[4])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- HTW1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4])
  d4.2 <- HTW1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4] + 10e-6)
  d4 <- (d4.2 - d4.1)/10e-6

  cbind(d1, d2, d3, d4)
}

returnList <- list(fct=fct, ssfct=ss, names=names, text=text, edfct=GR, deriv1=deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
