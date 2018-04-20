##From: Mesgaran et al., 2013
HTW2.fun <- function(time, Psi, thetaH, delta, mu, sigma){
  .germ2 <- Psi - thetaH/time + delta
  .germ2 <- ifelse(.germ2 < 0, 0.000001, .germ2)
  .germ3 <- .germ2/(mu + delta)
  germ <- 1 - exp(-exp((log(.germ3)/sigma)) )
  germ
}
"HTW2" <- function(){
fct <- function(x, parm){
  time <- x[,1]; Psi <- x[,2]
  thetaH <- parm[,1]; delta <- parm[,2]
  mu <- parm[,3]; sigma <- parm[,4]
  HTW2.fun(time, Psi, thetaH, delta, mu, sigma)
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
  Psib50 <- -coef(mod)[1]*sigmaPsib
  thetaH <- -coef(mod)[2]*sigmaPsib
  return(c(thetaH, delta, Psib50, sigmaPsib))
}
text <- "Hydrotime model with Weibull type II distribution of psib"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}
