# Halotime model with logistic distribution of base salt concentration
# 9/1/2024
HaloTL.fun <- function(time, SConc, ThetaHalo, SConcb50, sigma){
  plogis((SConc + (ThetaHalo/time) - SConcb50)/sigma, lower.tail = F) }

"HaloTL" <- function(){
fct <- function(x, parm){
  HaloTL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3]) }
names <- c("ThetaHalo", "SConcb50", "sigma")
text <- "Halotime model with logistic distribution of SConcb"
name <- "HaloTLL"
ss <- function(data){
  x1 <- data[, 1]
  x2 <- data[, 2]
  y <- data[, 3]
  pseudoY <- qnorm((y+10e-6)*0.99, lower.tail = F)
  mod <- lm(pseudoY ~ I(1/x1) + x2)
  sigma <- 1/coef(mod)[3]
  SConcb50 <- -coef(mod)[1]*sigma
  ThetaHalo <- coef(mod)[2]*sigma
  return(c(ThetaHalo, SConcb50, sigma))
}
GR <- function(parms, respl, reference="control", type="relative", SConc){
  HaloTL.gra <- function(ThetaHalo, SConcb50, sigma, SConc, g) {
    GR <- (sigma * qlogis(g, lower.tail = F) - SConc + SConcb50 )/ ThetaHalo #returns rate
    GR <- ifelse(GR > 0, GR, 0)
    1/GR # returns time
  }
  ThetaHalo <- as.numeric(parms[1])
  SConcb50 <- as.numeric(parms[2])
  sigma <- as.numeric(parms[3])
  # g <- respl/100 # bug corrected 9/1/25
  g <- respl
  if(type=="absolute"){

    EDp <- HaloTL.gra(ThetaHalo, SConcb50, sigma, SConc, g)

    #Approximation of derivatives(finite differences)
    d1.1 <- HaloTL.gra(ThetaHalo, SConcb50, sigma, SConc, g)
    d1.2 <- HaloTL.gra(ThetaHalo + 10e-6, SConcb50, sigma, SConc, g)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HaloTL.gra(ThetaHalo, SConcb50, sigma, SConc, g)
    d2.2 <- HaloTL.gra(ThetaHalo, SConcb50  + 10e-6, sigma, SConc, g)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HaloTL.gra(ThetaHalo, SConcb50, sigma, SConc, g)
    d3.2 <- HaloTL.gra(ThetaHalo, SConcb50, sigma + 10e-6, SConc, g)
    d3 <- (d3.2 - d3.1)/10e-6
    EDder <- c(d1, d2, d3)
  } else{ if(type=="relative") {
    .Pmax <- plogis((SConc - SConcb50 )/sigma, lower.tail = F)
    grel <- .Pmax*g
    EDp <- HaloTL.gra(ThetaHalo, SConcb50, sigma, SConc, grel)

    #Approximation of derivatives(finite differences)
    d1.1 <- HaloTL.gra(ThetaHalo, SConcb50, sigma, SConc, grel)
    d1.2 <- HaloTL.gra(ThetaHalo + 10e-6, SConcb50, sigma, SConc, grel)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HaloTL.gra(ThetaHalo, SConcb50, sigma, SConc, grel)
    d2.2 <- HaloTL.gra(ThetaHalo, SConcb50  + 10e-6, sigma, SConc, grel)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HaloTL.gra(ThetaHalo, SConcb50, sigma, SConc, grel)
    d3.2 <- HaloTL.gra(ThetaHalo, SConcb50, sigma + 10e-6, SConc, grel)
    d3 <- (d3.2 - d3.1)/10e-6
    EDder <- c(d1, d2, d3)
  } }
  return(list(EDp, EDder))
}

deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HaloTL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3])
  d1.2 <- HaloTL.fun(x[,1], x[,2], (parm[,1] + 10e-6), parm[,2], parm[,3])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HaloTL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3])
  d2.2 <- HaloTL.fun(x[,1], x[,2], parm[,1], (parm[,2] + 10e-6), parm[,3])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HaloTL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3])
  d3.2 <- HaloTL.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6))
  d3 <- (d3.2 - d3.1)/10e-6

  cbind(d1, d2, d3)
}
returnList <- list(fct=fct, ssfct=ss, name = name, names=names, text=text, edfct=GR, deriv1=deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
