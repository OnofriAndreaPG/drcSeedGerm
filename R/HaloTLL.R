# Halotime model with log-logistic distribution of base salt concentration
# 9/1/2024
HaloTLL.fun <- function(time, SConc, thetaHalo, SConcb50, sigma){
  .germ2 <- SConc + thetaHalo/time
  # .germ2 <- ifelse(.germ2 < 0, 0.000001, .germ2)
  .germ3 <- .germ2/SConcb50
  germ <- 1 - 1/(1 + exp(-(log(.germ3)/sigma)))
  germ
  # plogis(log(.germ3)/sigma, lower.tail = F)
}
"HaloTLL" <- function(){
fct <- function(x, parm){
  time <- x[,1]; SConc <- x[,2]
  thetaHalo <- parm[,1]
  SConcb50 <- parm[,2]; sigma <- parm[,3]
  HaloTLL.fun(time, SConc, thetaHalo, SConcb50, sigma)
}
text <- "Halotime model with log-logistic distribution of SConcb"
names <- c("thetaHalo", "SConcb50", "sigma")
name <- "HaloTLL"
ss <- function(data){
  x1 <- data[, 1]
  x2 <- data[, 2]
  y <- data[, 3]
  # delta <- - (min(x2) - 0.05)
  pseudoY <- qnorm((y+10e-6)*0.99, lower.tail = F)
  mod <- lm(pseudoY ~ I(1/x1) + x2)
  sigma <- 1/coef(mod)[3]
  SConcb50 <- -coef(mod)[1]*sigma
  thetaHalo <- coef(mod)[2]*sigma
  return(c(thetaHalo, SConcb50, sigma))
}

GT <- function(parms, respl, reference="control", type="relative", SConc){
    # This function produces the quantiles for times-to-event
    # Respl is the quantile on a relative value (from 0 to 1)
    HaloTLL.gra <- function(thetaHalo, SConcb50, sigma, SConc, g) {
    .temp1 <- sigma*(-log( - (g)/(g - 1) )) + log(SConcb50)
    .temp2 <- exp(.temp1) - SConc
    GR <- .temp2 / thetaHalo
    GR <- ifelse(GR > 0, GR, 0)
    1/GR # returns time
  }
  thetaHalo <- as.numeric(parms[1])
  SConcb50 <- as.numeric(parms[2])
  sigma <- as.numeric(parms[3])
  # g <- respl #/100
  g <- respl
  if(type=="absolute"){

    EDp <- HaloTLL.gra(thetaHalo, SConcb50, sigma, SConc, g)

    #Approximation of derivatives(finite differences)
    d1.1 <- HaloTLL.gra(thetaHalo, SConcb50, sigma, SConc, g)
    d1.2 <- HaloTLL.gra(thetaHalo + 10e-6, SConcb50, sigma, SConc, g)
    d1 <- (d1.2 - d1.1)/10e-6

    # d2.1 <- HaloTLL.gra(thetaHalo, delta, SConcb50, sigma, SConc, g)
    # d2.2 <- HaloTLL.gra(thetaHalo, delta  + 10e-6, SConcb50, sigma, SConc, g)
    # d2 <- (d2.2 - d2.1)/10e-6

    d2.1 <- HaloTLL.gra(thetaHalo, SConcb50, sigma, SConc, g)
    d2.2 <- HaloTLL.gra(thetaHalo, SConcb50  + 10e-6, sigma, SConc, g)
    d2<- (d2.2 - d2.1)/10e-6

    d3.1 <- HaloTLL.gra(thetaHalo, SConcb50, sigma, SConc, g)
    d3.2 <- HaloTLL.gra(thetaHalo, SConcb50, sigma + 10e-6, SConc, g)
    d3 <- (d3.2 - d3.1)/10e-6

    EDder <- c(d1, d2, d3)
  } else{ if(type=="relative") {
    .Pmax <- HaloTLL.fun(Inf, SConc, thetaHalo, delta, SConcb50, sigma)
    grel <- .Pmax*g
    EDp <- HaloTLL.gra(thetaHalo, SConcb50, sigma, SConc, grel)

    #Approximation of derivatives(finite differences)
    d1.1 <- HaloTLL.gra(thetaHalo, SConcb50, sigma, SConc, grel)
    d1.2 <- HaloTLL.gra(thetaHalo + 10e-6, SConcb50, sigma, SConc, grel)
    d1 <- (d1.2 - d1.1)/10e-6

    # d2.1 <- HaloTLL.gra(thetaHalo, delta, SConcb50, sigma, SConc, grel)
    # d2.2 <- HaloTLL.gra(thetaHalo, delta  + 10e-6, SConcb50, sigma, SConc, grel)
    # d2 <- (d2.2 - d2.1)/10e-6

    d2.1 <- HaloTLL.gra(thetaHalo, SConcb50, sigma, SConc, grel)
    d2.2 <- HaloTLL.gra(thetaHalo, SConcb50  + 10e-6, sigma, SConc, grel)
    d2<- (d2.2 - d2.1)/10e-6

    d3.1 <- HaloTLL.gra(thetaHalo, SConcb50, sigma, SConc, grel)
    d3.2 <- HaloTLL.gra(thetaHalo, SConcb50, sigma + 10e-6, SConc, grel)
    d3 <- (d3.2 - d3.1)/10e-6

    EDder <- c(d1, d2, d3)
  } }
  return(list(EDp, EDder))
}

deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HaloTLL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3]) #, parm[4])
  d1.2 <- HaloTLL.fun(x[,1], x[,2], (parm[,1] + 10e-6), parm[,2], parm[,3]) #, parm[4])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HaloTLL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3]) #, parm[4])
  d2.2 <- HaloTLL.fun(x[,1], x[,2], parm[,1], (parm[,2] + 10e-6), parm[,3]) #, parm[4])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HaloTLL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3]) #, parm[4])
  d3.2 <- HaloTLL.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6)) #, parm[4])
  d3 <- (d3.2 - d3.1)/10e-6

  # d4.1 <- HaloTLL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4])
  # d4.2 <- HaloTLL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[4] + 10e-6)
  # d4 <- (d4.2 - d4.1)/10e-6

  cbind(d1, d2, d3) #, d4)
}

returnList <- list(fct=fct, ssfct=ss, names=names, text=text,
                   edfct = GT, deriv1=deriv1, name = name)
class(returnList) <- "drcMean"
invisible(returnList)
}
