#Models of GR vs water potential ####################

GRPsiLin.fun <- function(Psi, Psib, thetaH) {
  #Linear hydrotime model (Bradford, 2002)
  GR50 <- ( Psi - Psib ) / thetaH
  return(ifelse(GR50 < 0, 0, GR50)) }

"GRPsiLin" <- function() {
  fct <- function(x, parm)
  { GR50 <- GRPsiLin.fun(x, parm[,1], parm[,2])
    return(ifelse(GR50<0, 0, GR50)) }
  ssfct <- function(data)
  { x <- data[, 1]
  y <- data[, 2]
  mod <- lm( y ~ x )
  thetaH <- 1/coef(mod)[2]
  Psib <- - coef(mod)[2]*thetaH
  return(c(Psib, thetaH))
  }
  names <- c("Psib", "thetaH")
  text <- "Linear hydrotime model (Bradford, 2002)"
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

# second order polynomial
GRPsiPol.fun <- function(Psi, Psib, thetaH) {
  #"Polynomial hydrotime model - Convex up"
  GR50 <- - ( Psi^2 - Psib^2) / (thetaH)
  return(ifelse(GR50 < 0, 0, GR50)) }

"GRPsiPol" <- function(){
   GR502.fct <- function(x, parm) {
    Psi <- x; Psib <- parm[,1]; thetaH <- parm[,2]
    GR50 <- GRPsiPol.fun(Psi, Psib, thetaH)
    return(ifelse(GR50 < 0, 0, GR50)) }
  GR502.names <- c("Psib", "thetaH")
  GR502.ss <- function(data){
    x <- data[, 1]
    y <- data[, 2]
    isPositive <- y > 0
    x1 <- x[isPositive]
    y1 <- y[isPositive]
    mod <- lm(y1 ~ I(x1^2))
    thetaH <- - 1/coef(mod)[2]
    Psib <- - sqrt(thetaH * coef(mod)[1])
    return(c(Psib,thetaH))}
  GR502.text <- "Polynomial hydrotime model - Convex up"
  GR502 <- list(fct = GR502.fct, ssfct=GR502.ss, names=GR502.names, text = GR502.text)
  class(GR502) <- "drcMean"
  invisible(GR502)
}

# second order polynomial - 2
GRPsiPol2.fun <- function(Psi, Psib, thetaH) {
  #Polynomial hydrotime model - Convex down
  GR50 <- (( Psi - Psib)^2)/ (thetaH)
  return(ifelse(GR50 < 0, 0, GR50)) }
"GRPsiPol2" <- function(){
  fct <- function(x, parm) {
    Psi <- x; Psib <- parm[,1]; thetaH <- parm[,2]
    GR50 <- GRPsiPol2.fun(Psi, Psib, thetaH)
    return(ifelse(Psi < Psib, 0, GR50)) }
  names <- c("Psib", "thetaH")
  ss <- function(data){
    x <- data[, 1]
    y <- data[, 2]
    isPositive <- y > 0
    x1 <- x[isPositive]
    y1 <- y[isPositive]
    Psib <- x1[which( y1==min(y1) )]
    mod <- lm(y1 ~ I((x1 - Psib)^2) - 1)
    thetaH <- 1/coef(mod)[1]
    return(c(Psib,thetaH))}
  text <- "Polynomial hydrotime model - Convex down"
  returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

#Models for GR vs Temperature #########################
GRT.GH.fun <- function(Temp, Tb, ThetaT){
  #Garcia-Huidobro, 1982 - Linear model for suboptimal temperatures
  t2 <- ifelse(Temp < Tb, Tb, Temp)
  GR <- (t2 - Tb)/ThetaT
  GR }
"GRT.GH" <- function() {
  fct <- function(x, parm) {
    Tb <- parm[,1]; ThetaT <- parm[,2]
    GR <- GRT.GH.fun(x, Tb, ThetaT)
    return(GR) }
  names <- c("Tb", "ThetaT")
  ssfct <- function(data){
    x <- data[, 1]
    y <- data[, 2]
    ss1 <- coef( lm(y ~ x) )
    ThetaT <- 1/ss1[2]
    Tb <- - ss1[1] * ThetaT
    return(c(Tb, ThetaT))}
  text <- "Linear model for suboptimal temperatures (Garcia-Huidobro, 1982)"
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}



#From Bradford 2002 - Broken-Stick Model
GRT.BS.fun <- function(Temp, k, Tb, To, ThetaT){
  t2 <- ifelse(Temp < Tb, Tb, ifelse(Temp > To, To, Temp))
  t1 <- ifelse(Temp < To, To, Temp)
  psival <- ifelse(1 - k*(t1 - To) > 0, 1 - k*(t1 - To), 0)
  GR <- psival * (t2 - Tb)/ThetaT
  GR }
"GRT.BS" <- function(){
fct <- function(x, parm) {
  k <- parm[,1]; Tb <- parm[,2]; To <- parm[,3]; ThetaT <- parm[,4]
  GR <- GRT.BS.fun(x, k, Tb, To, ThetaT)
  return(ifelse(GR < 0 , 0 , GR)) }
names <- c("k", "Tb", "To", "ThetaT")
ss <- function(data){
  pos <- which( data[,2]==max(data[,2]) )
  len <- length( data[,2] )

  reg1 <- data[1:pos, ]
  reg2 <- data[pos:len, ]
  x1 <- reg1[,1]; y1 <- reg1[, 2]
  x2 <- reg2[,1]; y2 <- reg2[, 2]

  ss1 <- coef( lm(y1 ~ x1) )
  ThetaT <- 1/ss1[2]
  Tb <- - ss1[1] * ThetaT
  ss2 <- coef( lm((1-y2) ~ x2) )
  k <- ss2[2]
  To <- - ss2[1] / k

  #k <- 0.1; Tb <- 2; To <- 20; ThetaT <- 35
  return(c(k, Tb, To, ThetaT))}
text <- "Broken-stick model (Bradford, 2002)"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}

#Yield loss function - derived from Kroppf, competition
GRT.YL.fun <- function(Temp, q, tmin, tmax, theta) {
  t <- Temp
  dt <- ifelse(t > tmin & t < tmax, t - tmin, 0)
  GR <- (dt/theta) * (1 - q*( dt/(tmax - tmin))/(1 + (q - 1)*(dt/(tmax - tmin))))
  return(ifelse(GR < 0 , 0 , GR)) }
"GRT.YL" <- function(){
fct <- function(x, parm) {
  t <- x
  q <- parm[,1]; tmin <- parm[,2]; tmax <- parm[,3]; theta <- parm[,4]
  GR <- GRT.YL.fun(t, q, tmin, tmax, theta)
  return(GR)}

names <- c("k", "Tb", "Tc", "ThetaT")
ss <- function(data){
  pos <- which( data[,2]==max(data[,2]) )
  len <- length( data[,2] )

  reg1 <- data[1:pos, ]
  reg2 <- data[pos:len, ]
  x1 <- reg1[,1]; y1 <- reg1[, 2]
  x2 <- reg2[,1]; y2 <- reg2[, 2]

  ss1 <- coef( lm(y1 ~ x1) )
  ThetaT <- 1/ss1[2]
  Tb <- - ss1[1] * ThetaT
  ss2 <- coef( lm((1-y2) ~ x2) )
  k <- ss2[2]
  td <- - ss2[1] / k
  Tc <- td + 1/k
  return(c(k, Tb, Tc, ThetaT))}
text <- "Yield-loss function (Kropff and van Laar, 1993)"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}
