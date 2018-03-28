GRPsiLin.fun <- function(Psi, Psib, thetaH) {
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

#From Mesgaran et al., 2017
GRT.M.fun <- function(predictor, Tc, Tb, ThetaT) {
  Temp <- predictor
  t2 <- ifelse(Temp < Tb, Tb, Temp)
  psival <- ifelse(1 - (Temp - Tb)/(Tc - Tb) > 0, 1 - (Temp - Tb)/(Tc - Tb), 0)
  GR <- psival * (t2 - Tb)/ThetaT
  return(ifelse(GR < 0 , 0 , GR)) }
"GRT.M" <- function(){
fct <- function(x, parm) {
  GR <- GRT.M.fun(x, parm[,1], parm[,2], parm[,3])
  return(GR)  }
names <- c("Tc", "Tb", "ThetaT")
ss <- function(data){
  x <- data[, 1]; y <- data[, 2]
  coefs <- coef( lm(y ~ x + I(x^2)))
  a <- coefs[3]; b <- coefs[2]; c <- coefs[1]
  Tb <- (-b + sqrt(b^2 - 4*a*c))/(2 * a)
  Tc <- (-b - sqrt(b^2 - 4*a*c))/(2 * a)
  ThetaT <- 1/b
  return(c(Tc, Tb, ThetaT))}
text <- "Polynomial temperature effect (Mohsen et al., 2017)"
returnList <- list(fct = fct, ssfct = ss, names = names, text = text)
class(returnList) <- "drcMean"
invisible(returnList)
}

#Rowse and Finch-Savage (with Tc and Td explicit, instead of k)
GRT.RFb.fun <- function(Temp, Tb, Td, Tc, ThetaT) {
  T2 <- ifelse(Temp < Tb, Tb, Temp)
  T1 <- ifelse(Temp < Td, Td, Temp)
  psival <- ifelse(1 - (T1 - Td)/(Tc - Td) > 0, 1 - (T1 - Td)/(Tc - Td), 0)
  GR <- psival * (T2 - Tb)/ThetaT
  return(ifelse(GR < 0 , 0 , GR)) }
"GRT.RFb" <- function(){
fct <- function(x, parm) {
  GR <- GRT.RFb.fun(x, parm[,1], parm[,2], parm[,3], parm[,4])
  return(ifelse(GR < 0 , 0 , GR))}
names <- c("Tb", "Td", "Tc", "ThetaT")
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
  Td <- - ss2[1] / k
  Tc <- Td + 1/k
  return(c(Tb, Td, Tc, ThetaT))}
text <- "Rowse - Finch-Savage model (derived from Rowse & Finch-Savage, 2003)"
returnList <- list(fct = fct, ssfct=ss, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}

# Exponential with switch-off (From Masin et al., 2017 - modified)
GRT.Exb.fun <- function(Temp, Tb, ThetaT, k, Tc) {
  GR50 <- ((Temp - Tb)/ThetaT) * ((1 - exp(k * (Temp - Tc)))/(1 - exp(k * (Tb - Tc))))
  return(ifelse(GR50 < 0 , 0 , GR50)) }
"GRT.Exb" <- function(){
fct <- function(x, parm) {
  GR50 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3], parm[,4])
  return(GR50) }
names <- c("Tb", "ThetaT", "k", "Tc")
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
  Tc <- (1 - ss2[1])/ss2[2]

  return(c(Tb, ThetaT, k, Tc))}
text <- "Exponential effect of temperature on GR50 (Type II - Masin et al., 2017)"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}
