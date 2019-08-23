#From Bradford 2002 - Broken-Stick Model
GRT.BS.fun <- function(Temp, Tc, Tb, To, ThetaT){
  t2 <- ifelse(Temp < Tb, Tb, ifelse(Temp > To, To, Temp))
  t1 <- ifelse(Temp < To, To, Temp)
  psival <- ifelse(1 - (t1 - To)/(Tc - To) > 0, 1 - (t1 - To)/(Tc - To), 0)
  GR <- psival * (t2 - Tb)/ThetaT
  GR }

"GRT.BS" <- function(){
fct <- function(x, parm) {
  Tc <- parm[,1]; Tb <- parm[,2]; To <- parm[,3]; ThetaT <- parm[,4]
  GR <- GRT.BS.fun(x, Tc, Tb, To, ThetaT)
  return(ifelse(GR < 0 , 0 , GR)) }
names <- c("Tc", "Tb", "To", "ThetaT")
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
  Tc <- To + 1/k

  #k <- 0.1; Tb <- 2; To <- 20; ThetaT <- 35
  return(c(Tc, Tb, To, ThetaT))}
deriv1 <- function(x, parms){

    #Approximation by using finite differences
    Temp <- x
    Tc <-  as.numeric(parms[,1]); Tb <- as.numeric(parms[,2]); To <- as.numeric(parms[,3])
    ThetaT <- as.numeric(parms[,4])

    d1.1 <- GRT.BS.fun(Temp, Tc, Tb, To, ThetaT)
    d1.2 <- GRT.BS.fun(Temp, (Tc + 10e-6), Tb, To, ThetaT)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- GRT.BS.fun(Temp, Tc, Tb, To, ThetaT)
    d2.2 <- GRT.BS.fun(Temp, Tc, (Tb + 10e-6), To, ThetaT)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- GRT.BS.fun(Temp, Tc, Tb, To, ThetaT)
    d3.2 <- GRT.BS.fun(Temp, Tc, Tb, (To + 10e-6), ThetaT)
    d3 <- (d3.2 - d3.1)/10e-6

    d4.1 <- GRT.BS.fun(Temp, Tc, Tb, To, ThetaT)
    d4.2 <- GRT.BS.fun(Temp, Tc, Tb, To, (ThetaT + 10e-6))
    d4 <- (d4.2 - d4.1)/10e-6

    cbind(d1, d2, d3, d4)
    }

text <- "Broken-stick model (Alvarado and Bradford, 2002) - Reparameterised"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text, deriv1 = deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}

#From Bradford 2002 - Broken-Stick Model - Original
GRT.BS2.fun <- function(Temp, k, Tb, To, ThetaT){
  t2 <- ifelse(Temp < Tb, Tb, ifelse(Temp > To, To, Temp))
  t1 <- ifelse(Temp < To, To, Temp)
  psival <- ifelse(1 - k*(t1 - To) > 0, 1 - k*(t1 - To), 0)
  GR <- psival * (t2 - Tb)/ThetaT
  GR }

"GRT.BS2" <- function(){
fct <- function(x, parm) {
  k <- parm[,1]; Tb <- parm[,2]; To <- parm[,3]; ThetaT <- parm[,4]
  GR <- GRT.BS2.fun(x, k, Tb, To, ThetaT)
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
deriv1 <- function(x, parms){

    #Approximation by using finite differences
    Temp <- x
    k <-  as.numeric(parms[,1]); Tb <- as.numeric(parms[,2]); To <- as.numeric(parms[,3])
    ThetaT <- as.numeric(parms[,4])

    d1.1 <- GRT.BS2.fun(Temp, k, Tb, To, ThetaT)
    d1.2 <- GRT.BS2.fun(Temp, (k + 10e-6), Tb, To, ThetaT)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- GRT.BS2.fun(Temp, k, Tb, To, ThetaT)
    d2.2 <- GRT.BS2.fun(Temp, k, (Tb + 10e-6), To, ThetaT)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- GRT.BS2.fun(Temp, k, Tb, To, ThetaT)
    d3.2 <- GRT.BS2.fun(Temp, k, Tb, (To + 10e-6), ThetaT)
    d3 <- (d3.2 - d3.1)/10e-6

    d4.1 <- GRT.BS2.fun(Temp, k, Tb, To, ThetaT)
    d4.2 <- GRT.BS2.fun(Temp, k, Tb, To, (ThetaT + 10e-6))
    d4 <- (d4.2 - d4.1)/10e-6

    cbind(d1, d2, d3, d4)
    }

text <- "Broken-stick model (Alvarado and Bradford, 2002) - Original"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text, deriv1 = deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}

