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

deriv1 <- function(x, parms){

    #Approximation by using finite differences
    Temp <- x
    k <-  as.numeric(parms[,1]); Tb <- as.numeric(parms[,2]); Tc <- as.numeric(parms[,3])
    ThetaT <- as.numeric(parms[,4])

    d1.1 <- GRT.YL.fun(Temp, k, Tb, Tc, ThetaT)
    d1.2 <- GRT.YL.fun(Temp, (k + 10e-6), Tb, Tc, ThetaT)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- GRT.YL.fun(Temp, k, Tb, Tc, ThetaT)
    d2.2 <- GRT.YL.fun(Temp, k, (Tb + 10e-6), Tc, ThetaT)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- GRT.YL.fun(Temp, k, Tb, Tc, ThetaT)
    d3.2 <- GRT.YL.fun(Temp, k, Tb, (Tc + 10e-6), ThetaT)
    d3 <- (d3.2 - d3.1)/10e-6

    d4.1 <- GRT.YL.fun(Temp, k, Tb, Tc, ThetaT)
    d4.2 <- GRT.YL.fun(Temp, k, Tb, Tc, (ThetaT + 10e-6))
    d4 <- (d4.2 - d4.1)/10e-6

    cbind(d1, d2, d3, d4)
    }

text <- "Yield-loss function (Kropff and van Laar, 1993)"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text, deriv1 = deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
