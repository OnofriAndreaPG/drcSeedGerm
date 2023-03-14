# From Bradford 2002 - Broken-Stick HTT Model based on GR
# GRTPsi.BS.fun <- function(Temp, Psi, Tc, Tb, To, ThetaT){
#   t2 <- ifelse(Temp < Tb, Tb, ifelse(Temp > To, To, Temp))
#   t1 <- ifelse(Temp < To, To, Temp)
#   psival <- ifelse(1 - (t1 - To)/(Tc - To) > 0, 1 - (t1 - To)/(Tc - To), 0)
#   GR <- psival * (t2 - Tb)/ThetaT
#   GR }
#
# "GRT.BS" <- function(){
# fct <- function(x, parm) {
#   Tc <- parm[,1]; Tb <- parm[,2]; To <- parm[,3]; ThetaT <- parm[,4]
#   GR <- GRT.BS.fun(x, Tc, Tb, To, ThetaT)
#   return(ifelse(GR < 0 , 0 , GR)) }
# names <- c("Tc", "Tb", "To", "ThetaT")
# ss <- function(data){
#   pos <- which( data[,2]==max(data[,2]) )
#   len <- length( data[,2] )
#
#   reg1 <- data[1:pos, ]
#   reg2 <- data[pos:len, ]
#   x1 <- reg1[,1]; y1 <- reg1[, 2]
#   x2 <- reg2[,1]; y2 <- reg2[, 2]
#
#   ss1 <- coef( lm(y1 ~ x1) )
#   ThetaT <- 1/ss1[2]
#   Tb <- - ss1[1] * ThetaT
#   ss2 <- coef( lm((1-y2) ~ x2) )
#   k <- ss2[2]
#   To <- - ss2[1] / k
#   Tc <- To + 1/k
#
#   #k <- 0.1; Tb <- 2; To <- 20; ThetaT <- 35
#   return(c(Tc, Tb, To, ThetaT))}
# deriv1 <- function(x, parms){
#
#     #Approximation by using finite differences
#     Temp <- x
#     Tc <-  as.numeric(parms[,1]); Tb <- as.numeric(parms[,2]); To <- as.numeric(parms[,3])
#     ThetaT <- as.numeric(parms[,4])
#
#     d1.1 <- GRT.BS.fun(Temp, Tc, Tb, To, ThetaT)
#     d1.2 <- GRT.BS.fun(Temp, (Tc + 10e-6), Tb, To, ThetaT)
#     d1 <- (d1.2 - d1.1)/10e-6
#
#     d2.1 <- GRT.BS.fun(Temp, Tc, Tb, To, ThetaT)
#     d2.2 <- GRT.BS.fun(Temp, Tc, (Tb + 10e-6), To, ThetaT)
#     d2 <- (d2.2 - d2.1)/10e-6
#
#     d3.1 <- GRT.BS.fun(Temp, Tc, Tb, To, ThetaT)
#     d3.2 <- GRT.BS.fun(Temp, Tc, Tb, (To + 10e-6), ThetaT)
#     d3 <- (d3.2 - d3.1)/10e-6
#
#     d4.1 <- GRT.BS.fun(Temp, Tc, Tb, To, ThetaT)
#     d4.2 <- GRT.BS.fun(Temp, Tc, Tb, To, (ThetaT + 10e-6))
#     d4 <- (d4.2 - d4.1)/10e-6
#
#     cbind(d1, d2, d3, d4)
#     }
#
# text <- "Broken-stick model (Alvarado and Bradford, 2002) - Reparameterised"
# returnList <- list(fct=fct, ssfct=ss, names=names, text=text, deriv1 = deriv1)
# class(returnList) <- "drcMean"
# invisible(returnList)
# }

# From Bradford 2002 - Broken-Stick Model - Original
GRTPsi.BS.fun <- function(Temp, Psi, k, Tb, To, ThetaHT, Psib){
  t2 <- ifelse(Temp < Tb, Tb, ifelse(Temp > To, To, Temp))
  t1 <- ifelse(Temp < To, To, Temp)
  psival <- ifelse(Psi - (Psib + k*(t1 - To)) > 0,
                   Psi - (Psib + k*(t1 - To)), 0)
  GR <- psival * (t2 - Tb)/ThetaHT
  GR }

"GRTPsi.BS" <- function(){
fct <- function(x, parm) {
  k <- parm[,1]; Tb <- parm[,2]; To <- parm[,3];
  ThetaHT <- parm[,4]; Psib <- parm[,5]

  GR <- GRTPsi.BS.fun(x[,1], x[,2], k, Tb, To, ThetaHT, Psib)

  return(ifelse(GR < 0 , 0 , GR)) }
names <- c("k", "Tb", "To", "ThetaHT", "Psib")
ss <- function(data){
  k <- 0.1
  Tb <- 2
  To <- 20
  ThetaHT <- 350
  Psib <- -1
  return(c(k, Tb, To, ThetaHT, Psib))}
deriv1 <- function(x, parms){

    #Approximation by using finite differences
    Temp <- x[,1]
    Psi <- x[,2]
    k <-  as.numeric(parms[,1]); Tb <- as.numeric(parms[,2]);
    To <- as.numeric(parms[,3])
    ThetaHT <- as.numeric(parms[,4])
    Psib <- as.numeric(parms[,5])

    d1.1 <- GRTPsi.BS.fun(Temp, Psi, k, Tb, To, ThetaHT, Psib)
    d1.2 <- GRTPsi.BS.fun(Temp, Psi, (k + 10e-6), Tb, To, ThetaHT, Psib)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- GRTPsi.BS.fun(Temp, Psi, k, Tb, To, ThetaHT, Psib)
    d2.2 <- GRTPsi.BS.fun(Temp, Psi, k, (Tb + 10e-6), To, ThetaHT, Psib)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- GRTPsi.BS.fun(Temp, Psi, k, Tb, To, ThetaHT, Psib)
    d3.2 <- GRTPsi.BS.fun(Temp, Psi, k, Tb, (To + 10e-6), ThetaHT, Psib)
    d3 <- (d3.2 - d3.1)/10e-6

    d4.1 <- GRTPsi.BS.fun(Temp, Psi, k, Tb, To, ThetaHT, Psib)
    d4.2 <- GRTPsi.BS.fun(Temp, Psi, k, Tb, To, (ThetaHT + 10e-6), Psib)
    d4 <- (d4.2 - d4.1)/10e-6

    d5.1 <- GRTPsi.BS.fun(Temp, Psi, k, Tb, To, ThetaHT, Psib)
    d5.2 <- GRTPsi.BS.fun(Temp, Psi, k, Tb, To, ThetaHT, (Psib + 10e-6))
    d5 <- (d5.2 - d5.1)/10e-6

    cbind(d1, d2, d3, d4, d5)
    }

text <- "Broken-stick HTT model (Alvarado and Bradford, 2002)"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text) #,
                   # deriv1 = deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}

