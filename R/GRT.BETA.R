curve(GRT.BETA.fun(x, 6.5, 5, 20, 28), from=0, to=30)
curve(GRT.BETA.fun(x, 3.5, 5, 20, 28), add=T)
curve(GRT.BETA.fun(x, 1.5, 5, 20, 28), add=T)
curve(GRT.BETA.fun(x, 0.5, 5, 20, 28), add=T)

# Beta function
GRT.BETA.fun <- function(Temp, a, Tb, To, Tc) {
  GR <- ifelse(Temp > Tb & Temp < Tc, (((Temp - Tb)/(To - Tb)*(Tc - Temp)/(Tc - To))^((Tc - To)/(To - Tb)))^a, 0)
  return(ifelse(GR < 0 , 0 , GR)) }
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

deriv1 <- function(x, parm){
  #Approximation by using finite differences
  d1.1 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3],
                   parm[,4])
  d1.2 <- GRT.Exb.fun(x, (parm[,1] + 10e-6), parm[,2], parm[,3],
                   parm[,4])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3],
                   parm[,4])
  d2.2 <- GRT.Exb.fun(x, parm[,1], (parm[,2] + 10e-6), parm[,3],
                   parm[,4])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3],
                   parm[,4])
  d3.2 <- GRT.Exb.fun(x, parm[,1], parm[,2], (parm[,3] + 10e-6),
                   parm[,4])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3],
                   parm[,4])
  d4.2 <- GRT.Exb.fun(x, parm[,1], parm[,2], parm[,3],
                   (parm[,4] + 10e-6))
  d4 <- (d4.2 - d4.1)/10e-6

  cbind(d1, d2, d3, d4)
}

text <- "Exponential effect of temperature on GR50 (Type II - Masin et al., 2017)"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text, deriv1 = deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}

