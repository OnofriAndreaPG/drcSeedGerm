#From Mesgaran et al., 2017
GRTPsi.M.fun <- function(Temp, Psi, k, Tb, ThetaHT, Psib) {
  t2 <- ifelse(Temp < Tb, Tb, Temp)
  psival <- ifelse(Psi - Psib - k * (Temp - Tb) > 0,
                   Psi - Psib - k * (Temp - Tb), 0)
  GR <- psival * (t2 - Tb)/ThetaHT
  return(ifelse(GR < 0 , 0 , GR)) }

"GRTPsi.M" <- function(){

#Mean function
fct <- function(x, parm) {
  GR <- GRTPsi.M.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[,4])
  return(GR)  }
#Names
names <- c("k", "Tb", "ThetaHT", "Psib")

#Self starter
ss <- function(data){
  # x <- data[, 1]; y <- data[, 2]
  # coefs <- coef( lm(y ~ x + I(x^2)))
  # a <- coefs[3]; b <- coefs[2]; c <- coefs[1]
  # Tb <- (-b + sqrt(b^2 - 4*a*c))/(2 * a)
  # Tc <- (-b - sqrt(b^2 - 4*a*c))/(2 * a)
  # ThetaT <- 1/b
  k <- 0.008
  Tb <- 8
  ThetaHT <- 200
  Psib <- -1
  return(c(k, Tb, ThetaHT, Psib))}

## Defining derivatives
deriv1 <- function(x, parms){

    #Approximation by using finite differences
    Temp <- x[,1]
    Psi <- x[,2]
    k <-  as.numeric(parms[,1]); Tb <- as.numeric(parms[,2])
    ThetaHT <- as.numeric(parms[,3]); Psib <- as.numeric(parms[,4])

    d1.1 <- GRTPsi.M.fun(Temp, Psi, k, Tb, ThetaHT, Psib)
    d1.2 <- GRTPsi.M.fun(Temp, Psi, (k + 10e-7), Tb, ThetaHT, Psib)
    d1 <- (d1.2 - d1.1)/10e-7

    d2.1 <- GRTPsi.M.fun(Temp, Psi, k, Tb, ThetaHT, Psib)
    d2.2 <- GRTPsi.M.fun(Temp, Psi, k, (Tb + 10e-7), ThetaHT, Psib)
    d2 <- (d2.2 - d2.1)/10e-7

    d3.1 <- GRTPsi.M.fun(Temp, Psi, k, Tb, ThetaHT, Psib)
    d3.2 <- GRTPsi.M.fun(Temp, Psi, k, Tb, (ThetaHT + 10e-7), Psib)
    d3 <- (d3.2 - d3.1)/10e-7

    d4.1 <- GRTPsi.M.fun(Temp, Psi, k, Tb, ThetaHT, Psib)
    d4.2 <- GRTPsi.M.fun(Temp, Psi, k, Tb, ThetaHT, (Psib + 10e-7))
    d4 <- (d4.2 - d4.1)/10e-7

    cbind(d1, d2, d3, d4)
    }
## Defining descriptive text
text <- "HTT model with polynomial temperature effect (Mohsen et al., 2017)"
returnList <- list(fct = fct, ssfct = ss, names = names, text = text, deriv1 = deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
