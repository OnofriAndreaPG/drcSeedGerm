#From Mesgaran et al., 2017
GRT.M.fun <- function(predictor, Tc, Tb, ThetaT) {
  Temp <- predictor
  t2 <- ifelse(Temp < Tb, Tb, Temp)
  psival <- ifelse(1 - (Temp - Tb)/(Tc - Tb) > 0, 1 - (Temp - Tb)/(Tc - Tb), 0)
  GR <- psival * (t2 - Tb)/ThetaT
  return(ifelse(GR < 0 , 0 , GR)) }

"GRT.M" <- function(){

#Mean function  
fct <- function(x, parm) {
  GR <- GRT.M.fun(x, parm[,1], parm[,2], parm[,3])
  return(GR)  }
#Names
names <- c("Tc", "Tb", "ThetaT")

#Self starter
ss <- function(data){
  x <- data[, 1]; y <- data[, 2]
  coefs <- coef( lm(y ~ x + I(x^2)))
  a <- coefs[3]; b <- coefs[2]; c <- coefs[1]
  Tb <- (-b + sqrt(b^2 - 4*a*c))/(2 * a)
  Tc <- (-b - sqrt(b^2 - 4*a*c))/(2 * a)
  ThetaT <- 1/b
  return(c(Tc, Tb, ThetaT))}

## Defining derivatives
deriv1 <- function(x, parms){
    
    #Approximation by using finite differences
    Temp <- x
    Tc <-  as.numeric(parms[,1]); Tb <- as.numeric(parms[,2])
    ThetaT <- as.numeric(parms[,3])
    
    d1.1 <- GRT.M.fun(Temp, Tc, Tb, ThetaT)
    d1.2 <- GRT.M.fun(Temp, (Tc + 10e-7), Tb, ThetaT)
    d1 <- (d1.2 - d1.1)/10e-7
    
    d2.1 <- GRT.M.fun(Temp, Tc, Tb, ThetaT)
    d2.2 <- GRT.M.fun(Temp, Tc, (Tb + 10e-7), ThetaT)
    d2 <- (d2.2 - d2.1)/10e-7
    
    d3.1 <- GRT.M.fun(Temp, Tc, Tb, ThetaT)
    d3.2 <- GRT.M.fun(Temp, Tc, Tb, (ThetaT + 10e-7))
    d3 <- (d3.2 - d3.1)/10e-7
  
    cbind(d1, d2, d3)
    }
## Defining descriptive text    
text <- "Polynomial temperature effect (Mohsen et al., 2017)"
returnList <- list(fct = fct, ssfct = ss, names = names, text = text, deriv1 = deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
