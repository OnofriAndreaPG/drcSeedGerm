#Rowse and Finch-Savage, 2003 - Original formulation
GRT.RF.fun <- function(Temp, k, Tb, Td, ThetaT) {
  t2 <- ifelse(Temp < Tb, Tb, Temp)
  t1 <- ifelse(Temp < Td, Td, Temp)
  psival <- ifelse(1 - k*(t1 - Td) > 0, 1 - k*(t1 - Td), 0)
  GR <- psival * (t2 - Tb)/ThetaT
  return(ifelse(GR < 0 , 0 , GR)) }

"GRT.RF" <- function(){

## Defining the non-linear function    
fct <- function(x, parm) {
  GR <- GRT.RF.fun(x, parm[,1], parm[,2], parm[,3], parm[,4])
  return(ifelse(GR < 0 , 0 , GR))}

## Defining names
names <- c("k", "Tb", "Td", "ThetaT")

## Defining self starter function
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
  return(c(k, Tb, Td, ThetaT))}

## Defining derivatives
    deriv1 <- function(x, parms){
    
    #Approximation by using finite differences
    Temp <- x
    k <-  as.numeric(parms[,1]); Tb <- as.numeric(parms[,2]); Td <- as.numeric(parms[,3])
    ThetaT <- as.numeric(parms[,4])
  
    d1.1 <- GRT.RF.fun(Temp, k, Tb, Td, ThetaT)
    d1.2 <- GRT.RF.fun(Temp, (k + 10e-6), Tb, Td, ThetaT)
    d1 <- (d1.2 - d1.1)/10e-6
    
    d2.1 <- GRT.RF.fun(Temp, k, Tb, Td, ThetaT)
    d2.2 <- GRT.RF.fun(Temp, k, (Tb + 10e-6), Td, ThetaT)
    d2 <- (d2.2 - d2.1)/10e-6
    
    d3.1 <- GRT.RF.fun(Temp, k, Tb, Td, ThetaT)
    d3.2 <- GRT.RF.fun(Temp, k, Tb, (Td + 10e-6), ThetaT)
    d3 <- (d3.2 - d3.1)/10e-6
    
    d4.1 <- GRT.RF.fun(Temp, k, Tb, Td, ThetaT)
    d4.2 <- GRT.RF.fun(Temp, k, Tb, Td, (ThetaT + 10e-6))
    d4 <- (d4.2 - d4.1)/10e-6
  
    cbind(d1, d2, d3, d4)
    }
## Defining descriptive text    
text <- "Rowse - Finch-Savage model (Rowse & Finch-Savage, 2003)"

## Returning the function with self starter and names
returnList <- list(fct = fct, ssfct = ss, names = names, text = text, deriv1 = deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
