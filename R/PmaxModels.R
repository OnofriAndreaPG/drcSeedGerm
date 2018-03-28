#Pmax vs Psi (shifted exponential, with asymptote)
PmaxPsi1.fun <- function(Psi, G, Psib, sigma) {
  Pmax <- G * (1 - exp( - (Psi - Psib)/sigma))
  Pmax <- ifelse(Pmax < 0 , 0, Pmax)
  return(Pmax) }

"PmaxPsi1" <- function(){
Pmax1.fct <- function(x, parm) {
  Pmax <- PmaxPsi1.fun(x, parm[,3], parm[,1], parm[,2])
  return(Pmax) }
Pmax1.names <- c("Psib", "sigma", "G")
Pmax1.text <- "Shifted exponential distribution of base osmotic potential"
Pmax1.ss <- function(data){
  data <- subset(data, data[,2]!=0)
  x <- data[, 1]
  y <- data[, 2]
  G <- max(y) * 1.05
  pseudoY <- -log( (G - y)/G)
  pseudoX <- x
  coefs <- coef( lm(pseudoY ~ pseudoX) )
  a <- coefs[1]
  b <- coefs[2]
  sigma <- 1/b
  Psib <- - a * sigma
  return(c(Psib, sigma, G))}
Pmax1 <- list(fct = Pmax1.fct, ssfct = Pmax1.ss, names = Pmax1.names, text = Pmax1.text)
class(Pmax1) <- "drcMean"
invisible(Pmax1)
}

PmaxT1.fun <- function(Temp, G, Tc, sigmaTc) {
  Pmax <- G * (1 - exp( - (Tc - Temp) / sigmaTc))
  Pmax <- ifelse(Pmax < 0, 0, Pmax)
  return(Pmax)}
"PmaxT1" <- function(){
fct <- function(x, parm) {
  Pmax <- PmaxT1.fun(x, parm[,1], parm[,2], parm[,3])
  return(Pmax)
}
names <- c("G", "Tc", "sigmaTc")
ss <- function(data){
  x <- data[,1]; y <- data[,2]
  y[y == 0] <- 10e-6
  G <- max(y) * 1.00005
  pseudoY <- log( (G - y)/G )
  coefs <- coef( lm(pseudoY ~ x) )
  sigmaTc <- 1/coefs[2]
  Tc <- - coefs[1] * sigmaTc
  #G=0.9; Tc=36; sigmaTc=1.5
return(c(G, Tc, sigmaTc))}
text <- "Truncated shifted exponential distribution for cut-off temperature"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}
