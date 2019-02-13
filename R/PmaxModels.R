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

# Psi <- seq(-1.5, 0, 0.1)
# Pmax <- PmaxPsi1na.fun(Psi, -0.7, 0.2)
# Pmax <-Pmax + rnorm(16, 0, 0.05)
# Pmax <- ifelse(Pmax < 0, 0, Pmax)
# Pmax <- ifelse(Pmax > 1, 0.99, Pmax)
# Pmax
#
# mod <- drm(Pmax ~ Psi, fct = PmaxPsi1na())
# summary(mod)
#
# mod <- drm(Pmax ~ Psi, fct = PmaxPsi1())
# summary(mod)

#Pmax vs Psi (shifted exponential, with asymptote)
PmaxPsi1na.fun <- function(Psi, Psib, sigma) {
  Pmax <- 1 - exp( - (Psi - Psib)/sigma)
  Pmax <- ifelse(Pmax < 0 , 0, Pmax)
  return(Pmax) }

"PmaxPsi1na" <- function(){
Pmax1na.fct <- function(x, parm) {
  Pmax <- PmaxPsi1na.fun(x, parm[,1], parm[,2])
  return(Pmax) }
Pmax1na.names <- c("Psib", "sigma")
Pmax1na.text <- "Shifted exponential distribution of base osmotic potential"
Pmax1na.ss <- function(data){
  data <- subset(data, data[,2]!=0)
  x <- data[, 1]
  y <- data[, 2]
  G <- 1.05
  pseudoY <- -log( (G - y)/G)
  pseudoX <- x
  coefs <- coef( lm(pseudoY ~ pseudoX) )
  a <- coefs[1]
  b <- coefs[2]
  sigma <- 1/b
  Psib <- - a * sigma
  return(c(Psib, sigma))}
Pmax1na <- list(fct = Pmax1na.fct, ssfct = Pmax1na.ss, names = Pmax1na.names, text = Pmax1na.text)
class(Pmax1na) <- "drcMean"
invisible(Pmax1na)
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
