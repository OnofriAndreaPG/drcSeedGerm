#Hydrotime models ######################################
HTex.fun <- function(time, Psi, ThetaH, mu, sigma){
  1 - exp(-exp( (Psi - (ThetaH/time) - mu)/sigma) ) }
"HTex" <- function(){
fct <- function(x, parm){
  HTex.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3]) }
names <- c("ThetaH", "mu", "sigma")
name <- "HTex"
text <- "Hydrotime model with Type II extreme distribution of Psib"
ss <- function(data){
  x1 <- data[, 1]
  x2 <- data[, 2]
  y <- data[, 3]
  pseudoY <- qnorm((y+10e-6)*0.99)
  mod <- lm(pseudoY ~ I(1/x1) + x2)
  sigmaPsib <- 1/coef(mod)[3]
  Psib50 <- -coef(mod)[1]*sigmaPsib
  ThetaH <- -coef(mod)[2]*sigmaPsib
  return(c(ThetaH, Psib50, sigmaPsib))
}
returnList <- list(fct=fct, ssfct=ss, name = name, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}
