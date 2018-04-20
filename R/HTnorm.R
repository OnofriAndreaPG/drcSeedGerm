#Hydrotime model with normal distribution of base water potential
HTnorm.fun <- function(time, Psi, ThetaH, Psib50, sigmaPsib){
  pnorm((Psi - (ThetaH/time) - Psib50)/sigmaPsib) }
"HTnorm" <- function(){
fct <- function(x, parm){
  HTnorm.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3]) }
names <- c("ThetaH", "Psib50", "sigmaPsib")
text <- "Hydrotime model with normal distribution of Psib (Bradford et al., 2002)"
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
returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}
