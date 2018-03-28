#Hydrotime models ######################################
HTL.fun <- function(time, Psi, ThetaH, Psib50, sigmaPsib){
  plogis((Psi - (ThetaH/time) - Psib50)/sigmaPsib) }
"HTL" <- function(){
fct <- function(x, parm){
  HTL.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3]) }
names <- c("ThetaH", "Psib50", "sigmaPsib")
text <- "Hydrotime model with logistic distribution of Psib (Mesgaran et al., 2013)"
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
