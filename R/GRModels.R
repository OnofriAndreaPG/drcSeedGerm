#Models of GR vs water potential ####################

GRPsiLin.fun <- function(Psi, Psib, thetaH) {
  #Linear hydrotime model (Bradford, 2002)
  GR50 <- ( Psi - Psib ) / thetaH
  return(ifelse(GR50 < 0, 0, GR50)) }

"GRPsiLin" <- function() {
  fct <- function(x, parm)
  { GR50 <- GRPsiLin.fun(x, parm[,1], parm[,2])
    return(ifelse(GR50<0, 0, GR50)) }
  ssfct <- function(data)
  { x <- data[, 1]
  y <- data[, 2]
  mod <- lm( y ~ x )
  thetaH <- 1/coef(mod)[2]
  Psib <- - coef(mod)[2]*thetaH
  return(c(Psib, thetaH))
  }
  names <- c("Psib", "thetaH")
  text <- "Linear hydrotime model (Bradford, 2002)"
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

# second order polynomial
GRPsiPol.fun <- function(Psi, Psib, thetaH) {
  #"Polynomial hydrotime model - Convex up"
  GR50 <- - ( Psi^2 - Psib^2) / (thetaH)
  return(ifelse(GR50 < 0, 0, GR50)) }

"GRPsiPol" <- function(){
   GR502.fct <- function(x, parm) {
    Psi <- x; Psib <- parm[,1]; thetaH <- parm[,2]
    GR50 <- GRPsiPol.fun(Psi, Psib, thetaH)
    return(ifelse(GR50 < 0, 0, GR50)) }
  GR502.names <- c("Psib", "thetaH")
  GR502.ss <- function(data){
    x <- data[, 1]
    y <- data[, 2]
    isPositive <- y > 0
    x1 <- x[isPositive]
    y1 <- y[isPositive]
    mod <- lm(y1 ~ I(x1^2))
    thetaH <- - 1/coef(mod)[2]
    Psib <- - sqrt(thetaH * coef(mod)[1])
    return(c(Psib,thetaH))}
  GR502.text <- "Polynomial hydrotime model - Convex up"
  GR502 <- list(fct = GR502.fct, ssfct=GR502.ss, names=GR502.names, text = GR502.text)
  class(GR502) <- "drcMean"
  invisible(GR502)
}

# second order polynomial - 2
GRPsiPol2.fun <- function(Psi, Psib, thetaH) {
  #Polynomial hydrotime model - Convex down
  GR50 <- (( Psi - Psib)^2)/ (thetaH)
  return(ifelse(GR50 < 0, 0, GR50)) }
"GRPsiPol2" <- function(){
  fct <- function(x, parm) {
    Psi <- x; Psib <- parm[,1]; thetaH <- parm[,2]
    GR50 <- GRPsiPol2.fun(Psi, Psib, thetaH)
    return(ifelse(Psi < Psib, 0, GR50)) }
  names <- c("Psib", "thetaH")
  ss <- function(data){
    x <- data[, 1]
    y <- data[, 2]
    isPositive <- y > 0
    x1 <- x[isPositive]
    y1 <- y[isPositive]
    Psib <- x1[which( y1==min(y1) )]
    mod <- lm(y1 ~ I((x1 - Psib)^2) - 1)
    thetaH <- 1/coef(mod)[1]
    return(c(Psib,thetaH))}
  text <- "Polynomial hydrotime model - Convex down"
  returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}
