GRT.GH2.fun <- function(Temp, Tb, beta){
  #Linear model for suboptimal temperatures
  t2 <- ifelse(Temp < Tb, Tb, Temp)
  GR <- (t2 - Tb) * beta
  GR }
"GRT.GH2" <- function() {
  fct <- function(x, parm) {
    Tb <- parm[,1]; beta <- parm[,2]
    GR <- GRT.GH2.fun(x, Tb, beta)
    return(GR) }
  names <- c("Tb", "beta")
  ssfct <- function(data){
    x <- data[, 1]
    y <- data[, 2]
    ss1 <- coef( lm(y ~ x) )
    beta <- ss1[2]
    Tb <- - ss1[1] / beta
    return(c(Tb, beta))}
deriv1 <- function(x, parm){
  #Approximation by using finite differences
  d1.1 <- GRT.GH2.fun(x, parm[,1], parm[,2])
  d1.2 <- GRT.GH2.fun(x, (parm[,1] + 10e-6), parm[,2])
  d1 <- (d1.2 - d1.1)/10e-6
  d2.1 <- GRT.GH2.fun(x, parm[,1], parm[,2])
  d2.2 <- GRT.GH2.fun(x, parm[,1], (parm[,2] + 10e-6))
  d2 <- (d2.2 - d2.1)/10e-6
  cbind(d1, d2)
}

  text <- "Linear model for suboptimal temperatures"
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text, deriv1 = deriv1)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

