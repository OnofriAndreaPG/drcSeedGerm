# Pmax vs Psi (shifted exponential, with asymptote)
PmaxPsi1.fun <- function(Psi, G, Psib, sigma) {
  Pmax <- G * (1 - exp( - (Psi - Psib)/sigma))
  Pmax <- ifelse(Pmax < 0 , 0, Pmax)
  return(Pmax) }

"PmaxPsi1" <- function(fixed = c(NA, NA, NA),
                       names = c("G", "Psib", "sigma"))
  {

  ## Checking arguments
  numParm <- 3
  if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
  if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

  ## Fixing parameters (using argument 'fixed')
  notFixed <- is.na(fixed)
  parmVec <- rep(0, numParm)
  parmVec[!notFixed] <- fixed[!notFixed]

  Pmax1.fct <- function(x, parm) {
    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
    parmMat[, notFixed] <- parm
    Pmax <- PmaxPsi1.fun(x, parmMat[,1], parmMat[,2], parmMat[,3])
    return(Pmax)
    }

  Pmax1.names <- names[notFixed]
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
  return(c(G, Psib, sigma)[notFixed])}

deriv1 <- function(x, parm){
  #Approximation by using finite differences
  # derivate parziali sui parametri
  parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
  parmMat[, notFixed] <- parm

  d1.1 <- PmaxPsi1.fun(x, parmMat[,1], parmMat[,2], parmMat[,3])
  d1.2 <- PmaxPsi1.fun(x, (parmMat[,1] + 10e-6), parmMat[,2], parmMat[,3])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- PmaxPsi1.fun(x, parmMat[,1], parmMat[,2], parmMat[,3])
  d2.2 <- PmaxPsi1.fun(x, parmMat[,1], (parmMat[,2] + 10e-6), parmMat[,3])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- PmaxPsi1.fun(x, parmMat[,1], parmMat[,2], parmMat[,3])
  d3.2 <- PmaxPsi1.fun(x, parmMat[,1], parmMat[,2], (parmMat[,3] + 10e-6))
  d3 <- (d3.2 - d3.1)/10e-6

  cbind(d1, d2, d3)[, notFixed]
}

derivx <- function(x, parm){
  d1.1 <- PmaxPsi1.fun(x, parm[,1], parm[,2], parm[,3])
  d1.2 <- PmaxPsi1.fun(x + 10e-6, parm[,1], parm[,2],
                      parm[,3])
  d1 <- (d1.2 - d1.1)/10e-6
  d1
}

Pmax1 <- list(fct = Pmax1.fct, ssfct = Pmax1.ss, names = Pmax1.names,
              text = Pmax1.text, deriv1 = deriv1, derivx = derivx)
class(Pmax1) <- "drcMean"
invisible(Pmax1)
}


PmaxT1.fun <- function(Temp, G, Tc, sigmaTc) {
  Pmax <- G * (1 - exp( - (Tc - Temp) / sigmaTc))
  Pmax <- ifelse(Pmax < 0, 0, Pmax)
  return(Pmax)}

"PmaxT1" <- function(fixed = c(NA, NA, NA),
                       names = c("G", "Tc", "sigmaTc")){

    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    fct <- function(x, parm) {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      Pmax <- PmaxT1.fun(x, parmMat[,1], parmMat[,2], parmMat[,3])
      return(Pmax)
    }
    names <- c("G", "Tc", "sigmaTc")[notFixed]
    ss <- function(data){
      x <- data[,1]; y <- data[,2]
      y[y == 0] <- 10e-6
      G <- max(y) * 1.00005
      pseudoY <- log( (G - y)/G )
      coefs <- coef( lm(pseudoY ~ x) )
      sigmaTc <- 1/coefs[2]
      Tc <- - coefs[1] * sigmaTc
      #G=0.9; Tc=36; sigmaTc=1.5
    return(c(G, Tc, sigmaTc)[notFixed])}

    deriv1 <- function(x, parm){
      #Approximation by using finite differences
      # derivate parziali sui parametri
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm

      d1.1 <- PmaxT1.fun(x, parmMat[,1], parmMat[,2], parmMat[,3])
      d1.2 <- PmaxT1.fun(x, (parmMat[,1] + 10e-6), parmMat[,2], parmMat[,3])
      d1 <- (d1.2 - d1.1)/10e-6

      d2.1 <- PmaxT1.fun(x, parmMat[,1], parmMat[,2], parmMat[,3])
      d2.2 <- PmaxT1.fun(x, parmMat[,1], (parmMat[,2] + 10e-6), parmMat[,3])
      d2 <- (d2.2 - d2.1)/10e-6

      d3.1 <- PmaxT1.fun(x, parmMat[,1], parmMat[,2], parmMat[,3])
      d3.2 <- PmaxT1.fun(x, parmMat[,1], parmMat[,2], (parmMat[,3] + 10e-6))
      d3 <- (d3.2 - d3.1)/10e-6

      cbind(d1, d2, d3)[, notFixed]
    }

    derivx <- function(x, parm){
      d1.1 <- PmaxT1.fun(x, parm[,1], parm[,2], parm[,3])
      d1.2 <- PmaxT1.fun(x + 10e-6, (parm[,1]), parm[,2],
                          parm[,3])
      d1 <- (d1.2 - d1.1)/10e-6
      d1
    }

    text <- "Truncated shifted exponential distribution for cut-off temperature"
    returnList <- list(fct=fct, ssfct=ss, names=names, text=text,
                   deriv1 = deriv1, derivx = derivx)
    class(returnList) <- "drcMean"
    invisible(returnList)
}
