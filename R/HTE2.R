# HT-to-event. With shifted exponential + GR50polynomial (convex down)
HTE2.fun <- function(time, Psi, G, Psib, sigmaPsib, thetaH, b){
  Pmax <- PmaxPsi1.fun(Psi, G, Psib, sigmaPsib)
  GR50 <- GRPsiPol2.fun(Psi, Psib, thetaH)
  plogis(b * (log(time) - log(1/GR50)))*Pmax
}
"HTE2" <- function(fixed = c(NA, NA, NA, NA, NA),
                   names = c("G", "Psib", "sigmaPsib", "thetaH", "b")){
  ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if ( !(length(fixed) == numParm) ) {stop("Not correct 'fixed' argument")}
    # Only G can be constrained
    if (any(!is.na(fixed[2:5]))) {stop("Only the G parameter can be constrained, at the moment")}

    ## Handling 'fixed' argument
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

fct <- function(x, parm){
  parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
  parmMat[, notFixed] <- parm
  parm <- parmMat

  S <- HTE2.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                parm[,4], parm[,5])
  return(S)
}
names <- names[notFixed]
name <- "HTE2"

ss <- function(data){
  data <- subset(data, is.finite(data[,1])==T)
  result <- c()
  PsiF <- factor(data[,2])
  for(i in 1:length(levels(PsiF))){
    temp <- subset(data, data[,2] == levels(PsiF)[i])
    x <- temp[,1] + 0.0001; y <- temp[,3]

    # self-start
    d <- max(y) * 1.05
    pseudoY <- log((d - y)/(y + 0.00001)  + 1e-06 )
    coefs <- coef( lm(pseudoY ~ log(x+0.000001)))
    k <- -coefs[1]; b <- coefs[2]
    ED50 <- exp(k/b)
    modT <- try(nls(y ~ d/(1 + exp(b * (log(x + 0.000001) - log(ED50)))),
                    start = list(d = d, b = b, ED50 = ED50)), silent=T)
    if(inherits(modT, "try-error")) {
      res <- as.numeric(levels(PsiF)[i])
    result <- c(result, res)}
  }
  result
  dataset_cum <- subset(data, is.finite(data[,1])==T)

  if(is.null(result)!=T){
    for(i in 1:length(result)) dataset_cum <- subset(dataset_cum, dataset_cum[,2] != result[i])}
  PsiF <- factor(dataset_cum[,2])
  x1 <- dataset_cum[,1]
  x2 <- dataset_cum[,2]
  y <- dataset_cum[,3]

  modI <- drm(y ~ x1, fct=LL.3(), curveid=PsiF, pmodels=list(~1,~PsiF-1,~PsiF-1), data=dataset_cum)
  psiLevels <- as.numeric(levels(PsiF))
  b <- - coef(modI)[1]
  Pmax <- coef(modI)[2:(length(psiLevels)+1)]
  modPmax <- drm(Pmax ~ psiLevels, fct=PmaxPsi1())

  if(is.na(fixed[1])){
    G <- coef(modPmax)[1]; Psib <- coef(modPmax)[2]; sigmaPsib <- coef(modPmax)[3]
  } else {
    G <- fixed[1]; Psib <- coef(modPmax)[1]; sigmaPsib <- coef(modPmax)[2]
  }
  # G <- coef(modPmax)[3]; Psib <- coef(modPmax)[1]; sigmaPsib <- coef(modPmax)[2]

  GR50 <- 1/coef(modI)[(length(psiLevels)+2):length(coef(modI))]
  modGR <- drm(GR50 ~ psiLevels, fct=GRPsiPol2())
  thetaH <- coef(modGR)[2]; Psib2 <- coef(modGR)[1]
  psib <- mean(Psib, Psib2)
  return(c(G, Psib, sigmaPsib, thetaH, b)[is.na(fixed)]) }

GR <- function(parms, respl, reference="control", type="relative", Psi){
   G <- as.numeric(parms[1]); Psib<- as.numeric(parms[2])
   sigmaPsib<- as.numeric(parms[3]); thetaH<- as.numeric(parms[4])
   b <- as.numeric(parms[5])
   g <- respl/100
  if(type=="absolute"){
    .Pmax <- PmaxPsi1.fun(Psi, G, Psib, sigmaPsib)
    .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
    .temp2 <- (.Pmax - g)/g
    .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
    .GR50 <- GRPsiPol2.fun(Psi, Psib, thetaH)
    .GR50 <- ifelse(.GR50>0, .GR50, 0)
    res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
    EDp <- 1/res
    .GR <- EDp

    #Beginning of derivatives ###########################
    d1 <- exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
    g)) + log(thetaH/((Psi - Psib)^2))) * ((1/b) * ((1/g) * (1 -
    exp(-(Psi - Psib) * (1/sigmaPsib)))/((1/g) * (G * (1 - exp(-(Psi -
    Psib) * (1/sigmaPsib))) - g))))/exp(-(1/b) * log((1/g) *
    (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) - g)) + log(thetaH/((Psi -
    Psib)^2)))^2
    d2 <- -(exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
    g)) + log(thetaH/((Psi - Psib)^2))) * ((1/b) * ((1/g) * (G *
    (exp(-(Psi - Psib) * (1/sigmaPsib)) * (1/sigmaPsib)))/((1/g) *
    (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) - g))) + thetaH *
    (2 * (Psi - Psib))/((Psi - Psib)^2)^2/(thetaH/((Psi - Psib)^2)))/exp(-(1/b) *
    log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
        g)) + log(thetaH/((Psi - Psib)^2)))^2)
    d3 <- -(exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
    g)) + log(thetaH/((Psi - Psib)^2))) * ((1/b) * ((1/g) * (G *
    (exp(-(Psi - Psib) * (1/sigmaPsib)) * ((Psi - Psib) * (1/sigmaPsib^2))))/((1/g) *
    (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) - g))))/exp(-(1/b) *
    log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
        g)) + log(thetaH/((Psi - Psib)^2)))^2)
    d4 <- -(exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
    g)) + log(thetaH/((Psi - Psib)^2))) * (1/((Psi - Psib)^2)/(thetaH/((Psi -
    Psib)^2)))/exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi -
    Psib) * (1/sigmaPsib))) - g)) + log(thetaH/((Psi - Psib)^2)))^2)
    d5 <- -(exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
    g)) + log(thetaH/((Psi - Psib)^2))) * (1/b^2 * log((1/g) *
    (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) - g)))/exp(-(1/b) *
    log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
        g)) + log(thetaH/((Psi - Psib)^2)))^2)
    #End of derivatives #########################

    EDder <- c(d1,d2,d3,d4,d5)
  } else{ if(type=="relative") {
    .Pmax <- PmaxPsi1.fun(Psi, G, Psib, sigmaPsib)
    .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
    .temp2 <- (1 - g)/g
    .GR50 <- GRPsiPol2.fun(Psi, Psib, thetaH)
    .GR50 <- ifelse(.GR50>0, .GR50, 0)
    res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
    EDp <- 1/res
    .GR <- EDp
    d1 <- 0
    d2 <- -(exp(-(1/b) * log(((1 - g)/g)) + log(thetaH/((Psi - Psib)^2))) *
    (thetaH * (2 * (Psi - Psib))/((Psi - Psib)^2)^2/(thetaH/((Psi -
        Psib)^2)))/exp(-(1/b) * log(((1 - g)/g)) + log(thetaH/((Psi -
    Psib)^2)))^2)
    d3 <- 0
    d4 <- -(exp(-(1/b) * log(((1 - g)/g)) + log(thetaH/((Psi - Psib)^2))) *
    (1/((Psi - Psib)^2)/(thetaH/((Psi - Psib)^2)))/exp(-(1/b) *
    log(((1 - g)/g)) + log(thetaH/((Psi - Psib)^2)))^2)
    d5 <- -(exp(-(1/b) * log(((1 - g)/g)) + log(thetaH/((Psi - Psib)^2))) *
    (1/b^2 * log(((1 - g)/g)))/exp(-(1/b) * log(((1 - g)/g)) +
    log(thetaH/((Psi - Psib)^2)))^2)
    EDder <- c(d1,d2,d3,d4,d5)
  } }
return(list(EDp, EDder))
}
deriv1 <- function(x, parm){
  parmMat <- matrix(parmVec, nrow(parm),
                        numParm, byrow = TRUE)
  parmMat[, notFixed] <- parm
  parm <- parmMat
  #Approximation by using finite differences

  d1.1 <- HTE2.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d1.2 <- HTE2.fun(x[,1], x[,2], (parm[,1] + 10e-6), parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTE2.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d2.2 <- HTE2.fun(x[,1], x[,2], parm[,1], (parm[,2] + 10e-6), parm[,3],
                   parm[,4], parm[,5])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTE2.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d3.2 <- HTE2.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6),
                   parm[,4], parm[,5])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- HTE2.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d4.2 <- HTE2.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   (parm[,4] + 10e-6), parm[,5])
  d4 <- (d4.2 - d4.1)/10e-6

  d5.1 <- HTE2.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d5.2 <- HTE2.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], (parm[,5] + 10e-6))
  d5 <- (d5.2 - d5.1)/10e-6

  cbind(d1, d2, d3, d4, d5)[,notFixed]
}
text <- "Hydro-time model with shifted exponential for Pmax and polynomial model for GR50"
returnList <- list(fct=fct, ssfct=ss, name = name, names=names, text=text, edfct=GR, deriv1=deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
