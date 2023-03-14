HTE1.fun <- function(time, Psi, G, Psib, sigmaPsib, thetaH, b){
  Pmax <- PmaxPsi1.fun(Psi, G, Psib, sigmaPsib)
  GR50 <- GRPsiLin.fun(Psi, Psib, thetaH)
  GR50 <- ifelse(GR50==0, 1e-08, GR50)
  plogis(b * (log(time) - log(1/GR50)))*Pmax
}

"HTE1"<- function(fixed = c(NA, NA, NA, NA, NA),
                  names = c("G", "Psib", "sigmaPsib", "thetaH", "b")){
# HT-to-event. With shifted exponential + GR50 linear

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

  S <- HTE1.fun(x[,1], x[,2], parmMat[,1], parmMat[,2], parmMat[,3],
                parmMat[,4], parmMat[,5])
  return(S)
}

names <- names[notFixed]

ss <- function(data){

  data <- subset(data, is.finite(data[,1])==T)
  PsiF <- factor(data[,2])

  result <- c()

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
  # result contiene i livelli da escludere dal fitting
  dataset_cum <- data # subset(data, is.finite(data[,1])==T)
  if(is.null(result) != T){
    for(i in 1:length(result)) dataset_cum <- subset(dataset_cum, dataset_cum[,2] != result[i])}
  PsiF <- factor(dataset_cum[,2])
  x1 <- dataset_cum[,1]
  x2 <- dataset_cum[,2]
  y <- dataset_cum[,3]

  modI <- drm(y ~ x1, fct=LL.3(), curveid=PsiF,
              pmodels=list(~1,~PsiF-1,~PsiF-1),
              data=dataset_cum)

  psiLevels <- as.numeric(levels(PsiF))
  b <- - coef(modI)[1]
  Pmax <- coef(modI)[2:(length(psiLevels)+1)]
  modPmax <- drm(Pmax ~ psiLevels, fct=PmaxPsi1(fixed = fixed[1:3]))

  if(is.na(fixed[1])){
    G <- coef(modPmax)[1]; Psib <- coef(modPmax)[2]; sigmaPsib <- coef(modPmax)[3]
  } else {
    G <- fixed[1]; Psib <- coef(modPmax)[1]; sigmaPsib <- coef(modPmax)[2]
  }

  GR50 <- 1/coef(modI)[(length(psiLevels)+2):length(coef(modI))]
  modGR <- drm(GR50 ~ psiLevels, fct=GRPsiLin())
  thetaH <- coef(modGR)[2]; Psib2 <- coef(modGR)[1]
  psib <- mean(Psib, Psib2)
  retVal <- c(G, psib, sigmaPsib, thetaH, b)

  return(retVal[is.na(fixed)]) }

GR <- function(parms, respl, reference="control", type="relative", Psi){
  # Questa funzione restituisce il germination rate, not time
  # respl è su una scala relativa ]0,1[
  HTE1.gra <- function(G, Psib, sigmaPsib, thetaH, b, Psi, g) {
    .Pmax <- PmaxPsi1.fun(Psi, G, Psib, sigmaPsib)
    .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
    .temp2 <- (.Pmax - g)/g
    .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
    .GR50 <- GRPsiLin.fun(Psi, Psib, thetaH)
    .GR50 <- ifelse(.GR50>0, .GR50, 0)
    res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
    res
  }
  HTE1.graRel <- function(G, Psib, sigmaPsib, thetaH, b, Psi, g) {
    .temp2 <- (1 - g)/g
    .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
    .GR50 <- GRPsiLin.fun(Psi, Psib, thetaH)
    .GR50 <- ifelse(.GR50>0, .GR50, 0)
    res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
    res
  }
   G <- as.numeric(parms[1]); Psib <- as.numeric(parms[2])
   sigmaPsib <- as.numeric(parms[3]); thetaH <- as.numeric(parms[4])
   b <- as.numeric(parms[5])
   g <- respl #/100
  if(type=="absolute"){
    EDp <- HTE1.gra(G, Psib, sigmaPsib, thetaH, b, Psi, g)

    #Approximation of derivatives(finite differences)
    d1.1 <- HTE1.gra(G, Psib, sigmaPsib, thetaH, b, Psi, g)
    d1.2 <- HTE1.gra(G + 10e-6, Psib, sigmaPsib, thetaH, b, Psi, g)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HTE1.gra(G, Psib, sigmaPsib, thetaH, b, Psi, g)
    d2.2 <- HTE1.gra(G, Psib + 10e-6, sigmaPsib, thetaH, b, Psi, g)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HTE1.gra(G, Psib, sigmaPsib, thetaH, b, Psi, g)
    d3.2 <- HTE1.gra(G, Psib, sigmaPsib + 10e-6, thetaH, b, Psi, g)
    d3<- (d3.2 - d3.1)/10e-6

    d4.1 <- HTE1.gra(G, Psib, sigmaPsib, thetaH, b, Psi, g)
    d4.2 <- HTE1.gra(G, Psib, sigmaPsib, thetaH+ 10e-6, b, Psi, g)
    d4 <- (d4.2 - d4.1)/10e-6

    d5.1 <- HTE1.gra(G, Psib, sigmaPsib, thetaH, b, Psi, g)
    d5.2 <- HTE1.gra(G, Psib, sigmaPsib, thetaH, b + 10e-6, Psi, g)
    d5 <- (d5.2 - d5.1)/10e-6

    EDder <- c(d1, d2, d3, d4, d5)

  } else{ if(type == "relative") {
    EDp <- HTE1.graRel(G, Psib, sigmaPsib, thetaH, b, Psi, g)

    #Approximation of derivatives(finite differences)
    d1.1 <- HTE1.graRel(G, Psib, sigmaPsib, thetaH, b, Psi, g)
    d1.2 <- HTE1.graRel(G + 10e-6, Psib, sigmaPsib, thetaH, b, Psi, g)
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- HTE1.graRel(G, Psib, sigmaPsib, thetaH, b, Psi, g)
    d2.2 <- HTE1.graRel(G, Psib + 10e-6, sigmaPsib, thetaH, b, Psi, g)
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- HTE1.graRel(G, Psib, sigmaPsib, thetaH, b, Psi, g)
    d3.2 <- HTE1.graRel(G, Psib, sigmaPsib + 10e-6, thetaH, b, Psi, g)
    d3<- (d3.2 - d3.1)/10e-6

    d4.1 <- HTE1.graRel(G, Psib, sigmaPsib, thetaH, b, Psi, g)
    d4.2 <- HTE1.graRel(G, Psib, sigmaPsib, thetaH+ 10e-6, b, Psi, g)
    d4 <- (d4.2 - d4.1)/10e-6

    d5.1 <- HTE1.graRel(G, Psib, sigmaPsib, thetaH, b, Psi, g)
    d5.2 <- HTE1.graRel(G, Psib, sigmaPsib, thetaH, b + 10e-6, Psi, g)
    d5 <- (d5.2 - d5.1)/10e-6

    EDder <- c(d1, d2, d3, d4, d5)

  } }
return(list(EDp, EDder))
}

deriv1 <- function(x, parm){
  #Approximation by using finite differences
  parmMat <- matrix(parmVec, nrow(parm),
                        numParm, byrow = TRUE)
  parmMat[, notFixed] <- parm

  a <- as.numeric(parmMat[,1])
  b <- as.numeric(parmMat[,2])
  ci <- as.numeric(parmMat[,3])
  di <- as.numeric(parmMat[,4])
  e <- as.numeric(parmMat[,5])

  d1.1 <- HTE1.fun(x[,1], x[,2], a, b, ci,
                   di, e)
  d1.2 <- HTE1.fun(x[,1], x[,2], (a + 10e-6), b, ci,
                   di, e)
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTE1.fun(x[,1], x[,2], a, b, ci,
                   di, e)
  d2.2 <- HTE1.fun(x[,1], x[,2], a, (b + 10e-6), ci,
                   di, e)
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTE1.fun(x[,1], x[,2], a, b, ci,
                   di, e)
  d3.2 <- HTE1.fun(x[,1], x[,2], a, b, (ci + 10e-6),
                   di, e)
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- HTE1.fun(x[,1], x[,2], a, b, ci,
                   di, e)
  d4.2 <- HTE1.fun(x[,1], x[,2], a, b, ci,
                   (di + 10e-6), e)
  d4 <- (d4.2 - d4.1)/10e-6

  d5.1 <- HTE1.fun(x[,1], x[,2], a, b, ci,
                   di, e)
  d5.2 <- HTE1.fun(x[,1], x[,2], a, b, ci,
                   di, (e + 10e-6))
  d5 <- (d5.2 - d5.1)/10e-6

  cbind(d1, d2, d3, d4, d5)[,notFixed]
}
name <- "HTE1"
text <- "Hydro-time model with shifted exponential for Pmax and linear model for GR50"
returnList <- list(fct=fct, ssfct=ss, name = name, names=names, text=text, edfct=GR, deriv1=deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}

