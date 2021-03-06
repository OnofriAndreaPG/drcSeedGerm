HTE1.fun <- function(time, Psi, G, Psib, sigmaPsib, thetaH, b){
  Pmax <- PmaxPsi1.fun(Psi, G, Psib, sigmaPsib)
  GR50 <- GRPsiLin.fun(Psi, Psib, thetaH)
  GR50 <- ifelse(GR50==0, 1e-08, GR50)
  plogis(b * (log(time) - log(1/GR50)))*Pmax
}

"HTE1"<- function(){
#HT-to-event. With shifted exponential + GR50 linear
fct <- function(x, parm){
  S <- HTE1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                parm[,4], parm[,5])
  return(S)
}
names <- c("G", "Psib", "sigmaPsib", "thetaH", "b")
ss <- function(data){
  data <- subset(data, is.finite(data[,1])==T)
  result <- c()
  PsiF <- factor(data[,2])
  for(i in 1:length(levels(PsiF))){
    temp <- subset(data, data[,2] == levels(PsiF)[i])
    #modT <- try(drm(temp[,3] ~ temp[,1], fct=LL.3()), silent=T)
    x <- temp[,1]; y <- temp[,3]
    modT <- try(nls(y ~ NLSLL.3(x, a, b, c)), silent=T)
    if(class(modT) == "try-error") {res <- as.numeric(levels(PsiF)[i])
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
  G <- coef(modPmax)[3]; Psib <- coef(modPmax)[1]; sigmaPsib <- coef(modPmax)[2]

  GR50 <- 1/coef(modI)[(length(psiLevels)+2):length(coef(modI))]
  modGR <- drm(GR50 ~ psiLevels, fct=GRPsiLin())
  thetaH <- coef(modGR)[2]; Psib2 <- coef(modGR)[1]
  psib <- mean(Psib, Psib2)
  return(c(G, psib, sigmaPsib, thetaH, b)) }
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
    .GR50 <- GRPsiLin.fun(Psi, Psib, thetaH)
    .GR50 <- ifelse(.GR50>0, .GR50, 0)
    res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
    EDp <- 1/res
    .GR <- EDp
      EDder <- c(
        1/.GR * ((1/b) * ((1/g) * (.Pmax/G)/((1/g) * (.Pmax - g))))/(1/.GR^2),
        -( ((1/b) * ((1/g) * (G * (exp(-(Psi - Psib) * (1/sigmaPsib)) * (1/sigmaPsib)))/((.Pmax - g)/g)) + 1/(Psi - Psib))/(1/.GR)),
        -(((1/b) * ((1/g) * (G * (exp(-(Psi - Psib) * (1/sigmaPsib)) * ((Psi - Psib) * (1/sigmaPsib^2))))/((.Pmax - g)/g)))/(1/.GR)),
        -(1/(Psi - Psib)/(1/.GR50))/(1/.GR),
        -(1/.GR * (1/b^2 * log((1/g) * (.Pmax - g)))/exp(-(1/b) * log((1/g) * (.Pmax - g)) + log(1/.GR50))^2)
      )
  } else{ if(type=="relative") {
    .Pmax <- PmaxPsi1.fun(Psi, G, Psib, sigmaPsib)
    .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
    .temp2 <- (1 - g)/g
    .GR50 <- GRPsiLin.fun(Psi, Psib, thetaH)
    .GR50 <- ifelse(.GR50>0, .GR50, 0)
    res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
    EDp <- 1/res
    .GR <- EDp
    EDder <- c(0,
        -(1/.GR * ((1/.GR50)/(Psi - Psib)/(1/.GR50))/(1/.GR50)^2),
        0,
        -(1/.GR * (1/(Psi - Psib)/(1/.GR50))/(1/.GR)^2),
        -(.GR * log((1 - g)/g) )/ (b^2)
      )
  } }
return(list(EDp, EDder))
}

deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HTE1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d1.2 <- HTE1.fun(x[,1], x[,2], (parm[,1] + 10e-6), parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTE1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d2.2 <- HTE1.fun(x[,1], x[,2], parm[,1], (parm[,2] + 10e-6), parm[,3],
                   parm[,4], parm[,5])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTE1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d3.2 <- HTE1.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6),
                   parm[,4], parm[,5])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- HTE1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d4.2 <- HTE1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   (parm[,4] + 10e-6), parm[,5])
  d4 <- (d4.2 - d4.1)/10e-6

  d5.1 <- HTE1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5])
  d5.2 <- HTE1.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], (parm[,5] + 10e-6))
  d5 <- (d5.2 - d5.1)/10e-6

  cbind(d1, d2, d3, d4, d5)
}

text <- "Hydro-time model with shifted exponential for Pmax and linear model for GR50"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text, edfct=GR, deriv1=deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}

HTE1na.fun <- function(time, Psi, Psib, sigmaPsib, thetaH, b){
  Pmax <- PmaxPsi1na.fun(Psi, Psib, sigmaPsib)
  GR50 <- GRPsiLin.fun(Psi, Psib, thetaH)
  GR50 <- ifelse(GR50==0, 1e-08, GR50)
  plogis(b * (log(time) - log(1/GR50)))*Pmax
}

"HTE1na"<- function(){
#HT-to-event. With shifted exponential + GR50 linear
fct <- function(x, parm){
  S <- HTE1na.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                parm[,4])
  return(S)
}
names <- c("Psib", "sigmaPsib", "thetaH", "b")
ss <- function(data){
  data <- subset(data, is.finite(data[,1])==T)
  result <- c()
  PsiF <- factor(data[,2])
  for(i in 1:length(levels(PsiF))){
    temp <- subset(data, data[,2] == levels(PsiF)[i])
    #modT <- try(drm(temp[,3] ~ temp[,1], fct=LL.3()), silent=T)
    x <- temp[,1]; y <- temp[,3]
    modT <- try(nls(y ~ NLSLL.3(x, a, b, c)), silent=T)
    if(class(modT) == "try-error") {res <- as.numeric(levels(PsiF)[i])
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
  modPmax <- drm(Pmax ~ psiLevels, fct=PmaxPsi1na())
  Psib <- coef(modPmax)[1]; sigmaPsib <- coef(modPmax)[2]

  GR50 <- 1/coef(modI)[(length(psiLevels)+2):length(coef(modI))]
  modGR <- drm(GR50 ~ psiLevels, fct=GRPsiLin())
  thetaH <- coef(modGR)[2]; Psib2 <- coef(modGR)[1]
  psib <- mean(Psib, Psib2)
  return(c(psib, sigmaPsib, thetaH, b)) }
GR <- function(parms, respl, reference="control", type="relative", Psi){
   Psib<- as.numeric(parms[1])
   sigmaPsib<- as.numeric(parms[2]); thetaH<- as.numeric(parms[3])
   b <- as.numeric(parms[4])
   g <- respl/100
  if(type=="absolute"){
     G <- 1
    .Pmax <- PmaxPsi1na.fun(Psi, Psib, sigmaPsib)
    .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
    .temp2 <- (.Pmax - g)/g
    .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
    .GR50 <- GRPsiLin.fun(Psi, Psib, thetaH)
    .GR50 <- ifelse(.GR50>0, .GR50, 0)
    res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
    EDp <- 1/res
    .GR <- EDp
      EDder <- c(
        -( ((1/b) * ((1/g) * (G * (exp(-(Psi - Psib) * (1/sigmaPsib)) * (1/sigmaPsib)))/((.Pmax - g)/g)) + 1/(Psi - Psib))/(1/.GR)),
        -(((1/b) * ((1/g) * (G * (exp(-(Psi - Psib) * (1/sigmaPsib)) * ((Psi - Psib) * (1/sigmaPsib^2))))/((.Pmax - g)/g)))/(1/.GR)),
        -(1/(Psi - Psib)/(1/.GR50))/(1/.GR),
        -(1/.GR * (1/b^2 * log((1/g) * (.Pmax - g)))/exp(-(1/b) * log((1/g) * (.Pmax - g)) + log(1/.GR50))^2)
      )
  } else{ if(type=="relative") {
    .Pmax <- PmaxPsi1na.fun(Psi, Psib, sigmaPsib)
    .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
    .temp2 <- (1 - g)/g
    .GR50 <- GRPsiLin.fun(Psi, Psib, thetaH)
    .GR50 <- ifelse(.GR50>0, .GR50, 0)
    res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
    EDp <- 1/res
    .GR <- EDp
    EDder <- c(-(1/.GR * ((1/.GR50)/(Psi - Psib)/(1/.GR50))/(1/.GR50)^2),
        0,
        -(1/.GR * (1/(Psi - Psib)/(1/.GR50))/(1/.GR)^2),
        -(.GR * log((1 - g)/g) )/ (b^2)
      )
  } }
return(list(EDp, EDder))
}

deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HTE1na.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4])
  d1.2 <- HTE1na.fun(x[,1], x[,2], (parm[,1] + 10e-6), parm[,2], parm[,3],
                   parm[,4])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTE1na.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4])
  d2.2 <- HTE1na.fun(x[,1], x[,2], parm[,1], (parm[,2] + 10e-6), parm[,3],
                   parm[,4])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTE1na.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4])
  d3.2 <- HTE1na.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6),
                   parm[,4])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- HTE1na.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4])
  d4.2 <- HTE1na.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   (parm[,4] + 10e-6))
  d4 <- (d4.2 - d4.1)/10e-6

  cbind(d1, d2, d3, d4)
}

text <- "Hydro-time model with shifted exponential for Pmax and linear model for GR50"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text, edfct=GR, deriv1=deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
