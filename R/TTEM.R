# TTEM - Based on PmaxT1 + GR50M
TTEM.fun <- function(time, Temp, G, Tc, sigmaTc, Tb, ThetaT, b) {
  Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
  GR50 <- GRT.Mb.fun(Temp, Tb, Tc, ThetaT)
  GR50 <- ifelse(GR50<=0, 1e-06, GR50)
  plogis(b * (log(time) - log(1/GR50)) )* Pmax
}

"TTEM" <- function(){
#TT-to-event. With shifted exponential + GR50 Polynomial
fct <- function(x, parm) {
  S <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[,4],
                parm[,5], parm[,6])
}
names <- c("G", "Tc", "sigmaTc", "Tb", "ThetaT", "b")
ss <- function(data){
  data <- subset(data, is.finite(data[,1])==T)
  result <- c()
  TempF <- factor(data[,2])
  for(i in 1:length(levels(TempF))){
    temp <- subset(data, data[,2] == levels(TempF)[i])
    x <- temp[,1]; y <- temp[,3]
    modT <- try(drm(y ~ x, fct=LL.3()), silent=T)
    #modT <- try(nls(y ~ NLSLL.3(x, a, b, c)), silent=T)
    if(inherits(modT, "try-error")) {res <- as.numeric(levels(TempF)[i])
    result <- c(result, res) }
  }
  dataset_cum <- subset(data, is.finite(data[,1])==T)
  if(is.null(result)!=T){
    for(i in 1:length(result)) dataset_cum <- subset(dataset_cum, dataset_cum[,2] != result[i])}
  TempF <- factor(dataset_cum[,2])

  x1 <- dataset_cum[,1]
  x2 <- dataset_cum[,2]
  y <- dataset_cum[,3]

  modI <- drm(y ~ x1, fct=LL.3(), curveid=TempF, pmodels=list(~1,~TempF-1,~TempF-1), data=dataset_cum)
  psiLevels <- as.numeric(levels(TempF))
  b <- -coef(modI)[1]
  Pmax <- coef(modI)[2:(length(psiLevels)+1)]
  modPmax <- drm(Pmax ~ psiLevels, fct=PmaxT1())
  G <- coef(modPmax)[1]; Tc1 <- coef(modPmax)[2]; sigmaTc <- coef(modPmax)[3]

  GR50 <- 1/coef(modI)[(length(psiLevels)+2):length(coef(modI))]
  modGR <- drm(GR50 ~ psiLevels, fct=GRT.Mb())
  Tc2 <- coef(modGR)[2]; Tb <- coef(modGR)[1]; ThetaT <- coef(modGR)[3]
  Tc <- (Tc1 + Tc2)/2
  # temp <- c(G, Tc, sigmaTc, Tb, ThetaT, b)
  # print(temp)
  return(c(G, Tc, sigmaTc, Tb, ThetaT, b)) }
text <- "Thermal-time model with shifted exponential for Pmax and Mesgaran model for GR50"
GR <- function(parms, respl=50, reference="control", type="relative", Temp){
   G <- as.numeric(parms[1]); Tc <- as.numeric(parms[2])
   sigmaTc <- as.numeric(parms[3]);
   Tb <- as.numeric(parms[4]); ThetaT <- as.numeric(parms[5])
   b <- as.numeric(parms[6])
   g <- respl/100
  if(type=="absolute"){
    .Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
    .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
    .temp2 <- (.Pmax - g)/g
    .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
    .GR50 <- GRT.Mb.fun(Temp, Tb, Tc, ThetaT)
    .GR50 <- ifelse(.GR50>0, .GR50, 0)
     res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
     EDp <- 1/res
     .GR <- EDp
     EDder1 <- 1/.GR * ((1/b) * ((1/g) * (.Pmax/G)/((1/g) * (.Pmax - g))))/exp(-(1/b) *log((1/g) * (.Pmax - g)) + log(1/.GR50))^2
     EDder2 <- 1/.GR * (ThetaT * (( (Temp - Tb)/(Tc - Tb))^2 )/(.GR50*ThetaT)^2/(1/.GR50) + (1/b) * ( (G * (exp(-(Tc - Temp)/sigmaTc) * (1/sigmaTc)))/((.Pmax - g))))/(1/.GR)^2
     EDder3 <- -(1/.GR * ((1/b) * ((1/g) * (G * (exp(-(Tc - Temp)/sigmaTc) * ((Tc - Temp)/sigmaTc^2)))/((1/g) * (.Pmax - g))))/(1/.GR)^2)
     EDder4 <- (1/.GR) * (ThetaT * ((1/(Tc - Tb) - (Temp - Tb)/(Tc - Tb)^2) * (Temp - Tb) - (1 - (Temp - Tb)/(Tc - Tb)))/((1 - (Temp - Tb)/(Tc - Tb)) * (Temp - Tb))^2/(1/.GR50))/(1/.GR)^2
     EDder5 <- -(1/.GR * (1/((1 - (Temp - Tb)/(Tc - Tb)) * (Temp - Tb))/(1/.GR50))/(1/.GR)^2)
     EDder6 <- -(1/.GR * (1/b^2 * log((1/g) * ((G * (1 - exp(-(Tc - Temp)/sigmaTc))) - g)))/(1/.GR)^2)
     EDder <- c(EDder1, EDder2, EDder3, EDder4, EDder5, EDder6)
  } else{ if(type=="relative") {
    .Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
    .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
    .temp2 <- (1 - g)/g
    .GR50 <- GRT.Mb.fun(Temp, Tb, Tc, ThetaT)
    .GR50 <- ifelse(.GR50>0, .GR50, 0)
     res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
     EDp <- 1/res
     .GR <- EDp
     EDder1 <- 0
     EDder2 <- (1/.GR) * (ThetaT * (((Temp - Tb)/(Tc - Tb))^2 )/(.GR50*ThetaT)^2/(1/.GR50))/(1/.GR)^2
     EDder3 <- 0
     EDder4 <- (1/.GR) * (ThetaT * ((1/(Tc - Tb) - (Temp - Tb)/(Tc - Tb)^2) * (Temp - Tb) - (1 - (Temp - Tb)/(Tc - Tb)))/(.GR50*ThetaT)^2/(1/.GR50))/(1/.GR)^2
     EDder5 <- -(1/.GR * (1/((1 - (Temp - Tb)/(Tc - Tb)) * (Temp - Tb))/(1/.GR50))/(1/.GR)^2)
     EDder6 <- -(1/.GR * (1/b^2 * log((1/g) * (1 - g)))/(1/.GR)^2)
     EDder <- c(EDder1, EDder2, EDder3, EDder4, EDder5, EDder6)
  } }
return(list(EDp, EDder))
}
deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6])
  d1.2 <- TTEM.fun(x[,1], x[,2], (parm[,1] + 10e-6), parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6])
  d2.2 <- TTEM.fun(x[,1], x[,2], parm[,1], (parm[,2] + 10e-6), parm[,3],
                   parm[,4], parm[,5], parm[,6])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6])
  d3.2 <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6),
                   parm[,4], parm[,5], parm[,6])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6])
  d4.2 <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   (parm[,4] + 10e-6), parm[,5], parm[,6])
  d4 <- (d4.2 - d4.1)/10e-6

  d5.1 <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6])
  d5.2 <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], (parm[,5] + 10e-6), parm[,6])
  d5 <- (d5.2 - d5.1)/10e-6

  d6.1 <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6])
  d6.2 <- TTEM.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], (parm[,6] + 10e-6))
  d6 <- (d6.2 - d6.1)/10e-6

  cbind(d1, d2, d3, d4, d5, d6)
}

returnList <- list(fct=fct, ssfct=ss, names=names, text=text, edfct=GR, deriv1 = deriv1)
class(returnList) <- "drcMean"
return(returnList)
}
