#TTERFc: exponential Pmax + GRT.RF + linear sigma
TTERFc.fun <- function(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s){
  Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
  GR50 <- GRT.RFb.fun(Temp, Tb, Td, Tc, ThetaT)
  sigma <- (1/b0) + s*(Temp - Tb)
  plogis( (1/sigma) * (log(time) - log(1/GR50)) )*Pmax
}

"TTERFc" <- function(){
fct <- function(x, parm){
  S <- TTERFc.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[,4],
                parm[,5], parm[,6], parm[,7], parm[,8])
  return(S)
}
names <- c("G", "Tc", "sigmaTc", "Td", "Tb", "ThetaT", "b0", "s")
ss <- function(data){
  data <- subset(data, is.finite(data[,1])==T)
  result <- c()
  TempF <- factor(data[,2])
  for(i in 1:length(levels(TempF))){
    temp <- subset(data, data[,2] == levels(TempF)[i])
    x <- temp[,1]; y <- temp[,3]
    modT <- try(drm(y ~ x, fct=LL.3()), silent=T)
    #modT <- try(nls(y ~ NLSLL.3(x, a, b, c)), silent=T)
    if(class(modT) == "try-error") {res <- as.numeric(levels(TempF)[i])
    result <- c(result, res) }
  }

  dataset_cum <- subset(data, is.finite(data[,1])==T)
  if(is.null(result)!=T){
    for(i in 1:length(result)) dataset_cum <- subset(dataset_cum, dataset_cum[,2] != result[i])}
  TempF <- factor(dataset_cum[,2])

  x1 <- dataset_cum[,1]
  x2 <- dataset_cum[,2]
  y <- dataset_cum[,3]

  modI <- drm(y ~ x1, fct=LL.3(), curveid=TempF, pmodels=list(~TempF-1,~TempF-1,~TempF-1), data=dataset_cum)
  psiLevels <- as.numeric(levels(TempF))
  sigma <- -1/coef(modI)[1:length(psiLevels)]
  Pmax <- coef(modI)[(length(psiLevels) + 1):(length(psiLevels)*2)]
  modPmax <- drm(Pmax ~ psiLevels, fct=PmaxT1())
  G <- coef(modPmax)[1]; Tc1 <- coef(modPmax)[2]; sigmaTc <- coef(modPmax)[3]

  GR50 <- 1/coef(modI)[(length(psiLevels)*2+1):length(coef(modI))]
  modGR <- drm(GR50 ~ psiLevels, fct=GRT.RFb())
  Tb <- coef(modGR)[1]; Td <- coef(modGR)[2];
  Tc2 <- coef(modGR)[3]; ThetaT <- coef(modGR)[4]
  Tc <- mean(c(Tc1, Tc2))

  modSigma <- lm(sigma ~ I(psiLevels - Tb))
  b0 <- 1/coef(modSigma)[1]; s <- coef(modSigma)[2]
  #temp <- c(G, Tc, sigmaTc, k, Tb, ThetaT, b0, s)
  #print(temp)
  return(c(G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)) }
text <- "Thermal-time model with shifted exponential for Pmax, Rowse-Finch-Savage model for GR50, linear for sigma"
deriv1 <- function(x, parms){
  #Approximation by using finite differences
  time <- x[,1]; Temp <- x[,2]
  G <-  as.numeric(parms[,1]); Tc <- as.numeric(parms[,2]); sigmaTc <- as.numeric(parms[,3])
  Td <- as.numeric(parms[,4]); Tb <- as.numeric(parms[,5]); ThetaT <-  as.numeric(parms[,6])
  b0 <- as.numeric(parms[,7]); s <-  as.numeric(parms[,8])

  d1.1 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
  d1.2 <- TTERFc.fun(time, Temp, (G + 10e-6), Tc , sigmaTc, Td, Tb, ThetaT, b0, s)
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
  d2.2 <- TTERFc.fun(time, Temp, G, (Tc  + 10e-6), sigmaTc, Td, Tb, ThetaT, b0, s)
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
  d3.2 <- TTERFc.fun(time, Temp, G, Tc, (sigmaTc + 10e-6), Td, Tb, ThetaT, b0, s)
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
  d4.2 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, (Td + 10e-6) , Tb, ThetaT, b0, s)
  d4 <- (d4.2 - d4.1)/10e-6

  d5.1 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
  d5.2 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, (Tb + 10e-6) , ThetaT, b0, s)
  d5 <- (d5.2 - d5.1)/10e-6

  d6.1 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
  d6.2 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, (ThetaT + 10e-6) , b0, s)
  d6 <- (d6.2 - d6.1)/10e-6

  d7.1 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
  d7.2 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, (b0 + 10e-6) , s)
  d7 <- (d7.2 - d7.1)/10e-6

  d8.1 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
  d8.2 <- TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, (s + 10e-6) )
  d8 <- (d8.2 - d8.1)/10e-6

  cbind(d1, d2, d3, d4, d5, d6, d7, d8)
}
GR <- function(parms, respl, reference="control", type="relative", Psi, Temp){
   G <- as.numeric(parms[1]); Tc <- as.numeric(parms[2])
   sigmaTc <- as.numeric(parms[3]);
   Td <- as.numeric(parms[4]); Tb <- as.numeric(parms[5])
   ThetaT <- as.numeric(parms[6]); b0 <- as.numeric(parms[7])
   s <- as.numeric(parms[8])
   #G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s
   g <- respl/100
  if(type=="absolute"){
        gra <- function(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s) {
       .Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
       .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
       .temp2 <- (.Pmax - g)/g
       .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
       .GR50 <- GRT.RFb.fun(Temp, Tb, Td, Tc, ThetaT)
       .GR50 <- ifelse(.GR50>0, .GR50, 0)
       .b <- 1 / ((1/b0) + s*(Temp - Tb))
       res <- as.numeric( exp( - (1/.b)*log(.temp2) + log(1/.GR50) ) )
       .GR <- 1/res
       .GR}
    EDp <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)

    #Beginning of derivatives (finite differences)
    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G + 10e-6, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc + 10e-6, sigmaTc, Td, Tb, ThetaT, b0, s)
    d2 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc + 10e-6, Td, Tb, ThetaT, b0, s)
    d3 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td + 10e-6, Tb, ThetaT, b0, s)
    d4 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb + 10e-6, ThetaT, b0, s)
    d5 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT + 10e-6, b0, s)
    d6 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0 + 10e-6, s)
    d7 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s + 10e-6)
    d8 <- (d1.2 - d1.1)/10e-6
    #End of derivatives
    EDder <- c(d1,d2,d3,d4,d5,d6,d7,d8)
  } else{ if(type=="relative") {
    gra <- function(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s) {
       .Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
       .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
       .temp2 <- (1 - g)/g
       .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
       .GR50 <- GRT.RFb.fun(Temp, Tb, Td, Tc, ThetaT)
       .GR50 <- ifelse(.GR50>0, .GR50, 0)
       .b <- 1 / ((1/b0) + s*(Temp - Tb))
       res <- as.numeric( exp( - (1/.b)*log(.temp2) + log(1/.GR50) ) )
       .GR <- 1/res
       .GR}
    EDp <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)

    #Beginning of derivatives (finite differences)
    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G + 10e-6, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc + 10e-6, sigmaTc, Td, Tb, ThetaT, b0, s)
    d2 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc + 10e-6, Td, Tb, ThetaT, b0, s)
    d3 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td + 10e-6, Tb, ThetaT, b0, s)
    d4 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb + 10e-6, ThetaT, b0, s)
    d5 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT + 10e-6, b0, s)
    d6 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0 + 10e-6, s)
    d7 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s + 10e-6)
    d8 <- (d1.2 - d1.1)/10e-6
    #End of derivatives
    EDder <- c(d1,d2,d3,d4,d5,d6,d7,d8)
  } }
return(list(EDp, EDder))
}
returnList <- list(fct=fct, ssfct=ss, names=names, deriv1=deriv1, text=text, edfct=GR)
class(returnList) <- "drcMean"
invisible(returnList)
}
