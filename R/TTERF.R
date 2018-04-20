#TTERF: exponential Pmax + GRT.RF
TTERF.fun <- function(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b){
  Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
  GR50 <- GRT.RFb.fun(Temp, Tb, Td, Tc, ThetaT)
  plogis( b * (log(time) - log(1/GR50)) )*Pmax
}

"TTERF" <- function(){
fct <- function(x, parm){
  S <- TTERF.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3], parm[,4],
                parm[,5], parm[,6], parm[,7])
  return(S)
}
names <- c("G", "Tc", "sigmaTc", "Td", "Tb", "ThetaT", "b")
ss <- function(data){
  #data <- data.frame(dataset$timeAf, dataset$Temp, dataset$Prop)
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

  modI <- drm(y ~ x1, fct=LL.3(), curveid=TempF, pmodels=list(~1,~TempF-1,~TempF-1), data=dataset_cum)
  psiLevels <- as.numeric(levels(TempF))
  b <- - coef(modI)[1]
  Pmax <- coef(modI)[2:(length(psiLevels)+1)]
  modPmax <- drm(Pmax ~ psiLevels, fct=PmaxT1())
  G <- coef(modPmax)[1]; Tc1 <- coef(modPmax)[2]; sigmaTc <- coef(modPmax)[3]

  GR50 <- 1/coef(modI)[(length(psiLevels)+2):length(coef(modI))]
  modGR <- drm(GR50 ~ psiLevels, fct=GRT.RFb())
  Tb <- coef(modGR)[1]; Td <- coef(modGR)[2];
  Tc2 <- coef(modGR)[3]; ThetaT <- coef(modGR)[4]
  Tc <- mean(c(Tc1, Tc2))

  #temp <- c(G, Tc, sigmaTc, Td, Tb, ThetaT, b)
  #print(temp)
  return(c(G, Tc, sigmaTc, Td, Tb, ThetaT, b)) }
text <- "Thermal-time model with shifted exponential for Pmax and Rowse-Finch-Savage model for GR50"
GR <- function(parms, respl, reference="control", type="relative", Psi, Temp){
   G <- as.numeric(parms[1]); Tc <- as.numeric(parms[2])
   sigmaTc <- as.numeric(parms[3]);
   Td <- as.numeric(parms[4]); Tb <- as.numeric(parms[5])
   ThetaT <- as.numeric(parms[6]); b <- as.numeric(parms[7])
   #G, Tc, sigmaTc, Td, Tb, ThetaT, b
   g <- respl/100
  if(type=="absolute"){
        gra <- function(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b) {
       .Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
       .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
       .temp2 <- (.Pmax - g)/g
       .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
       .GR50 <- GRT.RFb.fun(Temp, Tb, Td, Tc, ThetaT)
       .GR50 <- ifelse(.GR50>0, .GR50, 0)
       .b <- b
       res <- as.numeric( exp( - (1/.b)*log(.temp2) + log(1/.GR50) ) )
       .GR <- 1/res
       .GR}
    EDp <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)

    #Beginning of derivatives (finite differences)
    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G + 10e-6, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc + 10e-6, sigmaTc, Td, Tb, ThetaT, b)
    d2 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc + 10e-6, Td, Tb, ThetaT, b)
    d3 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td + 10e-6, Tb, ThetaT, b)
    d4 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb + 10e-6, ThetaT, b)
    d5 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT + 10e-6, b)
    d6 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b + 10e-6)
    d7 <- (d1.2 - d1.1)/10e-6
    #End of derivatives

    EDder <- c(d1,d2,d3,d4,d5,d6,d7)
  } else{ if(type=="relative") {
    gra <- function(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b) {
       .Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
       .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
       .temp2 <- (1 - g)/g
       .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
       .GR50 <- GRT.RFb.fun(Temp, Tb, Td, Tc, ThetaT)
       .GR50 <- ifelse(.GR50>0, .GR50, 0)
       .b <- b
       res <- as.numeric( exp( - (1/.b)*log(.temp2) + log(1/.GR50) ) )
       .GR <- 1/res
       .GR}
    EDp <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)

    #Beginning of derivatives (finite differences)
    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G + 10e-6, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc + 10e-6, sigmaTc, Td, Tb, ThetaT, b)
    d2 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc + 10e-6, Td, Tb, ThetaT, b)
    d3 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td + 10e-6, Tb, ThetaT, b)
    d4 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb + 10e-6, ThetaT, b)
    d5 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT + 10e-6, b)
    d6 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
    d1.2 <- gra(Temp, g, G, Tc, sigmaTc, Td, Tb, ThetaT, b + 10e-6)
    d7 <- (d1.2 - d1.1)/10e-6
    #End of derivatives

    EDder <- c(d1,d2,d3,d4,d5,d6,d7)
  } }
return(list(EDp, EDder))
}
deriv1 <- function(x, parms){
  #Approximation by using finite differences
  time <- x[,1]; Temp <- x[,2]
  G <-  as.numeric(parms[,1]); Tc <- as.numeric(parms[,2]); sigmaTc <- as.numeric(parms[,3])
  Td <- as.numeric(parms[,4]); Tb <- as.numeric(parms[,5]); ThetaT <-  as.numeric(parms[,6])
  b <- as.numeric(parms[,7])

  d1.1 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
  d1.2 <- TTERF.fun(time, Temp, (G + 10e-6), Tc , sigmaTc, Td, Tb, ThetaT, b)
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
  d2.2 <- TTERF.fun(time, Temp, G, (Tc  + 10e-6), sigmaTc, Td, Tb, ThetaT, b)
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
  d3.2 <- TTERF.fun(time, Temp, G, Tc, (sigmaTc + 10e-6), Td, Tb, ThetaT, b)
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
  d4.2 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, (Td + 10e-6) , Tb, ThetaT, b)
  d4 <- (d4.2 - d4.1)/10e-6

  d5.1 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
  d5.2 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, (Tb + 10e-6) , ThetaT, b)
  d5 <- (d5.2 - d5.1)/10e-6

  d6.1 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
  d6.2 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, (ThetaT + 10e-6) , b)
  d6 <- (d6.2 - d6.1)/10e-6

  d7.1 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
  d7.2 <- TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, (b + 10e-6))
  d7 <- (d7.2 - d7.1)/10e-6

  cbind(d1, d2, d3, d4, d5, d6, d7)
}

returnList <- list(fct=fct, ssfct=ss, names=names, text=text, edfct=GR, deriv1=deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}

TTERFgr <- function(Temp, respl, parms, relative=T) {
   G <- as.numeric(parms[1]); Tc <- as.numeric(parms[2])
   sigmaTc <- as.numeric(parms[3]);
   Td <- as.numeric(parms[4]); Tb <- as.numeric(parms[5])
   ThetaT <- as.numeric(parms[6]); b <- as.numeric(parms[7])
   g <- respl/100
   .Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
   .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
   .temp2 <- (1 - g)/g
   .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
   .GR50 <- GRT.RFb.fun(Temp, Tb, Td, Tc, ThetaT)
   .GR50 <- ifelse(.GR50>0, .GR50, 0)
   .b <- b
   res <- as.numeric( exp( - (1/.b)*log(.temp2) + log(1/.GR50) ) )
   .GR <- 1/res
   .GR}
TTERFga <- function(Temp, respl, parms, relative=T) {
   G <- as.numeric(parms[1]); Tc <- as.numeric(parms[2])
   sigmaTc <- as.numeric(parms[3]);
   Td <- as.numeric(parms[4]); Tb <- as.numeric(parms[5])
   ThetaT <- as.numeric(parms[6]); b <- as.numeric(parms[7])
   g <- respl/100
   .Pmax <- PmaxT1.fun(Temp, G, Tc, sigmaTc)
   .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
   .temp2 <- (.Pmax - g)/g
   .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
   .GR50 <- GRT.RFb.fun(Temp, Tb, Td, Tc, ThetaT)
   .GR50 <- ifelse(.GR50>0, .GR50, 0)
   .b <- b
   res <- as.numeric( exp( - (1/.b)*log(.temp2) + log(1/.GR50) ) )
   .GR <- 1/res
   .GR}
