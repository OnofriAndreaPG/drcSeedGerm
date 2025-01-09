#Pmax vs Temp and Psi ############################################
PmaxTP1.fun <- function(Psi, Temp, G, Psib50, kt, Tb, sigmaPsi) {
  Pmax <- G * (1 - exp( - (Psi - Psib50 - kt * (Temp - Tb))/sigmaPsi))
  Pmax <- ifelse(Temp < Tb | Psi < (Psib50 - kt * (Temp - Tb)), 0, Pmax)
  Pmax[Pmax<0] <- 0
  return(Pmax) }

GR50TPM.fun <- function(Psi, Temp, k, tb, ThetaHT, Psib) {
  t2 <- ifelse(Temp < tb, tb, Temp)
  psival <- ifelse(Psi - Psib - k*(Temp - tb) > 0, Psi - Psib - k*(Temp - tb), 0)
  GR <- psival * (t2 - tb)/ThetaHT
  return(ifelse(GR < 0 , 0 , GR)) }

HTTEMgr <- function(Psi, Temp, respl, parms) {
   G <- as.numeric(parms[1]); Psib <- as.numeric(parms[2])
   k <- as.numeric(parms[3]);
   Tb <- as.numeric(parms[4]); sigmaPsib <- as.numeric(parms[5])
   ThetaHT <- as.numeric(parms[6]); b <- as.numeric(parms[7])
   g <- respl/100
   .Pmax <- PmaxTP1.fun(Psi, Temp, G, Psib, k, Tb, sigmaPsib)
   .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
   .temp2 <- (1 - g)/g
   .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
   .GR50 <- GR50TPM.fun(Psi, Temp, k, Tb, ThetaHT, Psib)
   .GR50 <- ifelse(.GR50>0, .GR50, 0)
   res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
   .GR <- 1/res
   .GR}
HTTEMga <- function(Psi, Temp, respl, parms) {
   G <- as.numeric(parms[1]); Psib <- as.numeric(parms[2])
   k <- as.numeric(parms[3]);
   Tb <- as.numeric(parms[4]); sigmaPsib <- as.numeric(parms[5])
   ThetaHT <- as.numeric(parms[6]); b <- as.numeric(parms[7])
   g <- respl/100
   .Pmax <- PmaxTP1.fun(Psi, Temp, G, Psib, k, Tb, sigmaPsib)
   .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
   .temp2 <- (.Pmax - g)/g
   .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
   .GR50 <- GR50TPM.fun(Psi, Temp, k, Tb, ThetaHT, Psib)
   .GR50 <- ifelse(.GR50>0, .GR50, 0)
   res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
   .GR <- 1/res
   .GR}
HTTEM.fun <- function(time, Psi, Temp, G, Psib, kt, Tb, sigmaPsib, ThetaHT, b) {
  Pmax <- PmaxTP1.fun(Psi, Temp, G, Psib, kt, Tb, sigmaPsib)
  GR50 <- GR50TPM.fun(Psi, Temp, kt, Tb, ThetaHT,Psib)
  GR50 <- ifelse(GR50<=0, 1e-06, GR50)
  plogis(b * (log(time) - log(1/GR50)) )*Pmax
}

"HTTEM" <- function(){
fct <- function(x, parm) {
  S <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4],
                parm[,5], parm[,6], parm[,7])
}
names <- c("G", "Psib", "kt", "Tb", "sigmaPsib", "ThetaHT", "b")
ss <- function(data){
   G=0.8; Psib=-2; kt=0.05; Tb=3; sigmaPsib=0.2; ThetaHT=2000; b=0.5}
GR <- function(parms, respl, reference="control", type="relative", Psi, Temp){
   G <- as.numeric(parms[1]); Psib <- as.numeric(parms[2])
   k <- as.numeric(parms[3]);
   Tb <- as.numeric(parms[4]); sigmaPsib <- as.numeric(parms[5])
   ThetaHT <- as.numeric(parms[6]); b <- as.numeric(parms[7])
   g <- respl/100
  if(type=="absolute"){
        gra <- function(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b) {
       .Pmax <- PmaxTP1.fun(Psi, Temp, G, Psib, k, Tb, sigmaPsib)
       .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
       .temp2 <- (.Pmax - g)/g
       .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
       .GR50 <- GR50TPM.fun(Psi, Temp, k, Tb, ThetaHT, Psib)
       .GR50 <- ifelse(.GR50>0, .GR50, 0)
       res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
       .GR <- 1/res
       .GR}
    EDp <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)

    #Beginning of derivatives (finite differences) ###########################
    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G+ 10e-6, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib+ 10e-6, k, Tb, sigmaPsib, ThetaHT, b)
    d2 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib, k+ 10e-6, Tb, sigmaPsib, ThetaHT, b)
    d3 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib, k, Tb+ 10e-6, sigmaPsib, ThetaHT, b)
    d4 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib+ 10e-6, ThetaHT, b)
    d5 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT+ 10e-6, b)
    d6 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b+ 10e-6)
    d7 <- (d1.2 - d1.1)/10e-6
    #End of derivatives #########################
    EDder <- c(d1,d2,d3,d4,d5,d6,d7)
  } else{ if(type=="relative") {
    gra <- function(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b) {
       .Pmax <- PmaxTP1.fun(Psi, Temp, G, Psib, k, Tb, sigmaPsib)
       .Pmax <- ifelse(.Pmax > 0, .Pmax, 0)
       .temp2 <- (1 - g)/g
       .temp2 <- ifelse(.temp2 < 0, 0, .temp2)
       .GR50 <- GR50TPM.fun(Psi, Temp, k, Tb, ThetaHT, Psib)
       .GR50 <- ifelse(.GR50>0, .GR50, 0)
       res <- as.numeric( exp( - (1/b)*log(.temp2) + log(1/.GR50) ) )
       .GR <- 1/res
       .GR}
     EDp <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    #Beginning of derivatives (finite differences) ###########################
    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G+ 10e-6, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib+ 10e-6, k, Tb, sigmaPsib, ThetaHT, b)
    d2 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib, k+ 10e-6, Tb, sigmaPsib, ThetaHT, b)
    d3 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib, k, Tb+ 10e-6, sigmaPsib, ThetaHT, b)
    d4 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib+ 10e-6, ThetaHT, b)
    d5 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT+ 10e-6, b)
    d6 <- (d1.2 - d1.1)/10e-6

    d1.1 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b)
    d1.2 <- gra(Psi, Temp, g, G, Psib, k, Tb, sigmaPsib, ThetaHT, b+ 10e-6)
    d7 <- (d1.2 - d1.1)/10e-6
    #End of derivatives #########################
     EDder <- c(d1,d2,d3,d4,d5,d6,d7)
  } }
return(list(EDp, EDder))
}
deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6], parm[,7])
  d1.2 <- HTTEM.fun(x[,1], x[,2], x[,3], (parm[,1] + 10e-6), parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6], parm[,7])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6], parm[,7])
  d2.2 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], (parm[,2] + 10e-6), parm[,3],
                   parm[,4], parm[,5], parm[,6], parm[,7])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6], parm[,7])
  d3.2 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], (parm[,3] + 10e-6),
                   parm[,4], parm[,5], parm[,6], parm[,7])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6], parm[,7])
  d4.2 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                   (parm[,4] + 10e-6), parm[,5], parm[,6], parm[,7])
  d4 <- (d4.2 - d4.1)/10e-6

  d5.1 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                   parm[,4], parm[,5], parm[,6], parm[,7])
  d5.2 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                   parm[,4], (parm[,5] + 10e-6), parm[,6], parm[,7])
  d5 <- (d5.2 - d5.1)/10e-6

  d6.1 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                    parm[,4], parm[,5], parm[,6], parm[,7])
  d6.2 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                    parm[,4], parm[,5], (parm[,6] + 10e-6), parm[,7])
  d6 <- (d6.2 - d6.1)/10e-6

  d7.1 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                    parm[,4], parm[,5], parm[,6], parm[,7])
  d7.2 <- HTTEM.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3],
                    parm[,4], parm[,5], parm[,6], (parm[,7] + 10e-6))
  d7 <- (d7.2 - d7.1)/10e-6

  cbind(d1, d2, d3, d4, d5, d6, d7)
}
text <- "Hydro-thermal-time-model (Mesgaran et al., 2017)"
name <- "HTTEM"
returnList <- list(fct=fct, ssfct=ss, names=names, edfct=GR,
                   text=text, deriv1=deriv1, name = name)
class(returnList) <- "drcMean"
invisible(returnList)
}
