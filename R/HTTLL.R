#Hydro-Thermal time models
HTTLL.fun <- function(time, Psi, Temp, thetaHT, Tb, Psib50, Kt, delta, sigmaPsib){
  #Da Mohsen et al., 2017
  .germ1 <- thetaHT/((Temp - Tb)*time)
  .germ2 <- Psi - .germ1 - Kt*(Temp - Tb) + delta
  .germ2 <- ifelse(.germ2>0, .germ2, 0.0000001)
  .germ3 <- .germ2/(Psib50 + delta)
  .germ3 <- ifelse(.germ3>0, .germ3, 0.0000001)
  .germ4 <- log(.germ3)
  germ <- 1/(1 + exp(-(.germ4/sigmaPsib)))
  germ
  #ifelse(T > Tb, germ, 0)
}

"HTTLL" <- function(){
fct <- function(x, parm){
  time <- x[,1]; Psi <- x[,2]; Temp <- x[,3]
  thetaHT <- parm[,1]; Tb <- parm[,2]; Psib50 <- parm[,3]
  Kt <- parm[,4]; delta <- parm[,5]; sigmaPsib <- parm[,6]
  HTTLL.fun(time, Psi, Temp, thetaHT, Tb, Psib50, Kt, delta, sigmaPsib)
}
names <- c("thetaHT", "Tb", "Psib50", "Kt", "delta", "sigmaPsib")
ss <- function(data){
    # x1 <- phalaris$timeAf
    # x2 <- phalaris$water
    # x3 <- phalaris$temp
    # y <- phalaris$propCum
    # data <- data.frame(x1,x2,x3,y)
    data <- data[order(data[,3], data[,2], data[,1]), ]
    #
    x1 <- data[, 1]
    x2 <- data[, 2]
    x3 <- data[, 3]
    y <- data[, 4]
    temp <- drm(y ~ x1 + x2, fct=HTLL(), curveid=x3)
    #coef(temp)
    nLev <- length(levels(factor(x3)))
    sigmaPsib <- mean(coef(temp)[(3*nLev + 1):(4*nLev)])
    delta <- mean(coef(temp)[(1*nLev + 1):(2*nLev)])
    temp3 <- lm(1/coef(temp)[1:nLev] ~ as.numeric(levels(factor(x3))))
    Tb <- - coef(temp3)[1] / coef(temp3)[2]
    thetaHT <- - 1/coef(temp3)[1]
    temp2 <- lm(coef(temp)[(2*nLev + 1):(3*nLev)] ~ I(as.numeric(levels(factor(x3)))-Tb))
    temp2
    Kt <- coef(temp2)[2]
    Psib50 <- coef(temp2)[1]
    print(c(thetaHT, Tb, Psib50, Kt, delta, sigmaPsib))

  return(c(thetaHT, Tb, Psib50, Kt, delta, sigmaPsib))
}
text <- "Hydrothermal-time model (Mesgaran et al., 2017)"
returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}
