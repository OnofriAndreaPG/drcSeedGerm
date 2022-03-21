HTTnormRF.fun <- function(time, Psi, Temp, thetaHT, Tb, Td, Psib50, Kt, sigmaPsib){
  # Da Mohsen et al., 2017
  t1 <- ifelse(Temp < Td, Td, Temp)
  .germ1 <- thetaHT/((Temp - Tb)*time)
  .germ2 <- Psi - Psib50 - Kt*(t1 - Td) - .germ1
  .germ3 <- ifelse(Temp < Tb, 0, .germ2)
  .germ4 <- .germ2/sigmaPsib
  germ <- pnorm(.germ4)
  ifelse(Temp < Tb, 0, germ)
  }

HTTnormRF.fun <- function(time, Psi, Temp, thetaHT, Tb, Td, Psib50, Kt, sigmaPsib){
  # Da Mohsen et al., 2017
  t1 <- ifelse(Temp < Td, Td, Temp)
  .germ1 <- thetaHT/((Temp - Tb)*time)
  .germ2 <- Psi - Psib50 - Kt*(t1 - Td) - .germ1
  .germ3 <- ifelse(Temp < Tb, 0, .germ2)
  .germ4 <- .germ2/sigmaPsib
  germ <- pnorm(.germ4)
  ifelse(Temp < Tb, 0, germ)
  }

HTTnorm.fun <- function(time, Psi, Temp, thetaHT, Tb, Psib50, Kt, sigmaPsib){
  #Da Mohsen et al., 2017
  .germ1 <- thetaHT/((Temp - Tb)*time)
  .germ2 <- Psi - Psib50 - Kt*(Temp - Tb) - .germ1
  .germ3 <- ifelse(Temp < Tb, 0, .germ2)
  .germ4 <- .germ2/sigmaPsib
  germ <- pnorm(.germ4)
  ifelse(Temp < Tb, 0, germ)
  }

"HTTnorm" <- function(){
fct <- function(x, parm){
  time <- x[,1]; Psi <- x[,2]; Temp <- x[,3]
  thetaHT <- parm[,1]; Tb <- parm[,2]; Psib50 <- parm[,3]
  Kt <- parm[,4]; sigmaPsib <- parm[,5]
  HTTnorm.fun(time, Psi, Temp, thetaHT, Tb, Psib50, Kt, sigmaPsib)
}
text <- "Hydrothermal-time model with normal distribution of Psib (Bradford et al., 2002)"
names <- c("thetaHT", "Tb", "Psib50", "Kt", "sigmaPsib")
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

    temp <- drm(y ~ x1 + x2, fct=HTnorm(), curveid=x3)
    coef(temp)
    nLev <- length(levels(factor(x3)))
    sigmaPsib <- mean(coef(temp)[(2*nLev + 1):(3*nLev)])
    temp3 <- lm(1/coef(temp)[1:nLev] ~ as.numeric(levels(factor(x3))))
    Tb <- - coef(temp3)[1] / coef(temp3)[2]
    thetaHT <- - 1/coef(temp3)[1]
    temp2 <- lm(coef(temp)[(nLev + 1):(2*nLev)] ~ I(as.numeric(levels(factor(x3)))-Tb))
    temp2
    Kt <- coef(temp2)[2]
    Psib50 <- coef(temp2)[1]

  return(c(thetaHT, Tb, Psib50, Kt, sigmaPsib))
}
returnList <- list(fct=fct, ssfct=ss, names=names, text=text)
class(returnList) <- "drcMean"
invisible(returnList)
}
