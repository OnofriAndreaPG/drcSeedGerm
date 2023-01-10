HTTnorm.M.fun <- Vectorize(function(time, Psi, Temp, thetaHT, Tb, Psib50, Kt, sigmaPsib){
  # HTT model based on Psi and Temperature
  # Da Mohsen et al., 2017 - Psib50 decreases with Temperature
  # for any T > Tb
  .cond1 <- ifelse(Temp < Tb, 0, Temp - Tb)
  .germ1 <- thetaHT/(.cond1 * time)
  .germ2 <- Psi - Psib50 - Kt*(Temp - Tb) - .germ1
  .germ3 <- ifelse(Temp < Tb, 0, .germ2)
  .germ4 <- .germ2/sigmaPsib
  germ <- pnorm(.germ4)
  # print(germ)
  germ
  })

"HTTnorm.M" <- function(){
fct <- function(x, parm){
  time <- x[,1]; Psi <- x[,2]; Temp <- x[,3]
  thetaHT <- parm[,1]; Tb <- parm[,2]; Psib50 <- parm[,3]
  Kt <- parm[,4]; sigmaPsib <- parm[,5]
  HTTnorm.M.fun(time, Psi, Temp, thetaHT, Tb, Psib50, Kt, sigmaPsib)
}
text <- "Hydrothermal-time model with normal distribution of Psib (Bradford et al., 2002)"
names <- c("thetaHT", "Tb", "Psib50", "Kt", "sigmaPsib")
ss <- function(data){
    # x1 <- phalaris$timeAf
    # x2 <- phalaris$water
    # x3 <- phalaris$temp
    # y <- phalaris$propCum
    # data <- data.frame(x1,x2,x3,y)
    # data <- data[order(data[,3], data[,2], data[,1]), ]
    # #
    # x1 <- data[, 1]
    # x2 <- data[, 2]
    # x3 <- data[, 3]
    # y <- data[, 4]
    #
    # temp <- drm(y ~ x1 + x2, fct=HTnorm(), curveid=x3)
    # coef(temp)
    # nLev <- length(levels(factor(x3)))
    # sigmaPsib <- mean(coef(temp)[(2*nLev + 1):(3*nLev)])
    # temp3 <- lm(1/coef(temp)[1:nLev] ~ as.numeric(levels(factor(x3))))
    # Tb <- - coef(temp3)[1] / coef(temp3)[2]
    # thetaHT <- - 1/coef(temp3)[1]
    # temp2 <- lm(coef(temp)[(nLev + 1):(2*nLev)] ~ I(as.numeric(levels(factor(x3)))-Tb))
    # temp2
    # Kt <- coef(temp2)[2]
    # Psib50 <- coef(temp2)[1]
  thetaHT = 700; Tb = 3; Psib50 = -3; Kt = 0.06; sigmaPsib = 0.3

  return(c(thetaHT, Tb, Psib50, Kt, sigmaPsib))
}
deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HTTnorm.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5])
  d1.2 <- HTTnorm.M.fun(x[,1], x[,2], x[,3],(parm[,1] + 10e-6), parm[,2], parm[,3], parm[,4], parm[,5])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTTnorm.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5])
  d2.2 <- HTTnorm.M.fun(x[,1], x[,2], x[,3], parm[,1], (parm[,2] + 10e-6), parm[,3], parm[,4], parm[,5])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTTnorm.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5])
  d3.2 <- HTTnorm.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], (parm[,3] + 10e-6), parm[,4], parm[,5])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- HTTnorm.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5])
  d4.2 <- HTTnorm.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], (parm[,4] + 10e-6), parm[,5])
  d4 <- (d4.2 - d4.1)/10e-6

  d5.1 <- HTTnorm.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5])
  d5.2 <- HTTnorm.M.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], (parm[,5] + 10e-6))
  d5 <- (d5.2 - d5.1)/10e-6

  cbind(d1, d2, d3, d4, d5)
}
returnList <- list(fct=fct, names=names, text=text, deriv1 = deriv1) # ssfct=ss,
class(returnList) <- "drcMean"
invisible(returnList)
}
