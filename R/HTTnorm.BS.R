HTTnorm.BS.fun <- Vectorize(function(time, Psi, Temp, thetaHT, Tb, To, Psib50, Kt, sigmaPsib){
  # Da Mohsen et al., 2017 - Psib50 decreases with Temperature
  # for any T > To
  t1 <- ifelse(Temp < To, To, Temp)
  .germ1 <- thetaHT/((Temp - Tb)*time)
  .germ2 <- Psi - Psib50 - Kt*(t1 - To) - .germ1
  .germ3 <- ifelse(Temp < Tb, 0, .germ2)
  .germ4 <- .germ2/sigmaPsib
  germ <- pnorm(.germ4)
  ifelse(Temp < Tb, 0, germ)
  })

"HTTnorm.BS" <- function(){
fct <- function(x, parm){
  time <- x[,1]; Psi <- x[,2]; Temp <- x[,3]
  thetaHT <- parm[,1]; Tb <- parm[,2]; To <- parm[,3];
  Psib50 <- parm[,4]
  Kt <- parm[,5]; sigmaPsib <- parm[,6]
  HTTnorm.BS.fun(time, Psi, Temp, thetaHT, Tb, To, Psib50, Kt, sigmaPsib)
}
text <- "Hydrothermal-time model with normal distribution of Psib (Bradford et al., 2002)"
names <- c("thetaHT", "Tb", "To", "Psib50", "Kt", "sigmaPsib")
ss <- function(data){
    # x1 <- phalaris$timeAf
    # x2 <- phalaris$water
    # x3 <- phalaris$temp
    # y <- phalaris$propCum
    # data <- data.frame(x1,x2,x3,y)
    # data <- data[order(data[,3], data[,2], data[,1]), ]
    #
    # x1 <- data[, 1]
    # x2 <- data[, 2]
    # x3 <- data[, 3]
    # y <- data[, 4]

    # temp <- drm(y ~ x1 + x2, fct=HTnorm(), curveid=x3)
    # coef(temp)
    # nLev <- length(levels(factor(x3)))
    # sigmaPsib <- mean(coef(temp)[(2*nLev + 1):(3*nLev)])
    # temp3 <- lm(1/coef(temp)[1:nLev] ~ as.numeric(levels(factor(x3))))
    # Tb <- - coef(temp3)[1] / coef(temp3)[2]
    # thetaHT <- - 1/coef(temp3)[1]
    # temp2 <- lm(coef(temp)[(nLev + 1):(2*nLev)] ~ I(as.numeric(levels(factor(x3)))-Tb))
    # Kt <- coef(temp2)[2]
    # Psib50 <- coef(temp2)[1]
    # thetaHT <- 357;
    # Tb <- 4; To <- 33; Psib50 <- -1.1;
    # Kt <- 0.26; sigmaPsib <- 0.5

  # return(c(thetaHT, Tb, To, Psib50, Kt, sigmaPsib))
}

deriv1 <- function(x, parm){
  #Approximation by using finite differences

  d1.1 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d1.2 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3],(parm[,1] + 10e-6), parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d2.2 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], (parm[,2] + 10e-6), parm[,3], parm[,4], parm[,5], parm[,6])
  d2 <- (d2.2 - d2.1)/10e-6

  d3.1 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d3.2 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], (parm[,3] + 10e-6), parm[,4], parm[,5], parm[,6])
  d3 <- (d3.2 - d3.1)/10e-6

  d4.1 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d4.2 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], (parm[,4] + 10e-6), parm[,5], parm[,6])
  d4 <- (d4.2 - d4.1)/10e-6

  d5.1 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d5.2 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], (parm[,5] + 10e-6), parm[,6])
  d5 <- (d5.2 - d5.1)/10e-6

  d6.1 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], parm[,6])
  d6.2 <- HTTnorm.BS.fun(x[,1], x[,2], x[,3], parm[,1], parm[,2], parm[,3], parm[,4], parm[,5], (parm[,6] + 10e-6))
  d6 <- (d6.2 - d6.1)/10e-6

  cbind(d1, d2, d3, d4, d5, d6)
}
name <- "HTTnorm.BS"
returnList <- list(fct=fct, names=names, text=text, deriv1 = deriv1,
                   name = name) #ssfct=ss,
class(returnList) <- "drcMean"
invisible(returnList)
}
