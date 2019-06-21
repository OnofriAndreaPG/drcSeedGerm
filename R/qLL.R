# Quantile function for log-logistic distribution
qLL.2.fun <- function(x, b, e){
  e * ( ( (1 - x)/x ) ^ (-1/b) )  }

"qLL.2" <- function(){
fct <- function(x, parm){
  qLL.2.fun(x, parm[,1], parm[,2]) }
names <- c("b", "e")
text <- "Quantile function for log-logistic distribution"
ss <- function(data){
  x <- data[, 1]
  y <- data[, 2]
  mod <- drm(x ~ y, fct = LL.2())
  b <- -coef(mod)[1]
  e <- coef(mod)[2]
  return( c(b, e) )
}

deriv1 <- function(x, parm){
  #Approximation by using finite differences
  d1.1 <- qLL.2.fun(x, parm[,1], parm[,2])
  d1.2 <- qLL.2.fun(x, (parm[,1] + 10e-6), parm[,2])
  d1 <- (d1.2 - d1.1)/10e-6

  d2.1 <-qLL.2.fun(x, parm[,1], parm[,2])
  d2.2 <- qLL.2.fun(x, parm[,1], (parm[,2] + 10e-6))
  d2 <- (d2.2 - d2.1)/10e-6

  # d3.1 <- HTnorm.fun(x[,1], x[,2], parm[,1], parm[,2], parm[,3])
  # d3.2 <- HTnorm.fun(x[,1], x[,2], parm[,1], parm[,2], (parm[,3] + 10e-6))
  # d3 <- (d3.2 - d3.1)/10e-6

  cbind(d1, d2)
}

returnList <- list(fct=fct, ssfct=ss, names=names, text=text, deriv1=deriv1)
class(returnList) <- "drcMean"
invisible(returnList)
}
