LL.dist.fun <- function(x, b, d, e){
  # With link function ?
  # d <- 1 / (1 + exp (- d)) # invLink
  # b <- exp(b) # invLink
  resp <- plogis(b * (log(x) - log(e)) ) * d
  resp
}
LL.inv.fun <- function(y, b, d, e){
  # With link function ?
  # d <- 1 / (1 + exp (- d)) # invLink
  # b <- exp(b) # invLink
  e * ( ( (d - y)/y ) ^ (-1/b) )
}
LL.EDr.fun <- function(g, b, d, e){
  # With link function ?
  # d <- 1 / (1 + exp (- d)) # invLink
  # b <- exp(b) # invLink
  e * ( ( (1 - g)/g ) ^ (-1/b) )
}
LL.EDa.fun <- function(g, b, d, e){
  # With link function ?
  # d <- 1 / (1 + exp (- d)) # invLink
  # b <- exp(b) # invLink
  e * ( ( (d - g)/g ) ^ (-1/b) )
}

LL.dist <- function(fixed = c(NA, NA, NA), names = c("b", "d", "e")) {
    ## Checking arguments
    numParm <- 3
    if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

    ## Fixing parameters (using argument 'fixed')
    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]

    ## Defining the non-linear function
    fct <- function(x, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm
      b <- parmMat[, 1]; d <- parmMat[, 2]; e <- parmMat[, 3]
      resp <- LL.dist.fun(x, b, d, e)
      resp
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
      x <- dataf[, 1]
      y <- dataf[, 2]

      if(length(y[y != 0]) < 2){
         b <- 0; e <- 1; d <- 0.01
      } else {
    # print(max(y))
      data.srt <- sortedXyData(x, y)
      d <- NLSstRtAsymptote( data.srt )
      # d2 <- max(y) * 1.01
      e <- NLSstClosestX(data.srt, d/2)

      # print(c(d, d2, e))
      pseudoY <- log((d - y)/(y + 0.00001))
      pseudoX <- log( x + 0.000001 )
      coefs <- coef( lm(pseudoY ~ pseudoX ) )
      #k <- coefs[1];
      b <- -coefs[2]
      #e <- exp(k/b)
    # print(c(b, d, e))
      d <- ifelse(d > 0.96, 0.96, d)
      d <- ifelse(d < 0, 0, d)
      b <- ifelse(b < 0, 0.01, b)
    # print(c(b, d, e))

      # d <- log(d/(1 - d)) # link
      # b <- log(b) # Link
  # print(c(b, d, e))
      }
      return(c(b, d, e)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives
    deriv1 <- function(x, parm){
    #Approximation by using finite differences
      d1.1 <- LL.dist.fun(x, parm[,1], parm[,2], parm[,3])
      d1.2 <- LL.dist.fun(x, (parm[,1] + 10e-6), parm[,2], parm[,3])
      d1 <- (d1.2 - d1.1)/10e-6

      d2.1 <- LL.dist.fun(x, parm[,1], parm[,2], parm[,3])
      d2.2 <- LL.dist.fun(x, parm[,1], (parm[,2] + 10e-6), parm[,3])
      d2 <- (d2.2 - d2.1)/10e-6

      d3.1 <- LL.dist.fun(x, parm[,1], parm[,2], parm[,3])
      d3.2 <- LL.dist.fun(x, parm[,1], parm[,2], (parm[,3] + 10e-6))
      d3 <- (d3.2 - d3.1)/10e-6

      cbind(d1, d2, d3)
    }

    ## Defining the ED function
    EDfct <- function(parms, respl, reference="control", type="relative"){
    parm <- as.numeric(parms)
    # b <- as.numeric(parms[1])
    # d <- as.numeric(parms[2])
    # e <- as.numeric(parms[3])
    g <- respl/100
  if(type=="absolute"){

    EDp <- LL.EDa.fun(g, parm[1], parm[2], parm[3])

    #Approximation of derivatives(finite differences)
    d1.1 <- LL.EDa.fun(g, parm[1], parm[2], parm[3])
    d1.2 <- LL.EDa.fun(g, (parm[1] + 10e-6), parm[2], parm[3])
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- LL.EDa.fun(g, parm[1], parm[2], parm[3])
    d2.2 <- LL.EDa.fun(g, parm[1], (parm[2] + 10e-6), parm[3])
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- LL.EDa.fun(g, parm[1], parm[2], parm[3])
    d3.2 <- LL.EDa.fun(g, parm[1], parm[2], (parm[3] + 10e-6))
    d3<- (d3.2 - d3.1)/10e-6

    EDder <- c(d1, d2, d3)
  } else{ if(type=="relative") {
    EDp <- LL.EDr.fun(g, parm[1], parm[2], parm[3])

    #Approximation of derivatives(finite differences)
    d1.1 <- LL.EDr.fun(g, parm[1], parm[2], parm[3])
    d1.2 <- LL.EDr.fun(g, (parm[1] + 10e-6), parm[2], parm[3])
    d1 <- (d1.2 - d1.1)/10e-6

    d2.1 <- LL.EDr.fun(g, parm[1], parm[2], parm[3])
    d2.2 <- LL.EDr.fun(g, parm[1], (parm[2] + 10e-6), parm[3])
    d2 <- (d2.2 - d2.1)/10e-6

    d3.1 <- LL.EDr.fun(g, parm[1], parm[2], parm[3])
    d3.2 <- LL.EDr.fun(g, parm[1], parm[2], (parm[3] + 10e-6))
    d3<- (d3.2 - d3.1)/10e-6

    EDder <- c(d1, d2, d3)
    } }
    # print(EDder)
    return(list(EDp, EDder))
    }

    ## Defining the inverse function


    ## Defining descriptive text
    text <- "Log-logistic distribution for germination times"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text,
                       noParm = sum(is.na(fixed)), edfct=EDfct, deriv1=deriv1)

    class(returnList) <- "drcMean"
    invisible(returnList)
  }

