drmHTT <- function(formula, SE, fct,  data, curveid, pmodels, varPower = F){
# Function to fit HTT models to the observed GR values, weighting based on
# standard errorsd observed in the first step (Modified on 9/10/2019)

  MLfun <- function(formula, SE, fct, data, parms, indMat, curveid=NA){

    fr <- model.frame(formula, data)
    X <- model.matrix(fr, data)[,2]
    Y <- model.response(fr, "numeric")
    # if (missing(curveid))  # in case not supplied
    # {
    #     curveid <- rep(1, length(Y))
    # }
    gamma2i <- SE^2
    tau2 <- exp(parms[1])
    totVar <- tau2 + gamma2i

    nRow <- length(indMat[,1])
    nCol <- length(indMat[1,])
    parmVal <- parms[2:length(parms)]
    parmMat <- matrix(parmVal[indMat], nRow, nCol)
    colnames(parmMat) <- colnames(indMat)
    groupLevels <- as.character(curveid)
    pm <- t(parmMat[, groupLevels, drop = FALSE])
    expt <- fct$fct(X, pm)
    res <- Y - expt
    n <- length(res)
    ll <-  mvtnorm::dmvnorm(res, mean = rep(0, n), sigma = diag(totVar), log = T)
    -ll
}
  MLWfun <- function(formula, SE, fct, data, parms, indMat, curveid=NA){
    #Da rifare!!!!!
    fr <- model.frame(formula, data)
    X <- model.matrix(fr, data)[,2]
    Y <- model.response(fr, "numeric")
    w <- 1/vi
    parmMat <- matrix(parms, 1, length(parms))
    expt <- fct$fct(X, parmMat)
    res <- Y - expt
    #- ( -length(Y)/2*log(2*pi)-1/2*sum(log(totVar))-
    #      1/2*sum(res^2/totVar ) )
    wls <-   w * (res^2)
    - sum(wls)
  }
  MLfun2 <- function(formula, SE, fct, data, parms, indMat, curveid=NA){
    fr <- model.frame(formula, data)
    X <- model.matrix(fr, data)[,2]
    Y <- model.response(fr, "numeric")
    # if (missing(curveid))  # in case not supplied
    # {
    #     curveid <- rep(1, length(Y))
    # }
    gamma2i <- SE^2
    delta <- parms[1]

    nRow <- length(indMat[,1])
    nCol <- length(indMat[1,])
    parmVal <- parms[2:length(parms)]
    parmMat <- matrix(parmVal[indMat], nRow, nCol)
    colnames(parmMat) <- colnames(indMat)
    groupLevels <- as.character(curveid)

    pm <- t(parmMat[, groupLevels, drop = FALSE])
    expt <- fct$fct(X, pm)
    ll <-  mvtnorm::dmvnorm(Y, mean = expt, sigma = diag(gamma2i + delta*expt + 0.00001), log = T)
    # ll2 <- dnorm(Y, expt, sqrt(gamma2i+tau2*expt+0.00001), log = T)
    # sum(ll2)
    -ll
   }

  #Getting initial values and running a naive regression
  vi <- SE^2

  if (missing(curveid))  # in case not supplied
    {
        curveid <- rep(1, length(data[,1]))
        curveid <<- curveid
        naiveMod <- drm(formula, data=data, fct=fct)
  } else {
        curveid <<- curveid
        naiveMod <- drcSeedGerm::drm(formula, curveid=curveid, fct = fct, data = data,
               pmodels = pmodels )
    }
  MSE <- sum(residuals(naiveMod)^2)/naiveMod$df
  tau2 <- MSE - mean(vi)
  #print(coef(naiveMod))

  if(tau2 <= 0 & varPower == F){
    # Running a simple weighted regression
    parms <- c(naiveMod$coef)
    likfun <- optim(par = parms,
                      fn = MLWfun,
                      method = "BFGS", hessian=T,
                      formula = formula, SE = SE, fct = fct,
                  data = data)

  }else{ if(tau2 > 0 & varPower == F){
    parms <- c(log(tau2), naiveMod$coef)
    indMat <- naiveMod$indexMat
    names(parms) <- c("log-tau2", names(coef(naiveMod)))

    likfun <- optim(par = parms, #lower = c(0, rep(NA, length(parms))),
                      fn = MLfun,
                      method = "BFGS", hessian=T,
                      formula = formula, SE = SE, fct = fct,
                  data=data, curveid = curveid, indMat = indMat)
  } else { if(varPower == T) {
    parms <- c(0.01, naiveMod$coef)
    indMat <- naiveMod$indexMat
    names(parms) <- c("delta", names(coef(naiveMod)))
    likfun <- optim(par = parms,
                fn = MLfun2,
                method = "BFGS", hessian=T,
                formula = formula, SE = SE, fct = fct,
                data = data, curveid = curveid, indMat = indMat)
  }
  }}

  #Get the necessary infos from optim object
  parnames <- names(parms) #c("log-tau2", fct$names)
  #print(parnames)
  parvalues <- as.numeric(likfun$par)
  #print(parvalues)

  logLik <- likfun$value
  vcov <- solve(likfun$hessian)
  parSE1 <- as.numeric( sqrt(diag(vcov)) )

  if(varPower == F){
    tau2 <- exp(parvalues[1])
    tau2.SE <- exp(parvalues[1])*parSE1[1]
    tau <- sqrt(tau2)
    tau.SE <- 1/( 2 * sqrt(tau2) ) * tau2.SE
  } else {
    #niente per ora
  }
    #If I do not want to use the hessian to get SEs and
    # I do not want to fit a weighted regression I can
    # X <- fct$deriv1(naiveMod$dataList$dose, t(parvalues[2:length(parvalues)]))
    # colnames(X) <- fct$names
    # totVar <- exp(parvalues[1]) + vi
    # w <- 1/totVar
    # W <- diag(w)
    # vcovMeta <- solve(t(X) %*% W %*% X) %*%
    #            (t(X) %*% W %*% diag(totVar) %*%
    #             t(W) %*% X) %*% solve(t(X) %*% W %*% X)
    # parSE2 <- as.numeric( sqrt(diag(vcovMeta)) )
    # parSE2 <- c(NA, parSE2)

  coefficients <- data.frame("Estimate" = parvalues, "SE.ML"= parSE1)
  row.names(coefficients) <- parnames
  if(varPower == F){
    tauPar <- data.frame("Estimate" = c(tau2, tau), "SE.ML" = c(tau2.SE, tau.SE))
    row.names(tauPar) <- c("tau2", "tau")
  } else {
    tauPar <- NA
  }
  rm(curveid, envir = globalenv())
  list("coefficients" = coefficients, "vcov" = vcov,
         "logLik" = logLik, "tau" = tauPar, "convergence" = likfun$convergence)
}

meta.drm <- function(formula, SE, fct, data, varPower = F){
# Function to fit HTT models to the observed GR values, weighting based on
# standard errorsd observed in the first step. SUPERSEDED!!!!!!!!!!!!!!!!|

  MLfun <- function(formula, SE, fct, data, parms){
    fr <- model.frame(formula, data)
    X <- model.matrix(fr, data)[,2]
    Y <- model.response(fr, "numeric")
    vi <- SE^2
    tau <- exp(parms[1])
    totVar <- tau + vi
    parmMat <- matrix(parms[2:length(parms)], 1, length(parms[2:length(parms)]))
    expt <- fct$fct(X, parmMat)
    res <- Y - expt
    ll <-  -1/2 * log(2*pi) - 1/2 * log(totVar) - 1/2 * (res^2) / totVar
    -sum(ll)
  }
  MLWfun <- function(formula, SE, fct, data, parms){
    fr <- model.frame(formula, data)
    X <- model.matrix(fr, data)[,2]
    Y <- model.response(fr, "numeric")
    w <- 1/vi
    parmMat <- matrix(parms, 1, length(parms))
    expt <- fct$fct(X, parmMat)
    res <- Y - expt
    #- ( -length(Y)/2*log(2*pi)-1/2*sum(log(totVar))-
    #      1/2*sum(res^2/totVar ) )
    wls <-   w * (res^2)
    -sum(wls)
  }

  MLfun2 <- function(formula, SE, fct, data, parms){
    fr <- model.frame(formula, data)
    X <- model.matrix(fr, data)[,2]
    Y <- model.response(fr, "numeric")
    gamma2i <- SE^2
    delta <- parms[1]
    parmMat <- matrix(parms[2:length(parms)], 1,
                      length(parms[2:length(parms)]))
    expt <- fct$fct(X, parmMat)
    ll <-  mvtnorm::dmvnorm(Y, mean = expt, sigma = diag(gamma2i + delta*expt + 0.00001), log = T)
    # ll2 <- dnorm(Y, expt, sqrt(gamma2i+tau2*expt+0.00001), log = T)
    # sum(ll2)
    -ll
   }


  #Getting initial values
  vi <- SE^2
  naiveMod <- drm(formula, data=data, fct=fct)
  MSE <- sum(residuals(naiveMod)^2)/naiveMod$df
  tau2 <- MSE - mean(vi)


  if(tau2 <= 0 & varPower == F){
    # Running a simple weighted regression
    parms <- c(naiveMod$coef)
    likfun <- optim(par = parms,
                      fn = MLWfun,
                      method = "BFGS", hessian=T,
                      formula = formula, SE = SE, fct = fct,
                  data = data)

  }else{
    if(tau2 > 0 & varPower == F){
    parms <- c(log(tau2), naiveMod$coef)
    names(parms) <- c("log-tau2", fct$names)
    likfun <- optim(par = parms, #lower = c(0, rep(NA, length(parms))),
                      fn = MLfun,
                      method = "BFGS", hessian=T,
                      formula = formula, SE = SE, fct = fct,
                  data=data)
  } else { if(varPower == T) {
    parms <- c(0.1, naiveMod$coef)
    names(parms) <- c("delta", fct$names)
    likfun <- optim(par = parms,
                fn = MLfun2,
                method = "BFGS", hessian=T,
                formula = formula, SE = SE, fct = fct,
                data = data)
  }
  }}

  #Get the necessary infos from optim object
  parnames <- names(parms) #c("log-tau2", fct$names)
  #print(parnames)
  parvalues <- as.numeric(likfun$par)
  #print(parvalues)

  logLik <- likfun$value
  vcov <- solve(likfun$hessian)
  parSE1 <- as.numeric( sqrt(diag(vcov)) )

  if(varPower == F){
    tau2 <- exp(parvalues[1])
    tau2.SE <- exp(parvalues[1])*parSE1[1]
    tau <- sqrt(tau2)
    tau.SE <- 1/( 2 * sqrt(tau2) ) * tau2.SE
  } else {
    #niente per ora
  }
    #If I do not want to use the hessian to get SEs and
    # I do not want to fit a weighted regression I can
    # X <- fct$deriv1(naiveMod$dataList$dose, t(parvalues[2:length(parvalues)]))
    # colnames(X) <- fct$names
    # totVar <- exp(parvalues[1]) + vi
    # w <- 1/totVar
    # W <- diag(w)
    # vcovMeta <- solve(t(X) %*% W %*% X) %*%
    #            (t(X) %*% W %*% diag(totVar) %*%
    #             t(W) %*% X) %*% solve(t(X) %*% W %*% X)
    # parSE2 <- as.numeric( sqrt(diag(vcovMeta)) )
    # parSE2 <- c(NA, parSE2)

  coefficients <- data.frame("Estimate" = parvalues, "SE.ML"= parSE1)
  row.names(coefficients) <- parnames
  if(varPower == F){
    tauPar <- data.frame("Estimate" = c(tau2, tau), "SE.ML" = c(tau2.SE, tau.SE))
    row.names(tauPar) <- c("tau2", "tau")
  } else {
    tauPar <- NA
  }

  list("coefficients" = coefficients, "vcov" = vcov,
         "logLik" = logLik, "tau" = tauPar, "convergence" = likfun$convergence)
}
