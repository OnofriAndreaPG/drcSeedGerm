drmSG <- function(formula, data, curveid, fct="LL", min=0.1,
                       GRlev = c(10, 30, 50)) {

  fr <- model.frame(formula, data)
  timeBef <- model.matrix(fr, data)[,2]
  timeAf <- model.matrix(fr, data)[,3]
  nSeeds <- model.response(fr, "numeric")

  anName <- deparse(substitute(curveid))  # storing name for later use
  Group <- as.factor( subset(data, select = anName) [,1])

  DataC <- data.frame(Group, timeBef, timeAf, nSeeds)

  result <- data.frame()
  resultES <- data.frame()
  resultLow <- data.frame()
  resultUp <- data.frame()
  pMaxFin <- data.frame()
  coefMod <- data.frame()

  if(fct == "LL"){fct1 <- LL.3(); fct2<- LL.2()
  }else{ if(fct == "W1") { fct1 <- W1.3(); fct2<- W1.2()} }

  #Fitting group by group
  nLev <- length(levels(Group))
    for(i in 1:nLev){
      Gname <- levels(Group)[i]
      dataTemp <- subset(DataC, Group==levels(Group)[i])
      nTot <- sum(dataTemp$nSeeds)
      nGerm <- nTot - sum( dataTemp[dataTemp$timeAf==Inf,]$nSeeds )
      pMaxO <- nGerm/nTot

        #Return pMx0 #################################
        pMaxFin[i, 1] <- Gname; pMaxFin[i, 2] <- nGerm
        pMaxFin[i, 3] <- nTot; pMaxFin[i, 4] <- pMaxO
        colnames(pMaxFin) <- c("Group", "nGerm", "nTot", "pMaxFin")

      nFirst <- sum( dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$nSeeds )
      pFirst <- nFirst/nTot
      tFirst <- dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$timeAf
      tLast <- dataTemp[is.finite(dataTemp$timeAf)==F,]$timeBef

      if(pMaxO < min){
          #If minimum threshold of germination is not reached
          #I assume there is negligible germination (mod = 1)
          #min may be user-defined
          res <-   c(i, nGerm, nTot, pMaxO,1, 0,  10E-6, NA, rep(NA, 18), rep(10E-6, 18))
          resES <- c(i, nGerm, nTot, NA,   1, NA, 10E-6, NA, rep(NA, 18), rep(10E-6, 18))
          resLow<- c(i, nGerm, nTot, NA,   1, NA, NA   , tLast, rep(tLast, 18), rep(1/tLast, 18))
          resUp <- c(i, nGerm, nTot, NA,   1, NA, NA   , NA, rep(Inf, 18), rep(0, 18))

            coefModT <- c(1, NA, NA, NA, NA, NA, NA)
            coefMod <- rbind(coefMod, coefModT)

          result <- rbind(result, res)
          resultES <- rbind(resultES, resES)
          resultLow <- rbind(resultLow, resLow)
          resultUp <- rbind(resultUp, resUp)
          cat(paste("Group ", i, ": pMax lower than minimum threshold", "\n", sep=""))
          next;} else{
      if(pFirst > 0.95){

        #Pmax at first inspection > 0.95. Cannot fit germination model
        #Use midPoint imputation (mod = 2)
        res   <- c(i, nGerm, nTot, pMaxO, 2, NA, pMaxO, tFirst/2, rep(tFirst/2, 18), rep(1/(tFirst/2), 18))
        resES <- c(i, nGerm, nTot, NA,    2, NA, NA, NA,       rep(NA, 18), rep(1/(tFirst/2), 18))
        resLow<- c(i, nGerm, nTot, NA,    2, NA, NA, 0,        rep(0, 18), rep(Inf, 18))
        resUp <- c(i, nGerm, nTot, NA,    2, NA, NA, tFirst,    rep(tFirst, 18), rep(1/tFirst, 18))
        result <- rbind(result, res)
        resultES <- rbind(resultES, resES)
        resultLow <- rbind(resultLow, resLow)
        resultUp <- rbind(resultUp, resUp)
          coefModT <- c(2, NA, NA, NA, NA, NA, NA)
          coefMod <- rbind(coefMod, coefModT)
        cat(paste("Group ", i, ": Germination at first inspection almost complete", "\n", sep=""))
        next; } else{

        #Fit germination model: LL.2 and LL.3
        #A germination medel can be fit. Try LL2 and LL3
        #options(echo=F)
        #fct1 <- paste(fct, ".3()"); fct2 <- paste(fct, ".2()")
        if(fct == "LL") fct1 <- LL.3(); fct2<- LL.2()

        cureMod <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
                      fct = fct1, type = "event",
                      upperl = c(NA, 1, NA)), silent=T)
        cureMod2 <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
                      fct = fct2, type = "event"), silent=T )
        #options(echo=T)
          } }
      #Look at what model is OK
      #class(cureMod); class(cureMod2)
      if(class(cureMod) == "try-error" & class(cureMod2) == "try-error"){

        #No parameteric fit was possible (mod = 3). To be meditated
        cat(paste("Group ", i, ": No parametric fit was possible", "\n", sep=""))

        res <- c(i, nGerm, nTot, pMaxO, NA, NA, NA, NA, rep(NA, 18), rep(NA, 18))
        res <- c(i, NA, NA, NA, NA, NA, NA, NA, rep(NA, 18), rep(NA, 18))
        res <- c(i, NA, NA, NA, NA, NA, NA, NA, rep(NA, 18), rep(NA, 18))
        res <- c(i, NA, NA, NA, NA, NA, NA, NA, rep(NA, 18), rep(NA, 18))
        result <- rbind(result, res)
        resultES <- rbind(resultES, resES)
        resultLow <- rbind(resultLow, resLow)
        resultUp <- rbind(resultUp, resUp)
          coefModT <- c(3, NA, NA, NA, NA, NA, NA)
          coefMod <- rbind(coefMod, coefModT)
        } else {if(class(cureMod) == "try-error" & class(cureMod2) == "drc"){

        #LL.3 could not be fit, but LL.2 was ok
        #mod = 4
        cat(paste("Group ", i, ": ", fct, ".3() could not be fit. ", fct, ".2() is fit instead", "\n", sep=""))
        coefs <- coef(cureMod2)
        coefES <- summary(cureMod2)$coef[,2]
        coefL <- coefs - 2*coefES
        coefU <- coefs + 2*coefES
        RSS <- sum(residuals(cureMod2)^2)

        tg1o <- ED(cureMod2, seq(0.10, 0.90, 0.1), type = "absolute", display=F)
        tg2o <- ED(cureMod2, seq(0.10, 0.90, 0.1), type = "relative", display=F)
        tg1 <- tg1o[,1]; tg1es <- tg1o[,2]; tg1L <- tg1 - 2*tg1es; tg1U <- tg1 + 2*tg1es
        tg2 <- tg2o[,1]; tg2es <- tg2o[,2]; tg2L <- tg2 - 2*tg2es; tg2U <- tg2 + 2*tg2es
        GR1 <- 1/tg1
        GR2 <- 1/tg2
        GR1es <- sqrt( (tg1es^2) * ((-1/tg1)^2)^2 ) #Delta method
        GR2es <- sqrt( (tg2es^2) * ((-1/tg2)^2)^2 )
        GR1L <- GR1 - 2*GR1es; GR1U <- GR1 + 2*GR1es
        GR2L <- GR2 - 2*GR2es; GR2U <- GR2 + 2*GR2es
        tg1[is.na(tg1)] <- Inf
        GR1[is.na(GR1)] <- 1E-6
        tg1es[is.na(tg1es)] <- NA; tg1L[is.na(tg1L)] <- tLast; tg1U[is.na(tg1U)] <- Inf
        GR1es[is.na(GR1es)] <- 1E-6; GR1L[is.na(GR1L)] <- 1/tLast; GR1U[is.na(GR1U)] <- 0

        #Fit with LL.2()
        res   <- c(i, nGerm, nTot, pMaxO, 3, coefs[1],  1,  coefs[2], tg1, tg2, GR1, GR2)
        resES <- c(i, NA, NA, NA,         3, coefES[1], NA, coefES[2], tg1es, tg2es, GR1es, GR2es)
        resLow<- c(i, NA, NA, NA,         3, coefL[1], NA,  coefL[2], tg1L, tg2L, GR1L, GR2L)
        resUp <- c(i, NA, NA, NA,         3, coefU[1], NA,  coefU[2], tg1U, tg2U, GR1U, GR2U)

        result <- rbind(result, res)
        resultES <- rbind(resultES, resES)
        resultLow <- rbind(resultLow, resLow)
        resultUp <- rbind(resultUp, resUp)

            coefModT <- c(4, coefs[1], NA, coefs[2], coefES[1], NA, coefES[2])
            coefMod <- rbind(coefMod, coefModT)

        } else{
          coefs <- coef(cureMod)
        #LL.3 was ok
        #mod = 5
        coefES <- summary(cureMod)$coef[,2]
        coefL <- coefs - 2*coefES
        coefU <- coefs + 2*coefES
        RSS <- sum(residuals(cureMod)^2)
        #plot(cureMod, main=i, axes=F)
        tg1o <- ED(cureMod, seq(0.10, 0.90, 0.1), type = "absolute", display=F)
        tg2o <- ED(cureMod, seq(0.10, 0.90, 0.1), type = "relative", display=F)
        tg1 <- tg1o[,1]; tg1es <- tg1o[,2]; tg1L <- tg1 - 2*tg1es; tg1U <- tg1 + 2*tg1es
        tg2 <- tg2o[,1]; tg2es <- tg2o[,2]; tg2L <- tg2 - 2*tg2es; tg2U <- tg2 + 2*tg2es
        GR1 <- 1/tg1
        GR2 <- 1/tg2
        GR1es <- sqrt( (tg1es^2) * ((-1/tg1)^2)^2 ) #Delta method
        GR2es <- sqrt( (tg2es^2) * ((-1/tg2)^2)^2 )
        GR1L <- GR1 - 2*GR1es; GR1U <- GR1 + 2*GR1es
        GR2L <- GR2 - 2*GR2es; GR2U <- GR2 + 2*GR2es
        tg1[is.na(tg1)] <- NA
        GR1[is.na(GR1)] <- 1E-6
        tg1es[is.na(tg1es)] <- NA; tg1L[is.na(tg1L)] <- mean(tLast); tg1U[is.na(tg1U)] <- Inf
        GR1es[is.na(GR1es)] <- 1E-6; GR1L[is.na(GR1L)] <- 1/mean(tLast); GR1U[is.na(GR1U)] <- 0

        #Fit with LL.3()
        res    <- c(i, nGerm, nTot, pMaxO, 4, coefs, tg1, tg2, GR1, GR2)
        resES  <- c(i, NA, NA, NA, 4, coefES, tg1es, tg2es, GR1es, GR2es)
        resLow <- c(i, NA, NA, NA, 4, coefL, tg1L, tg2L, GR1L, GR2L)
        resUp  <- c(i, NA, NA, NA, 4, coefU, tg1U, tg2U, GR1U, GR2U)
          coefModT <- c(5, coefs, coefES)
          coefMod <- rbind(coefMod, coefModT)

        result <- rbind(result, res);
        resultES <- rbind(resultES, resES)
        resultLow <- rbind(resultLow, resLow)
        resultUp <- rbind(resultUp, resUp)
        cat(paste("Group ", i, ": ", fct, ".3() was fitted", "\n", sep=""))


          } }
    }

  names(coefMod) <- c("Code", "b", "d", "e", "es_b", "es_d", "es_e")
  names(result) <- c("Group", "nGerm", "nTot", "PmaxObs", "Model", "b", "d", "e", paste("GTA", seq(10,90, by=10), sep=""),
                         paste("GTR", seq(10,90, by=10), sep=""),
                         paste("GRA", seq(10,90, by=10), sep=""),
                         paste("GRR", seq(10,90, by=10), sep=""))
  names(resultUp) <- names(resultLow) <- names(resultES) <- names(result)
  result <- data.frame(Group=levels(Group), result[,2:length(result[1,])])
  resultES <- data.frame(Group=levels(Group), resultES[,2:length(resultES[1,])])
  resultLow <- data.frame(Group=levels(Group), resultLow[,2:length(resultLow[1,])])
  resultUp <- data.frame(Group=levels(Group), resultUp[,2:length(resultUp[1,])])
  print("Process successfully finished")
  coefs <- cbind(result[,c(1, 5, 6:8)], resultES[,6:8])
  results <- list("Estimates"=result, "SE"=resultES, "Low"=resultLow, "Up"=resultUp)
  returnList <- list(coefficients = coefs, results = results, pMaxFin = pMaxFin,
                     coefs = coefMod)
  }
