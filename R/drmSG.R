drmSG <- function(formula, data, curveid, fct="LL", min=0.1, probs = c(10, 30, 50)) {
  #formula <- nSeeds ~ timeBef + timeAf
  #curveid <- Dish
  #data <- rape; fct ="LL"
  #min=0.1; probs = c(10, 30, 50)
  #Group <- as.factor(rape$Dish)

  fr <- model.frame(formula, data)
  timeBef <- model.matrix(fr, data)[,2]
  timeAf <- model.matrix(fr, data)[,3]
  nSeeds <- model.response(fr, "numeric")

  anName <- deparse(substitute(curveid))  # storing name for later use
  Group <- as.factor( subset(data, select = anName) [,1])

  DataC <- data.frame(Group, timeBef, timeAf, nSeeds)

  pFinal <- data.frame()
  coefMod <- data.frame()
  GRest <- data.frame()
  Test <- data.frame()
  GRest.se <- data.frame()
  Gname <- c()

  if(fct == "LL"){fct1 <- LL.3(); fct2<- LL.2()
  }else if(fct == "W1") { fct1 <- W1.3(); fct2<- W1.2()}

  #Fitting group by group
  nLev <- length(levels(Group))
    for(i in 1:nLev){

      Gname[i] <- levels(Group)[i]
      dataTemp <- subset(DataC, Group==levels(Group)[i])
      nTot <- sum(dataTemp$nSeeds)
      nGerm <- nTot - sum( dataTemp[dataTemp$timeAf==Inf,]$nSeeds )
      #pMaxO <- nGerm/nTot
      final <- pMaxFin(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                         dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                       nTot)
      pMaxO <- final$pMaxFin

      #Return pMx0 #################################
      pFinal[i, 1] <- Gname[i]; pFinal[i, 2] <- nGerm
      pFinal[i, 3] <- nTot; pFinal[i, 4] <- final$pMaxFin
      pFinal[i, 5] <- final$se
      colnames(pFinal) <- c("Group", "nGerm", "nTot", "pFinal", "SE")

      nFirst <- sum( dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$nSeeds )
      pFirst <- nFirst/nTot
      tFirst <- dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$timeAf
      tLast <- dataTemp[is.finite(dataTemp$timeAf)==F,]$timeBef

      if(pMaxO < min){
        #If minimum threshold of germination is not reached
        #I assume there is negligible germination (mod = 1)
        #min may be user-defined

        coefModT <- c(1, NA, NA, NA, NA, NA, NA)
        coefMod <- rbind(coefMod, coefModT)

        GRest.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                              dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                              probs = probs, nSeeds=nTot, type=1,
                              rate=T, se.fit = T)
        Test.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                             dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                             probs = probs, nSeeds=nTot, type=1,
                             rate=F)
        GRest <- rbind(GRest, c(1, GRest.1$estimate) )
        GRest.se <- rbind(GRest.se, c(1, GRest.1$se) )
        Test <- rbind(Test, c(1, Test.1) )

        cat(paste("Group ", i, ": pMax lower than minimum threshold", "\n", sep=""))
        next;
        }else if(pFirst > 0.95){
          #Pmax at first inspection > 0.95. Cannot fit germination model
          #(mod = 2)

          coefModT <- c(2, NA, NA, NA, NA, NA, NA)
          coefMod <- rbind(coefMod, coefModT)

          GRest.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                                dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                                probs = probs, nSeeds=nTot, type=1,
                                rate=T, se.fit = T)
          Test.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                               dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                               probs = probs, nSeeds=nTot, type=1,
                               rate=F)

          GRest <- rbind(GRest, c(2, GRest.1$estimate))
          GRest.se <- rbind(GRest.se, c(2, GRest.1$se) )
          Test <- rbind(Test, c(1, Test.1) )

          cat(paste("Group ", i, ": Germination at first inspection almost complete", "\n", sep=""))
          next;
          } else {
          #Fit germination model: LL.2 and LL.3
          #A germination medel can be fit. Try LL2 and LL3
          #options(echo=F)
          #fct1 <- paste(fct, ".3()"); fct2 <- paste(fct, ".2()")
          #if(fct == "LL") fct1 <- LL.3(); fct2<- LL.2()

          cureMod <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
                              fct = fct1, type = "event",
                              upperl = c(NA, 1, NA)), silent=T)
          cureMod2 <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
                               fct = fct2, type = "event"), silent=T )
          #options(echo=T)
          }

      #Look at what model is OK
      #class(cureMod); class(cureMod2)
      if(class(cureMod) == "try-error" & class(cureMod2) == "try-error"){

        #No parameteric fit was possible (mod = 3). To be meditated
        cat(paste("Group ", i, ": No parametric fit was possible", "\n", sep=""))

        coefModT <- c(3, NA, NA, NA, NA, NA, NA)
        coefMod <- rbind(coefMod, coefModT)

        GRest.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                              dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                              probs = probs, nSeeds=nTot, type=1,
                              rate=T, se.fit = T)
        Test.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                             dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                             probs = probs, nSeeds=nTot, type=1,
                             rate=F)
        GRest <- rbind(GRest, c(3, GRest.1$estimate) )
        GRest.se <- rbind(GRest.se, c(3, GRest.1$se) )
        Test <- rbind(Test, c(3, Test.1) )


        }else if(class(cureMod) == "try-error" & class(cureMod2) == "drc"){

        #LL.3 could not be fit, but LL.2 was ok
        #mod = 4
        cat(paste("Group ", i, ": ", fct, ".3() could not be fit. ", fct, ".2() is fit instead", "\n", sep=""))

        #Fit LL.2()
        coefs <- coef(cureMod2)
        coefES <- summary(cureMod2)$coef[,2]
        coefModT <- c(4, coefs[1], 1, coefs[2], coefES[1], 0, coefES[2])
        coefMod <- rbind(coefMod, coefModT)

        GRest.1 <- GR(cureMod2, probs/100, type = "absolute")
        Test.1 <-  ED(cureMod2, probs/100, type = "absolute", display=F)
        GRest <- rbind(GRest, c(4, as.numeric(GRest.1$estimate)) )
        GRest.se <- rbind(GRest.se, c(4, as.numeric(GRest.1$se) ))
        Test <- rbind(Test, c(4, as.numeric(Test.1[,1]) ))
        }else{ coefs <- coef(cureMod)
      #LL.3 was ok
      #mod = 5
      #Fit LL.3()
      coefs <- coef(cureMod)
      coefES <- summary(cureMod)$coef[,2]

      coefModT <- c(5, coefs, coefES)
      coefMod <- rbind(coefMod, coefModT)

      GRest.1 <- GR(cureMod, probs/100, type = "absolute")
      Test.1 <-  ED(cureMod, probs/100, type = "absolute", display=F)
      GRest <- rbind(GRest, c(5, GRest.1$estimate) )
      GRest.se <- rbind(GRest.se, c(5, GRest.1$se) )
      Test <- rbind(Test, c(5, as.numeric(Test.1[,1]) ) )

      cat(paste("Group ", i, ": ", fct, ".3() was fitted", "\n", sep=""))
      }
    }

  names(coefMod) <- c("Code", "b", "d", "e", "es_b", "es_d", "es_e")
  names(GRest) <- c("Code", paste(probs, "%", sep="") )
  names(GRest.se) <- c("Code", paste(probs, "%", sep="") )
  names(Test) <- c("Code", paste(probs, "%", sep="") )

  coefMod <- data.frame(Group=Gname, coefMod, check.names=F)
  GRest <- data.frame(Group=Gname, GRest, check.names=F)
  GRest.se <- data.frame(Group=Gname, GRest.se, check.names=F)
  Test <- data.frame(Group=Gname, Test, check.names=F)
  report <- data.frame(Group = Gname, Code = Test$Code)
  print("Process successfully finished")

  returnList <- list(pFinal = pFinal,
                     coefficients = coefMod, GRest = GRest,
                     GRest.se = GRest.se,
                     Test = Test, report = report)
  returnList
}

GR <- function(mod, respLev, type="absolute"){
  GT <- ED(mod, respLev=respLev, type=type, display=F)
  GTval <- 1/GT[,1]
  GTes <- (GT[,2] * 1/GT[,1]^2)
  list(estimate=GTval, se=GTes)
}

