drmSG <- function(formula, data, curveid, fct="LL", min=0.04, maxIni = 0.95,
                  g = c(10, 30, 50), type = "absolute") {
  # formula <- nSeeds ~ timeBef + timeAf
  # curveid <- rape$Dish
  # data <- rape; fct ="LL"
  # min=0.04; g = c(10, 30, 50)
  # Group <- as.factor(temp$Dish)
  # Return code. 0: no germinations; 1: negligible germinations
  # 2: complete germination at first inspection; 3: fit is impossible
  # 4: LL.2() and 5: LL.3()
  #Updated: 16/7/19
  probs <- g
  if(type != "absolute" & type != "relative") {
    cat(paste("ERROR: argument type can only be 'absolute or 'relative'", "\n", sep=""))
    stop() }

  fr <- model.frame(formula, data)
  timeBef <- model.matrix(fr, data)[,2]
  timeAf <- model.matrix(fr, data)[,3]
  nSeeds <- model.response(fr, "numeric")
  anName <- deparse(substitute(curveid))  # storing name for later use
  Group <- factor( subset(data, select = anName) [,1])
  DataC <- data.frame(Group, timeBef, timeAf, nSeeds)

  pFinal <- data.frame()
  coefMod <- data.frame()
  GRest <- data.frame()
  Test <- data.frame()
  GRest.se <- data.frame()
  Test.se <- data.frame()
  Gname <- c()

  if(fct == "LL"){fct1 <- LL.3(); fct2<- LL.2()
  }else if(fct == "W1") { fct1 <- W1.3(); fct2<- W1.2()
  }else if(fct == "W2") { fct1 <- W2.3(); fct2<- W2.2() }

  #Fitting group by group
  nLev <- length(levels(Group))
  for(i in 1:nLev){

      Gname[i] <- levels(Group)[i]
      dataTemp <- subset(DataC, Group==levels(Group)[i])
      nTot <- sum(dataTemp$nSeeds)
      nGerm <- nTot - sum( dataTemp[dataTemp$timeAf==Inf,]$nSeeds )

      if(type == "absolute"){
        sumSeeds <- nTot
        drmProbs <- probs/100
      } else {sumSeeds <- nGerm
        drmProbs <- probs}

      if(nGerm == 0){
      final <- pMaxFin(time=dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                       counts=dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                       nSeeds=nTot)
      }else{
      final <- pMaxFin(time=dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                       counts=dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                       nSeeds=nTot, se.fit = T)
      }
      pMax0 <- final$pMaxFin

      #Return pMax0
      pFinal[i, 1] <- Gname[i]; pFinal[i, 2] <- nGerm
      pFinal[i, 3] <- nTot; pFinal[i, 4] <- final$pMaxFin
      if(nGerm == 0){
      pFinal[i, 5] <- 0
      }else{
      pFinal[i, 5] <- final$se
      }
      colnames(pFinal) <- c("Group", "nGerm", "nTot", "pFinal", "SE")


      nFirst <- sum( dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$nSeeds )
      pFirst <- nFirst/nTot
      tFirst <- dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$timeAf
      tLast <- dataTemp[is.finite(dataTemp$timeAf)==F,]$timeBef

      if(nGerm == 0){
        #There are no germinations at all. Code = 0
        coefModT <- c(0, NA, NA, NA, NA, NA, NA)
        coefMod <- rbind(coefMod, coefModT)

        GReste <- rep(0, length(probs))
        GRestSE <- rep(0, length(probs))

        Teste <- rep(Inf, length(probs))
        TesteSE <- rep(NA, length(probs))

        GRest <- rbind(GRest, c(0, GReste) )
        GRest.se <- rbind(GRest.se, c(0, GRestSE) )
        Test <- rbind(Test, c(0, Teste ) )
        Test.se <- rbind(Test.se, c(0, TesteSE ) )
        message(paste("Group ", i, ": no germinations were observed", sep=""))
        next;
      }else if(pMax0 < min){
        #If minimum threshold of germination is not reached
        #I assume there is negligible germination (mod = 1)
        #min may be user-defined
        #Standard errors for GR are not calculated

        coefModT <- c(1, NA, NA, NA, NA, NA, NA)
        coefMod <- rbind(coefMod, coefModT)

        GRest.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                              dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                              probs = probs, nSeeds=sumSeeds, type=1,
                              rate=T, se.fit = F)
        Test.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                             dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                             probs = probs, nSeeds=sumSeeds, type=1,
                             rate=F)
        GRest <- rbind(GRest, c(1, GRest.1) )
        GRest.se <- rbind(GRest.se, c(1, rep(NA, length(probs))) )
        Test <- rbind(Test, c(1, Test.1) )
        Test.se <- rbind(Test.se, c(1, rep(NA, length(probs)) ) )

        message(paste("Group ", i, ": pMax lower than minimum threshold", sep=""))
        next;
        }else if(pFirst > maxIni){
          #Pmax at first inspection > 0.95. Cannot fit germination model
          #(mod = 2)

          coefModT <- c(2, NA, NA, NA, NA, NA, NA)
          coefMod <- rbind(coefMod, coefModT)

          GRest.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                                dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                                probs = probs, nSeeds=sumSeeds, type=1,
                                rate=T, se.fit = T)
          Test.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                               dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                               probs = probs, nSeeds=sumSeeds, type=1,
                               rate=F)

          GRest <- rbind(GRest, c(2, GRest.1$estimate))
          GRest.se <- rbind(GRest.se, c(2, GRest.1$se) )
          Test <- rbind(Test, c(2, Test.1) )
          Test.se <- rbind(Test.se, c(2, rep(NA, length(probs)) ) )

          message(paste("Group ", i, ": Germination at first inspection almost complete", sep=""))
          next;
          } else {
          #Try to fit a germination model with two or three-parameters
          cureMod <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
                              fct = fct1, type = "event"), silent=T)
          cureMod2 <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
                               fct = fct2, type = "event"), silent=T )
          #options(echo=T)
          }

      #Try whether, in spite of succesful fit, summary does not work
       if( class(cureMod) != "try-error") {
             p <- try( summary(cureMod), silent=T)
             if(class(p) == "try-error") class(cureMod) <- "try-error"
       }
       if( class(cureMod2) != "try-error") {
             p <- try( summary(cureMod2), silent=T)
             if(class(p) == "try-error") class(cureMod2) <- "try-error"
       }
      #Look at what model is OK
      #class(cureMod); class(cureMod2)
      if(class(cureMod) == "try-error" & class(cureMod2) == "try-error"){

        #No parametric fit was possible (mod = 3)
        message(paste("Group ", i, ": No parametric fit was possible", sep=""))

        coefModT <- c(3, NA, NA, NA, NA, NA, NA)
        coefMod <- rbind(coefMod, coefModT)

        GRest.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                              dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                              probs = probs, nSeeds=sumSeeds, type=1,
                              rate=T, se.fit = T)
        Test.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                             dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                             probs = probs, nSeeds=sumSeeds, type=1,
                             rate=F)
        GRest <- rbind(GRest, c(3, GRest.1$estimate) )
        GRest.se <- rbind(GRest.se, c(3, GRest.1$se) )
        Test <- rbind(Test, c(3, Test.1) )
        Test.se <- rbind(Test.se, c(3, rep(NA, length(probs)) ) )

        }else if( class(cureMod) == "try-error"  & class(cureMod2) == "drc" )
                {

        #LL.3 could not be fit, but LL.2 was ok
        #mod = 4
        message(paste("Group ", i, ": ", fct, ".3() could not be fitted. ", fct, ".2() was fitted instead", sep=""))

        #Fit LL.2()
        coefs <- coef(cureMod2)
        coefES <- summary(cureMod2)$coef[,2]
        coefModT <- c(4, coefs[1], 1, coefs[2], coefES[1], 0, coefES[2])
        coefMod <- rbind(coefMod, coefModT)

        GRest.1 <- GRate(cureMod2, respLev=c(drmProbs),  type = type, vcov. = vcov)
        Test.1 <-  GTime(cureMod2, respLev=c(drmProbs),  type = type, vcov. = vcov)
        GRest <- rbind(GRest, c(4, GRest.1$Estimate) )
        GRest.se <- rbind(GRest.se, c(4, GRest.1$SE) )
        Test <- rbind(Test, c(4, Test.1$Estimate ) )
        Test.se <- rbind(Test.se, c(4, Test.1$SE ) )

        }else if(class(cureMod) == "drc" & cureMod$coefficients[2] > 1)
                {
        #LL.3 could not be fit, but LL.2 was ok (same as above, different reason)
        #mod = 4
        message(paste("Group ", i, ": ", fct, ".3() could not be fitted. ", fct, ".2() was fitted instead", sep=""))

        #Fit LL.2()
        coefs <- coef(cureMod2)
        coefES <- summary(cureMod2)$coef[,2]
        coefModT <- c(4, coefs[1], 1, coefs[2], coefES[1], 0, coefES[2])
        coefMod <- rbind(coefMod, coefModT)

        GRest.1 <- GRate(cureMod2, respLev=c(drmProbs),  type = type, vcov. = vcov)
        Test.1 <-  GTime(cureMod2, respLev=c(drmProbs),  type = type, vcov. = vcov)
        GRest <- rbind(GRest, c(4, GRest.1$Estimate) )
        GRest.se <- rbind(GRest.se, c(4, GRest.1$SE) )
        Test <- rbind(Test, c(4, Test.1$Estimate ) )
        Test.se <- rbind(Test.se, c(4, Test.1$SE ) )

        }else{ coefs <- coef(cureMod)
      #LL.3 was ok
      #mod = 5
      #Fit LL.3()
      coefs <- coef(cureMod)
      coefES <- summary(cureMod)$coef[,2]

      coefModT <- c(5, coefs, coefES)
      coefMod <- rbind(coefMod, coefModT)

      GRest.1 <- GRate(cureMod, respLev=c(drmProbs),  type = type, vcov. = vcov)
      Test.1 <-  GTime(cureMod, respLev=c(drmProbs),  type = type, vcov. = vcov)
      GRest <- rbind(GRest, c(5, GRest.1$Estimate) )
      GRest.se <- rbind(GRest.se, c(5, GRest.1$SE) )
      Test <- rbind(Test, c(5, Test.1$Estimate ) )
      Test.se <- rbind(Test.se, c(5, Test.1$SE ) )

      message(paste("Group ", i, ": ", fct, ".3() was fitted", sep=""))
      }
    }

  names(coefMod) <- c("Code", "b", "d", "e", "es_b", "es_d", "es_e")
  names(GRest) <- c("Code", paste(probs, "%", sep="") )
  names(GRest.se) <- c("Code", paste(probs, "%", sep="") )
  names(Test) <- c("Code", paste(probs, "%", sep="") )
  names(Test.se) <- c("Code", paste(probs, "%", sep="") )

  coefMod <- data.frame(Group=Gname, coefMod, check.names=F)
  GRest <- data.frame(Group=Gname, GRest, check.names=F)
  GRest.se <- data.frame(Group=Gname, GRest.se, check.names=F)
  Test <- data.frame(Group=Gname, Test, check.names=F)
  Test.se <- data.frame(Group=Gname, Test.se, check.names=F)
  report <- data.frame(Group = Gname, Code = Test$Code)
  message("Process successfully finished")

  #coefMod <- coefMod[,-2]
  returnList <- list(pFinal = pFinal,
                     coefficients = coefMod[,-2], GRg = GRest[,-2],
                     GRg.se = GRest.se[,-2],
                     Tg = Test[,-2], Tg.se = Test.se[,-2], report = report)
  returnList
}

SGindices <- function(counts, dish, nViable, moniTimes, fct = "LL",
                      g = c(10, 30, 50), min = 0.04) {
    # Similar function, but working on the field-book data in wide form
    temp <- makeDrm(counts, treat = data.frame(tratt = dish), nViable, moniTimes)
    mod <- drmSG(count ~ timeBef + timeAf, curveid = Dish,
                       data = temp, fct = fct, g = g)

    return(mod)
}
# drm.sg <- function(formula, data, curveid, fct = "LL", min=0.04, maxIni = 0.95,
#                    g = c(10, 30, 50), type = "absolute") {
#   # formula <- nSeeds ~ timeBef + timeAf
#   # data(rape); data <- rape
#   # data <- dataset
#   # curveid <- data$Dish
#   # data <- dataset; fct ="LL"
#   # min=0.04; g = c(10, 30, 50)
#   # Group <- as.factor(curveid)
#   # Return code. 0: no germinations; 1: negligible germinations
#   # 2: complete germination at first inspection; 3: fit is impossible
#   # 4: LL.2() and 5: LL.3()
#   #Updated: 16/7/19
#   probs <- g
#   if(type != "absolute" & type != "relative") {
#     cat(paste("ERROR: argument type can only be 'absolute or 'relative'", "\n", sep=""))
#     stop() }
#
#   fr <- model.frame(formula, data)
#   timeBef <- model.matrix(fr, data)[,2]
#   timeAf <- model.matrix(fr, data)[,3]
#   nSeeds <- model.response(fr, "numeric")
#   anName <- deparse(substitute(curveid))  # storing name for later use
#   # Group <- factor( subset(data, select = Dish) [,1])
#   Group <- factor( subset(data, select = anName) [,1])
#   DataC <- data.frame(Group, timeBef, timeAf, nSeeds)
#
#   pFinal <- data.frame()
#   coefMod <- data.frame()
#   GRest <- data.frame()
#   Test <- data.frame()
#   GRest.se <- data.frame()
#   Test.se <- data.frame()
#   Gname <- c()
#
#   if(fct == "LL"){fct1 <- LL.3(); fct2<- LL.2()
#   }else if(fct == "W1") { fct1 <- W1.3(); fct2<- W1.2()
#   }else if(fct == "W2") { fct1 <- W2.3(); fct2<- W2.2() }
#
#   #Fitting group by group
#   nLev <- length(levels(Group))
#   for(i in 1:nLev){
#
#       Gname[i] <- levels(Group)[i]
#       dataTemp <- subset(DataC, Group==levels(Group)[i])
#       nTot <- sum(dataTemp$nSeeds)
#       nGerm <- nTot - sum( dataTemp[dataTemp$timeAf==Inf,]$nSeeds )
#       nUnGerm <- nTot - nGerm
#       pMax0 <- nGerm/nTot
#       pMax0se <- sqrt( nTot * pMax0 * (1 - pMax0) ) / sqrt(nTot)
#
#       # Individua i g level (relative or absolute)
#       if(type == "absolute"){
#         sumSeeds <- nTot
#         drmProbs <- probs/100
#       } else {sumSeeds <- nGerm
#         drmProbs <- probs}
#
#       pFinal[i, 1] <- Gname[i]; pFinal[i, 2] <- nGerm
#       pFinal[i, 3] <- nTot; pFinal[i, 4] <- pMax0
#       pFinal[i, 5] <- pMax0se
#
#       nFirst <- sum( dataTemp[dataTemp$timeBef == min(dataTemp$timeBef),]$nSeeds )
#       pFirst <- nFirst/nTot
#       tFirst <- dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$timeAf
#       tLast <- dataTemp[is.finite(dataTemp$timeAf)==F,]$timeBef
#
#       if(nGerm == 0){
#         #There are no germinations at all. Code = 0
#         coefModT <- c(0, NA, NA, 0, 0, NA, NA)
#         coefMod <- rbind(coefMod, coefModT)
#
#         GReste <- rep(0, length(probs))
#         GRestSE <- rep(0, length(probs))
#
#         Teste <- rep(Inf, length(probs))
#         TesteSE <- rep(NA, length(probs))
#
#         GRest <- rbind(GRest, c(0, GReste) )
#         GRest.se <- rbind(GRest.se, c(0, GRestSE) )
#         Test <- rbind(Test, c(0, Teste ) )
#         Test.se <- rbind(Test.se, c(0, TesteSE ) )
#         message(paste("Group ", i, ": no germinations were observed", sep=""))
#         next;
#       }else if(pMax0 < min){
#         # If minimum threshold of germination is not reached
#         # I assume there is negligible germination (mod = 1)
#         # min may be user-defined
#         # Standard errors for GR are not calculated
#
#         coefModT <- c(1, NA, NA, pMax0, pMax0se, NA, NA)
#         coefMod <- rbind(coefMod, coefModT)
#
#         GReste <- rep(0, length(probs))
#         GRestSE <- rep(0, length(probs))
#
#         Teste <- rep(Inf, length(probs))
#         TesteSE <- rep(NA, length(probs))
#
#         GRest <- rbind(GRest, c(0, GReste) )
#         GRest.se <- rbind(GRest.se, c(0, GRestSE) )
#         Test <- rbind(Test, c(0, Teste ) )
#         Test.se <- rbind(Test.se, c(0, TesteSE ) )
#
#         message(paste("Group ", i, ": pMax lower than minimum threshold", sep=""))
#         next;
#         } else {
#         #Try to fit a germination model with two or three-parameters
#         cureMod <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
#                             fct = fct1, type = "event"), silent=T)
#         cureMod2 <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
#                              fct = fct2, type = "event"), silent=T )
#         #options(echo=T)
#         }
#
#       #Try whether, in spite of succesful fit, summary does not work
#        if( class(cureMod) != "try-error") {
#              p <- try( summary(cureMod), silent=T)
#              if(class(p) == "try-error") class(cureMod) <- "try-error"
#        }
#        if( class(cureMod2) != "try-error") {
#              p <- try( summary(cureMod2), silent=T)
#              if(class(p) == "try-error") class(cureMod2) <- "try-error"
#        }
#       #Look at what model is OK
#       #class(cureMod); class(cureMod2)
#       if(class(cureMod) == "try-error" & class(cureMod2) == "try-error"){
#
#         #No parametric fit was possible (mod = 3)
#         message(paste("Group ", i, ": No parametric fit was possible", sep=""))
#
#         coefModT <- c(2, NA, NA, NA, NA, NA, NA)
#         coefMod <- rbind(coefMod, coefModT)
#         GRest <- rbind(GRest, c(2, rep(NA, length(probs))) )
#         GRest.se <- rbind(GRest.se, c(2, rep(NA, length(probs)) ))
#         Test <- rbind(Test, c(2, rep(NA, length(probs)) ))
#         Test.se <- rbind(Test.se, c(2, rep(NA, length(probs)) ) )
#
#         }else if( class(cureMod) == "try-error"  & class(cureMod2) == "drc" )
#                 {
#
#         #LL.3 could not be fit, but LL.2 was ok
#         #mod = 4
#         message(paste("Group ", i, ": ", fct, ".3() could not be fitted. ", fct, ".2() was fitted instead", sep=""))
#
#         coefs <- coef(cureMod2)
#         coefES <- summary(cureMod2)$coef[,2]
#         coefModT <- c(3, coefs[1], 1, coefs[2], coefES[1], 0, coefES[2])
#         coefMod <- rbind(coefMod, coefModT)
#
#         GRest.1 <- GRate(cureMod2, respLev=c(drmProbs),  type = type, vcov. = vcov)
#         Test.1 <-  GTime(cureMod2, respLev=c(drmProbs),  type = type, vcov. = vcov)
#         GRest <- rbind(GRest, c(3, GRest.1$Estimate) )
#         GRest.se <- rbind(GRest.se, c(3, GRest.1$SE) )
#         Test <- rbind(Test, c(3, Test.1$Estimate ) )
#         Test.se <- rbind(Test.se, c(3, Test.1$SE ) )
#
#         }else if(class(cureMod) == "drc" & cureMod$coefficients[2] > 1)
#                 {
#         #LL.3 could not be fit, but LL.2 was ok (same as above, different reason)
#         #mod = 4
#         message(paste("Group ", i, ": ", fct, ".3() could not be fitted. ", fct, ".2() was fitted instead", sep=""))
#
#         #Fit LL.2()
#         coefs <- coef(cureMod2)
#         coefES <- summary(cureMod2)$coef[,2]
#         coefModT <- c(4, coefs[1], 1, coefs[2], coefES[1], 0, coefES[2])
#         coefMod <- rbind(coefMod, coefModT)
#
#         GRest.1 <- GRate(cureMod2, respLev=c(drmProbs),  type = type, vcov. = vcov)
#         Test.1 <-  GTime(cureMod2, respLev=c(drmProbs),  type = type, vcov. = vcov)
#         GRest <- rbind(GRest, c(4, GRest.1$Estimate) )
#         GRest.se <- rbind(GRest.se, c(4, GRest.1$SE) )
#         Test <- rbind(Test, c(4, Test.1$Estimate ) )
#         Test.se <- rbind(Test.se, c(4, Test.1$SE ) )
#
#         }else{ coefs <- coef(cureMod)
#       #LL.3 was ok
#       #mod = 5
#       #Fit LL.3()
#       coefs <- coef(cureMod)
#       coefES <- summary(cureMod)$coef[,2]
#
#       coefModT <- c(5, coefs, coefES)
#       coefMod <- rbind(coefMod, coefModT)
#
#       GRest.1 <- GRate(cureMod, respLev=c(drmProbs),  type = type, vcov. = vcov)
#       Test.1 <-  GTime(cureMod, respLev=c(drmProbs),  type = type, vcov. = vcov)
#       GRest <- rbind(GRest, c(5, GRest.1$Estimate) )
#       GRest.se <- rbind(GRest.se, c(5, GRest.1$SE) )
#       Test <- rbind(Test, c(5, Test.1$Estimate ) )
#       Test.se <- rbind(Test.se, c(5, Test.1$SE ) )
#
#       message(paste("Group ", i, ": ", fct, ".3() was fitted", sep=""))
#       }
#     }
#
#   names(coefMod) <- c("Code", "b", "d", "e", "es_b", "es_d", "es_e")
#   names(GRest) <- c("Code", paste(probs, "%", sep="") )
#   names(GRest.se) <- c("Code", paste(probs, "%", sep="") )
#   names(Test) <- c("Code", paste(probs, "%", sep="") )
#   names(Test.se) <- c("Code", paste(probs, "%", sep="") )
#   colnames(pFinal) <- c("Group", "nGerm", "nTot", "pFinal", "SE")
#
#   coefMod <- data.frame(Group=Gname, coefMod, check.names=F)
#   GRest <- data.frame(Group=Gname, GRest, check.names=F)
#   GRest.se <- data.frame(Group=Gname, GRest.se, check.names=F)
#   Test <- data.frame(Group=Gname, Test, check.names=F)
#   Test.se <- data.frame(Group=Gname, Test.se, check.names=F)
#   report <- data.frame(Group = Gname, Code = Test$Code)
#   message("Process successfully finished")
#
#   #coefMod <- coefMod[,-2]
#   returnList <- list(pFinal = pFinal,
#                      coefficients = coefMod[,-2], GRg = GRest[,-2],
#                      GRg.se = GRest.se[,-2],
#                      Tg = Test[,-2], Tg.se = Test.se[,-2], report = report)
#   returnList
# }

# drm.sg.np <- function(formula, data, curveid, fct="LL", min=0.04, maxIni = 0.95,
#                   g = c(10, 30, 50), type = "absolute") {
#   # formula <- nSeeds ~ timeBef + timeAf
#   # curveid <- rape$Dish
#   # data <- rape; fct ="LL"
#   # min=0.04; g = c(10, 30, 50)
#   # Group <- as.factor(temp$Dish)
#   # Return code. 0: no germinations; 1: negligible germinations
#   # 2: complete germination at first inspection; 3: fit is impossible
#   # 4: LL.2() and 5: LL.3()
#   #Updated: 16/7/19
#   probs <- g
#   if(type != "absolute" & type != "relative") {
#     cat(paste("ERROR: argument type can only be 'absolute or 'relative'", "\n", sep=""))
#     stop() }
#
#   fr <- model.frame(formula, data)
#   timeBef <- model.matrix(fr, data)[,2]
#   timeAf <- model.matrix(fr, data)[,3]
#   nSeeds <- model.response(fr, "numeric")
#   anName <- deparse(substitute(curveid))  # storing name for later use
#   Group <- factor( subset(data, select = anName) [,1])
#   DataC <- data.frame(Group, timeBef, timeAf, nSeeds)
#
#   pFinal <- data.frame()
#   coefMod <- data.frame()
#   GRest <- data.frame()
#   Test <- data.frame()
#   GRest.se <- data.frame()
#   Test.se <- data.frame()
#   Gname <- c()
#
#   #Fitting group by group
#   nLev <- length(levels(Group))
#   for(i in 1:nLev){
#
#       Gname[i] <- levels(Group)[i]
#       dataTemp <- subset(DataC, Group==levels(Group)[i])
#       nTot <- sum(dataTemp$nSeeds)
#       nGerm <- nTot - sum( dataTemp[dataTemp$timeAf==Inf,]$nSeeds )
#       nUnGerm <- nTot - nGerm
#       pMax0 <- nGerm/nTot
#       pMax0se <- sqrt( nTot * pMax0 * (1 - pMax0) ) / sqrt(nTot)
#
#       # Individua i g level (relative or absolute)
#       if(type == "absolute"){
#         sumSeeds <- nTot
#         drmProbs <- probs/100
#       } else {sumSeeds <- nGerm
#         drmProbs <- probs}
#
#       pFinal[i, 1] <- Gname[i]; pFinal[i, 2] <- nGerm
#       pFinal[i, 3] <- nTot; pFinal[i, 4] <- pMax0
#       pFinal[i, 5] <- pMax0se
#
#       nFirst <- sum( dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$nSeeds )
#       pFirst <- nFirst/nTot
#       tFirst <- dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$timeAf
#       tLast <- dataTemp[is.finite(dataTemp$timeAf)==F,]$timeBef
#
#       if(nGerm == 0){
#         #There are no germinations at all. Code = 0
#         coefModT <- c(0, NA, NA, NA, NA, NA, NA)
#         coefMod <- rbind(coefMod, coefModT)
#
#         GReste <- rep(0, length(probs))
#         GRestSE <- rep(0, length(probs))
#
#         Teste <- rep(Inf, length(probs))
#         TesteSE <- rep(NA, length(probs))
#
#         GRest <- rbind(GRest, c(0, GReste) )
#         GRest.se <- rbind(GRest.se, c(0, GRestSE) )
#         Test <- rbind(Test, c(0, Teste ) )
#         Test.se <- rbind(Test.se, c(0, TesteSE ) )
#         message(paste("Group ", i, ": no germinations were observed", sep=""))
#         next;
#       }else if(pMax0 < min){
#         # If minimum threshold of germination is not reached
#         # I assume there is negligible germination (mod = 1)
#         # min may be user-defined
#         # Standard errors for GR are not calculated
#
#         coefModT <- c(1, NA, NA, NA, NA, NA, NA)
#         coefMod <- rbind(coefMod, coefModT)
#
#         GRest.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
#                               dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
#                               probs = probs, nSeeds=sumSeeds, type=1,
#                               rate=T, se.fit = F)
#         Test.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
#                              dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
#                              probs = probs, nSeeds=sumSeeds, type=1,
#                              rate=F)
#         GRest <- rbind(GRest, c(1, GRest.1) )
#         GRest.se <- rbind(GRest.se, c(1, rep(NA, length(probs))) )
#         Test <- rbind(Test, c(1, Test.1) )
#         Test.se <- rbind(Test.se, c(1, rep(NA, length(probs)) ) )
#
#         message(paste("Group ", i, ": pMax lower than minimum threshold", sep=""))
#         next;
#         } else {
#           #Pmax at first inspection > 0.95. Cannot fit germination model
#           #(mod = 2)
#
#           coefModT <- c(2, NA, NA, NA, NA, NA, NA)
#           coefMod <- rbind(coefMod, coefModT)
#
#           GRest.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
#                                 dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
#                                 probs = probs, nSeeds=sumSeeds, type=1,
#                                 rate=T, se.fit = T)
#           Test.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
#                                dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
#                                probs = probs, nSeeds=sumSeeds, type=1,
#                                rate=F)
#
#           GRest <- rbind(GRest, c(2, GRest.1$estimate))
#           GRest.se <- rbind(GRest.se, c(2, GRest.1$se) )
#           Test <- rbind(Test, c(2, Test.1) )
#           Test.se <- rbind(Test.se, c(2, rep(NA, length(probs)) ) )
#
#           message(paste("Group ", i, ": Germination at first inspection almost complete", sep=""))
#
#         }
#   } # End for
#   names(coefMod) <- c("Code", "b", "d", "e", "es_b", "es_d", "es_e")
#   names(GRest) <- c("Code", paste(probs, "%", sep="") )
#   names(GRest.se) <- c("Code", paste(probs, "%", sep="") )
#   names(Test) <- c("Code", paste(probs, "%", sep="") )
#   names(Test.se) <- c("Code", paste(probs, "%", sep="") )
#
#   coefMod <- data.frame(Group=Gname, coefMod, check.names=F)
#   GRest <- data.frame(Group=Gname, GRest, check.names=F)
#   GRest.se <- data.frame(Group=Gname, GRest.se, check.names=F)
#   Test <- data.frame(Group=Gname, Test, check.names=F)
#   Test.se <- data.frame(Group=Gname, Test.se, check.names=F)
#   report <- data.frame(Group = Gname, Code = Test$Code)
#   message("Process successfully finished")
#
#   #coefMod <- coefMod[,-2]
#   returnList <- list(pFinal = pFinal,
#                      coefficients = coefMod[,-2], GRg = GRest[,-2],
#                      GRg.se = GRest.se[,-2],
#                      Tg = Test[,-2], Tg.se = Test.se[,-2], report = report)
#   returnList
# }
#


drm.group <- function(formula, curveid, pmodels, weights, data = NULL,
                      subset, fct, type = c("continuous"),
                      bcVal = NULL, bcAdd = 0, start, naaction = na.omit,
                      robust = "mean", logDose = NULL, control = drmc(),
                      lowerl = NULL, upperl = NULL, separate = FALSE,
                      pshifts = NULL, varcov = NULL) {
  # formula <- GR ~ TempFix
  # data <- dataset; curveid <- Comb3
  # fr <- model.frame(formula, data)
  # timeBef <- model.matrix(fr, data)[,2]
  # timeAf <- model.matrix(fr, data)[,3]
  # nSeeds <- model.response(fr, "numeric")
  anName <- deparse(substitute(curveid)) # storing name for later use
  if(is.null(anName)){
    print("Curveid argument is not given. This function works only on grouped data")
    stop() }

  fitGroup <- function(df) {
    fitTry <- try(drm(formula, pmodels, weights, subset, fct, data = df,
                       type, bcVal, start, na.action, robust, logDose,
                       lowerl, upperl, separate,
                       pshifts, varcov), silent = T)
    if(class(fitTry) == "try-error"){
      return(class(fitTry))
    } else {
      return(fit) }
  }
  resultList <- dlply(data, c(anName), fitGroup(df), .progress = "text" )
  return(invisible(resultList))
}

nlmod <- drm.group(GR ~ TempFix, curveid = Comb3, data = dataset,
             fct = GRT.Exb())

drm.sg <- function(formula, data, curveid = NULL, fct = "NP", min = 0.04, maxIni = 0.95,
                    g = c(10, 30, 50), type = "absolute", se.q = F,
                    B = 500, seed = 1234) {
  # formula <- nSeeds ~ timeBef + timeAf
  # data <- dataset
  probs <- g
  if(type != "absolute" & type != "relative") {
    cat(paste("ERROR: argument type can only be 'absolute or 'relative'", "\n", sep=""))
    stop() }
  fr <- model.frame(formula, data)
  timeBef <- model.matrix(fr, data)[,2]
  timeAf <- model.matrix(fr, data)[,3]
  nSeeds <- model.response(fr, "numeric")
  anName <- deparse(substitute(curveid)) # storing name for later use
  if(is.null(anName)){
    DataC <- data.frame(Group = 1, timeBef, timeAf, nSeeds)
  } else {
    # anName <- deparse(substitute(curveid))
    Group <- factor( subset(data, select = anName) [,1])
    DataC <- data.frame(Group, timeBef, timeAf, nSeeds)
  }
  if(fct == "NP") {
    # Non parametric fit
    result <- dlply(DataC, c("Group"),
      function(df) npFit(formula,
                       data = df, g = probs),
      .progress = "text" )
  } else if (fct == "KM") {
    # KM estimators
    result <- dlply(DataC, c("Group"),
      function(df) kmFit(formula,
                       data = df, g = probs, se.q = se.q,
                       B = 500, seed = 1234),
      .progress = "text" )

  } else {
    # parametric fit
   result <- dlply(DataC, c("Group"),
      function(df) pFit(formula, fct = fct,
                       data = df, g = probs),
      .progress = "text" )
  }

# result
return(invisible(result))

}

npFit <- function(formula, data, min=0.04, maxIni = 0.95,
                       g = c(10, 30, 50), type = "absolute") {
  # Function to fit non-parametric survival models to one group
  # Updated: 12/5/2020
  probs <- g
  if(type != "absolute" & type != "relative") {
    cat(paste("ERROR: argument type can only be 'absolute or 'relative'", "\n", sep=""))
    stop() }

  fr <- model.frame(formula, data)
  timeBef <- model.matrix(fr, data)[,2]
  timeAf <- model.matrix(fr, data)[,3]
  nSeeds <- model.response(fr, "numeric")
  anName <- deparse(substitute(curveid))  # storing name for later use
  dataTemp <- data.frame(timeBef, timeAf, nSeeds)
  nTot <- sum(dataTemp$nSeeds)
  nGerm <- nTot - sum( dataTemp[dataTemp$timeAf==Inf,]$nSeeds )
  nUnGerm <- nTot - nGerm
  sumSeeds <- nTot
  pMax0 <- nGerm/nTot
  pMax0se <- sqrt( nTot * pMax0 * (1 - pMax0) ) / sqrt(nTot)

  nFirst <- sum( dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$nSeeds )
  pFirst <- nFirst/nTot
  tFirst <- dataTemp[dataTemp$timeBef==min(dataTemp$timeBef),]$timeAf
  tLast <- dataTemp[is.finite(dataTemp$timeAf)==F,]$timeBef

  if(nGerm == 0){
    #There are no germinations at all. Code = 0
    GReste <- rep(0, length(probs))
    GRestSE <- rep(0, length(probs))
    Teste <- rep(Inf, length(probs))
  }else if(pMax0 < min){
    # If minimum threshold of germination is not reached
    # I assume there is negligible germination (mod = 1)
    # min may be user-defined
    # Standard errors for GR are not calculated
    GReste <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                          dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                          probs = probs, nSeeds=sumSeeds, type=1,
                          rate=T, se.fit = F)
    Teste <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                         dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                         probs = probs, nSeeds=sumSeeds, type=1,
                         rate=F)
    GRestSE <- rep(NA, length(probs))
  } else {
    #Pmax at first inspection > 0.95. Cannot fit germination model
    #(mod = 2)
    GRest.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                          dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                          probs = probs, nSeeds=sumSeeds, type=1,
                          rate=T, se.fit = T)
    Test.1 <- quantileSG(dataTemp[dataTemp$timeAf!=Inf,]$timeAf,
                         dataTemp[dataTemp$timeAf!=Inf,]$nSeeds,
                         probs = probs, nSeeds=sumSeeds, type=1,
                         rate=F)

    GReste <- GRest.1$estimate
    GRestSE <- GRest.1$se
    Teste <- Test.1
   }
  return <- c(pMax = pMax0, pMax.ES = pMax0se, GR = GReste, GR.ES = GRestSE,
              gT = Teste)
  return
}

kmFit <- function(formula, data, min=0.04, maxIni = 0.95,
                       g = c(10, 30, 50), type = "absolute",
                  B = 1000, seed = 1234, se.q = se.q) {
  # Function to fit non-parametric survival models to one group
  # Updated: 12/5/2020
  # formula <- nSeeds ~ timeBef + timeAf
  # data <-  dataset
  probs <- g/100
  if(type != "absolute" & type != "relative") {
    cat(paste("ERROR: argument type can only be 'absolute or 'relative'", "\n", sep=""))
    stop() }

  fr <- model.frame(formula, data)
  timeBef <- model.matrix(fr, data)[,2]
  timeAf <- model.matrix(fr, data)[,3]
  nSeeds <- model.response(fr, "numeric")

  dataTemp <- data.frame(timeBef, timeAf, nSeeds)
  nTot <- sum(dataTemp$nSeeds)
  nGerm <- nTot - sum( dataTemp[dataTemp$timeAf==Inf,]$nSeeds )
  nUnGerm <- nTot - nGerm
  sumSeeds <- nTot
  pMax0 <- nGerm/nTot
  pMax0se <- sqrt( nTot * pMax0 * (1 - pMax0) ) / sqrt(nTot)

  if(nGerm == 0){
    #There are no germinations at all. Code = 0
    GReste <- rep(0, length(probs))
    GRestSE <- rep(0, length(probs))
    Teste <- rep(Inf, length(probs))
    TestSE <- rep(NA, length(probs))

  } else {
    # KM estimators
    dataTemp <- makeSurv2(nSeeds, data = dataTemp)
    survObj <- survival::survfit(survival::Surv(timeBef, timeAf, type = "interval2") ~ 1,
          data = dataTemp, conf.type = "none")
    Test <- c()
    TestSE <- c()
    for(j in 1:length(probs)){
      q <- 1 - probs[j]
      Test[j] <- as.numeric(unname(quantile(survObj, probs = 1 - q))[[1]])
      quant_est <- c()
      if(is.na(Test[j]) == F & se.q == T){
      for (i in 1:B) {
          btsp <- sample(c(1:length(dataTemp$timeBef)), replace = TRUE)
          tmp_timeBef <- dataTemp$timeBef[btsp]
          tmp_timeAf <- dataTemp$timeAf[btsp]
          tmp_df <- data.frame(tmp_timeBef, tmp_timeAf)
          tmp_df <- tmp_df[order(tmp_df$tmp_timeBef), ]
          fit1 <- survival::survfit(survival::Surv(tmp_timeBef, tmp_timeAf, type = "interval2") ~ 1,
                                    conf.type = "none", data = tmp_df)
          Finv <- as.numeric(unname(quantile(fit1, probs = 1 - q))[[1]])
          if (is.na(Finv)) Finv <- max(tmp_timeBef)
          quant_est[i] <- Finv
      }
      TestSE[j] <- sd(quant_est)
      } else { TestSE[j] <- NA}
    }
  GRest <- 1/Test
  GRestSE <- (TestSE * 1/Test^2)
  }
  names(Test) <- paste("T:", g, "%", sep = "")
  names(TestSE) <- paste("T.ES:", g, "%", sep = "")
  names(GRest) <- paste("GR:", g, "%", sep = "")
  names(GRestSE) <- paste("GR.ES:", g, "%", sep = "")
  return <- c(pMax = pMax0, pMax.ES = pMax0se, GRest,
              GRestSE, Test, TestSE)
  return
}

pFit <- function(formula, data, fct = "LL", min=0.04, maxIni = 0.95,
                   g = c(10, 30, 50), type = "absolute") {
  # Function to fit parametric survival models to one group
  # If the curve is degenerated, returns an error message
  # Updated: 12/5/2020
  probs <- g
  if(type != "absolute" & type != "relative") {
    cat(paste("ERROR: argument type can only be 'absolute or 'relative'", "\n", sep=""))
    stop() }

  fr <- model.frame(formula, data)
  timeBef <- model.matrix(fr, data)[,2]
  timeAf <- model.matrix(fr, data)[,3]
  nSeeds <- model.response(fr, "numeric")
  dataTemp <- data.frame(timeBef, timeAf, nSeeds)
  nTot <- sum(dataTemp$nSeeds)

  # Individua i g level (relative or absolute)
  if(type == "absolute"){
    sumSeeds <- nTot
    drmProbs <- probs/100
    } else {sumSeeds <- nGerm
    drmProbs <- probs}

  # Select the distribution
  if(fct == "LL"){fct1 <- LL.dist(); fct2 <- LL.2()
  } else if(fct == "W1") { fct1 <- W1.3(); fct2<- W1.2()
  } else if(fct == "W2") { fct1 <- W2.3(); fct2<- W2.2() }

  #sink("error_log.txt")
  # Fitting three parameters
  mod <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
                            fct = fct1, type = "event",
                  upperl = c(NA, 1, NA)), silent=T)
  if( class(mod) != "try-error") {
      p <- try( summary(mod), silent=T)
      if(class(p) == "try-error") class(mod) <- "try-error"
  }

  # Fitting two parameters
  mod2 <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
                               fct = fct2, type = "event"),
                    silent=T )
  if( class(mod2) != "try-error") {
      p <- try( summary(mod2), silent=T)
      if(class(p) == "try-error") class(mod2) <- "try-error"
  }

  #Look at what model is OK
  if(class(mod) == "try-error" & class(mod2) == "try-error"){
      #No parametric fit was possible (mod = 3)
      return <- "No parametric fit was possible"
  } else if( class(mod) == "try-error"  & class(mod2) == "drc" ){
      #three parameters could not be fit, but LL.2 was ok
      coefs <- c(coef(mod2)[1], 1, coef(mod2)[2])
      coefES <- c(summary(mod2)$coef[1,2], NA, summary(mod2)$coef[2,2])
      Test.1 <-  ED(mod2, respLev=c(drmProbs), type = type, vcov. = vcov, display = F)
      Test <- Test.1[,1]
      Test.se <- Test.1[,2]
      GRest <- 1/Test
      GRest.se <- (Test.se * 1/Test^2)
      return <- c(coefs, coefES, GRest, GRest.se, Test, Test.se)
  } else {
      #three parameters was ok
      coefs <- coef(mod)
      coefES <- summary(mod)$coef[,2]
      Test.1 <-  ED(mod, respLev=c(drmProbs), type = type, vcov. = vcov, display = F)
      Test <- Test.1[,1]
      Test.se <- Test.1[,2]
      GRest <- 1/Test
      GRest.se <- (Test.se * 1/Test^2)
      return <- c(coefs, coefES, GRest, GRest.se, Test, Test.se)
  }
  if(length(return) > 1) {
    names(return) <- c("b", "d", "e", "SE.b", "SE.d", "SE.e",
                       paste("GR:", g, "%", sep = ""),
                       paste("SE.GR:", g, "%", sep = ""),
                       paste("T:", g, "%", sep = ""),
                       paste("SE.T:", g, "%", sep = ""))
  }
  return
}

