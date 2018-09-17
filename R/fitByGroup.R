fitByGroup <- function(Group, timeBef, timeAf, nSeeds) {
  #Group <- Dish
  Group <- as.factor(Group)
  DataC <- data.frame(Group, timeBef, timeAf, nSeeds)
  head(DataC)
  result <- data.frame()
  nLev <- length(levels(Group))
    for(i in 1:nLev){
      print(i)
      dataTemp <- subset(DataC, Group==levels(Group)[i])
      #if(dataTemp[dataTemp$timeAf==Inf,]$nSeeds == 0)
      nTot <- sum(dataTemp$nSeeds)
      nGerm <- nTot - dataTemp[dataTemp$timeAf==Inf,]$nSeeds
      pMaxO <- nGerm/nTot
      nFirst <- dataTemp[dataTemp$timeBef==0,]$nSeeds
      pFirst <- nFirst/nTot
      tFirst <- dataTemp[dataTemp$timeBef==0,]$timeAf

      if(pMaxO < 0.1){
        #Pmax < 0.1
        res <- c(i, nGerm, nTot, pMaxO, NA, NA, NA, Inf, rep(Inf, 18), rep(0, 18))

        result <- rbind(result, res)
        next;}else{if(pFirst > 0.95){

          #Pmax at first inspection > 0.95
          res <- c(i, nGerm, nTot, pMaxO, NA, NA, NA, tFirst, rep(tFirst, 18), rep(1/tFirst, 18))

          result <- rbind(result, res)
          next; } else{
            cureMod <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
                          fct = LL.3(), type = "event",
                          upperl = c(NA, 1, NA)), silent=T )
            cureMod2 <- try( drm(nSeeds ~ timeBef + timeAf, data = dataTemp,
                          fct = LL.2(), type = "event"), silent=T )
            } }
      class(cureMod); class(cureMod2)
      if(class(cureMod) == "try-error" & class(cureMod2) == "try-error"){

        #No fit was possible
        res <- c(i, nGerm, nTot, pMaxO, NA, NA, NA, NA, rep(NA, 18), rep(NA, 18))

        result <- rbind(result, res)
        } else {if(class(cureMod) == "try-error" & class(cureMod2) == "drc"){
          coefs <- coef(cureMod2)
          RSS <- sum(residuals(cureMod2)^2)
          #plot(cureMod, main=i, axes=F)
          tg1 <- ED(cureMod2, seq(0.10, 0.90, 0.1), type = "absolute", display=F)
          tg2 <- ED(cureMod2, seq(0.10, 0.90, 0.1), type = "relative", display=F)
          GR1 <- 1/tg1
          GR2 <- 1/tg2
          tg1[is.na(tg1)] <- Inf
          GR1[is.na(GR1)] <- 0
          #GR1[,1]

          #Fit with LL.2()
          res <- c(i, nGerm, nTot, pMaxO, 2, coefs[1], 1, coefs[2], tg1[,1], tg2[,1], GR1[,1], GR2[,1])

          result <- rbind(result, res)
          #print("caso2")
          } else{ coefs <- coef(cureMod)
          RSS <- sum(residuals(cureMod)^2)
          #plot(cureMod, main=i, axes=F)
          tg1 <- ED(cureMod, seq(0.10, 0.90, 0.1), type = "absolute", display=F)
          tg2 <- ED(cureMod, seq(0.10, 0.90, 0.1), type = "relative", display=F)
          GR1 <- 1/tg1
          GR2 <- 1/tg2
          tg1[is.na(tg1)] <- Inf
          GR1[is.na(GR1)] <- 0
          GR1[,1]

          #Fit with LL.3()
          res <- c(i, nGerm, nTot, pMaxO, 3, coefs, tg1[,1], tg2[,1], GR1[,1], GR2[,1])

          result <- rbind(result, res); #print("caso3")
          } }
    }

  names(result) <- c("Group", "nGerm", "nTot", "PmaxObs", "Model", "b", "d", "e", paste("GTA", seq(10,90, by=10), sep=""),
                         paste("GTR", seq(10,90, by=10), sep=""),
                         paste("GRA", seq(10,90, by=10), sep=""),
                         paste("GRR", seq(10,90, by=10), sep=""))
  result <- data.frame(Group=levels(Group), result[,2:length(result[1,])])
  result
}
