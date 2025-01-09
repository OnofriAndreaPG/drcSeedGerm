makeDrm <- function(counts, treat, nViable, moniTimes) {
  # this function reshapes a common recording sheet for germination assays
  # into the kind of dataset required by the function drm() in the drc package

  # counts is a dataframe listing for each Petri dish (rows) the number of germinated seeds
  # at each assessment time (culumns);
  # treat is a dataframe listing for each Petri dish (rows) the levels of each treatment (culumns)
  # number of rows in dataset and treat.
  # moniTimes = vector of monitoring times. Same length as number of columns in dataset
  # nViable is a vector with the number of viable seeds per each Petri dish

  # As for 15/12/2022 it has been superseded and replaced by the
  # melt_te function in drcte. Left for compatibility reasons
  msgText <- paste("'makeDrm' is deprecated.", "\n",
                   "Use 'melt_te' instead.", "\n",
                   "See ?melt_te for help on the new function",
                   sep = "")
  .Deprecated("melt_te", package = "drcte", msgText)

  # Creating objects to store the results
  dati <- as.data.frame(counts)
  tratt <- treat
  tempi <- c(moniTimes, Inf)

  dati[is.na(dati)] <- 0
  dati2 <- t(apply(dati, 1, cumsum))

  numPetri <- length(dati[,1]); numTimes <- length(tempi); numTesi<-length(tratt[1,])
  numRecords <- numPetri*numTimes
  nSeeds <- numeric(numRecords); nCum <- numeric(numRecords); Prop <- numeric(numRecords); Dish <- numeric(numRecords)
  nGerm <- numeric(numRecords); nGermPetri <- numeric(numRecords)
  nGermPetri <- apply(dati, 1, sum)
  tempi2 <- c(0, tempi[1:length(tempi)-1])
  group <- data.frame()
  final.date <- max(tempi2)
  Petri <- rep(c(1:numPetri), each=numTimes)
  timeAf <- rep(tempi, numPetri)
  timeBef <- rep(tempi2, numPetri)

  cont <- 1
  for (j in 1:numPetri){ #Per ogni riga - Petri dish
    for(i in 1:(numTimes-1)) { #Per ogni campionamento
      Dish[cont] <- j
      nSeeds[cont] <- dati[j,i]
      nGerm[cont] <- nGermPetri[j]
      nCum[cont] <- dati2[j,i]
      Prop[cont] <- dati2[j,i]/nViable[j]
      cont= cont+1
    }
    Dish[cont] <- j
    nSeeds[cont] <- nViable[j] - nGermPetri[j]
    nGerm[cont] <- nGermPetri[j]
    nCum[cont] <- NA
    Prop[cont] <- NA
    cont=cont+1
  }
  group <- tratt[rep(row.names(tratt), rep(numTimes, length(tratt[,1]))), 1:length(tratt[1,])]
  datiFin <- data.frame(group, Dish=Dish, timeBef=timeBef, timeAf=timeAf, count=nSeeds, nCum=nCum, propCum=Prop)
  return(datiFin)
}

makeDrm2 <- function(counts, treat, nViable, moniTimes, Dish,  cumulative=TRUE){
  #This function transform a dataset as cumulative germinations in a dataset for
  #use with DRM type = "event. Dish is a factor that identifies seeds
  #Group is a data.frame
  count <- counts; treatGroups <- treat; moniTime <- moniTimes
  temp <- data.frame()
  result <- data.frame()
  Dish <- factor(Dish)
  temp <- data.frame(moniTime, count, nViable, Dish, treatGroups)
  for(i in 1:length(levels(factor(temp$Dish)))){
  temp2 <- temp[temp$Dish==levels(factor(temp$Dish))[i],]
  timeBef <- c(0, temp2$moniTime)
  timeAf <- c(temp2$moniTime, Inf)
  nViable <- c(temp2$nViable, temp2$nViable[1])
  if(cumulative==T){
  #Decumulate if necessary
    nSeeds <- c(temp2$count[1], diff(temp2$count), temp2$nViable[1] - tail(temp2$count, 1))
    nCum <- c(temp2$count, NA)
    #print("OK")
  }
  else{
    #Cumulate, if necessary
    nSeeds <- c(temp2$count, temp2$nViable[1] - sum(temp2$count))
    nCum <- cumsum(nSeeds)
    nCum[length(nCum)] <- NA
  }

  Dish <- c(temp2$Dish, i)
  GroupT <- temp2[,5:(length(temp2[1,]))]
  GroupT <- data.frame(GroupT)
  GroupT <- rbind(GroupT, tail(GroupT, 1))
  dataset_t <- data.frame()
  dataset_t <- data.frame(GroupT, Dish, timeBef, timeAf, nViable, nSeeds, nCum)
  result <- rbind(result, dataset_t)

  }

  result
}
