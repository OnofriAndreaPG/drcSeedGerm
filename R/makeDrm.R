makeDrm <- function(counts, treat, nViable, moniTimes) {
  # this function reshapes a common recording sheet for germination assays
  # into the kind of dataset required by the function drm() in the drc package

  #counts is a dataframe listing for each Petri dish (rows) the number of germinated seeds
  #at each assessment time (culumns);
  #treat is a dataframe listing for each Petri dish (rows) the levels of each treatment (culumns)
  #number of rows in dataset and treat.
  #moniTimes = vector of monitoring times. Seme length as number of columns in dataset
  #nViable is a vector with the number of viable seeds per each Petri dish

  #Creating objects to store the results
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

makeDrm.new <- function(counts, treat, nViable, moniTimes) {
  # this function reshapes a common recording sheet for germination assays
  # into the kind of dataset required by the function drm() in the drc package

  # counts is a dataframe listing for each Petri dish (rows) the number of germinated seeds
  # at each assessment time (culumns);
  # treat is a dataframe listing for each Petri dish (rows) the levels of each treatment (culumns)
  # number of rows in dataset and treat.
  # moniTimes = vector of monitoring times. Seme length as number of columns in dataset
  # nViable is a vector with the number of viable seeds per each Petri dish

  moniTimesAf <- c(moniTimes, Inf)
  moniTimesBef <- c(0, moniTimes)
  last <- nViable - apply(counts, 1, sum, na.rm = "na.omit")
  nDish <- length(counts[,1])
  Dish_tmp <- 1:nDish
  counts <- cbind(counts, last)
  colnames(counts) <- moniTimesAf
  dFrame <- cbind(Dish_tmp, treat, nViable, counts)
  firstCol <- length(treat[1,]) + 3
  lastCol <- length(dFrame[1,])
  # head(dFrame)
  #dFrame <- tibble::tibble(dFrame, .name_repair = c("unique"))
  dFrame2 <- tidyr::pivot_longer(dFrame, names_to = "timeAf",
                                  values_to = "nSeeds",
             cols = c(all_of(firstCol):all_of(lastCol)))
  # reshape::melt(dFrame, # id.vars = c("timeAf"),
  #               measure.vars = c(firstCol:lastCol),
  #               variable_name = "nSeeds")
  dFrame2 <- as.data.frame(dFrame2)

  dFrame2 <- na.omit(dFrame2)
  timeBef <- c(0, dFrame2$timeAf[-length(dFrame2$timeAf)])
  timeBef[timeBef == "Inf"] <- "0"
  dFrame2$timeBef <- timeBef
  dFrame2 <- plyr::ddply(dFrame2, .(Dish_tmp), .fun=transform,
                         nCum = cumsum(nSeeds))
  finData <- cbind(dFrame2[,2:(firstCol - 2)],
                   nViable = dFrame2$nViable,
                   timeBef = as.numeric(dFrame2$timeBef),
                   timeAf = as.numeric(dFrame2$timeAf),
                   nSeeds = as.numeric(dFrame2$nSeeds),
                   nCum = as.numeric(dFrame2$nCum))
  row.names(finData) <- 1:length(finData[,1])
  return(finData)
}

makeDrm2 <- function(counts, treat, nViable, moniTimes, Dish,  cumulative=T){
  #This function transform a dataset as cumulative germinations in a dataset for
  #use with DRM type = "event. Dish is a factor that identidfies seeds
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
