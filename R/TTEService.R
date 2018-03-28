#####################################################################
#This module contains several service functions for time-to-event
#analysis with drm
#Date of last revision: 18/01/2018
#Version 1.0
#####################################################################

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
  dati <- counts
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


makeSurv <- function(dataset, treat, seeds, moniTimes) {
#This function organises a dataset to be submitted to survival analysis
#i.e. one row per each seed.

#dataset is a dataframe listing for each Petri dish (rows) the number of germinated seeds
#at each assessment time (culumns);
#treat is a dataframe listing for each Petri dish (rows) the levels of each treatment (culumns)
#seeds = vector listing the number of viable seeds for each Petri dish (same length as
#number of rows in dataset and treat.
#moniTimes = vector of monitoring times. Seme length as number of columns in dataset

dati <- dataset
tratt <- treat
numSemi <- seeds
tempi <- moniTimes

#Sostituisce i NAs con 0 e prepara il file
dati[is.na(dati)] <- 0

#Dimensionamento
numPetri <- length(dati[,1]);numSeeds <- sum(numSemi);numTimes <- length(tempi);numTesi<-length(tratt[1,])
timeBef <- numeric(numSeeds);timeAf <- numeric(numSeeds);cens <- numeric(numSeeds)
tempi2 <- c(0, tempi[1:length(tempi)-1])
group <- data.frame()
final.date <- max(tempi)
seed_ind <- 1
j <- 1
for (j in 1:numPetri){ #Per ogni riga - Petri dish
germ <- 0
    for(i in 1:numTimes) { #Per ogni campionamento
	    	if(dati[j,i]!=0)
	    	for(z in 1:dati[j,i]){ #Per ogni seme germinato, scrive il tempo
	        	timeAf[seed_ind] <- tempi[i]
	        	timeBef[seed_ind] <- tempi2[i]
	        	cens[seed_ind] <- 1
	        	seed_ind <- seed_ind + 1
	        	germ <- germ + 1
	        }
	 }
	#Per i semi non germinati, scrive di conseguenza
	nonGerm <- numSemi[j] - germ
	if(nonGerm>0) {for(z in 1:nonGerm){timeAf[seed_ind] <- NA; timeBef[seed_ind] <- tempi[i]
			seed_ind <-  seed_ind + 1}}
  #print(paste(j,germ,nonGerm,numSemi[j], length(timeAf)))
	}
group <- tratt[rep(row.names(tratt), numSemi), 1:length(tratt[1,])]
datiFin <- data.frame(group, timeBef=timeBef, timeAf=timeAf, status=cens)#, event=event, time=time, time2=time2)
return(datiFin)
}

#ServiceFunctions
#Log-Logistic Function for bioassay work nlsLL.3
NLSLL.3mean <- function(predictor, b, d, ED50) {
                      x <- predictor
                      d/(1+exp(b*(log(x+0.000001)-log(ED50))))
}

NLSLL.3Init <- function(mCall, LHS, data) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          d <- max(y) * 1.05
          ## Linear regression on pseudo y values
          pseudoY <- log((d-y)/(y+0.00001))
          coefs <- coef( lm(pseudoY ~ log(x+0.000001)))
          k <- -coefs[1]; b <- coefs[2]
          ED50 <- exp(k/b)
          value <- c(b,d,ED50)
          names(value) <- mCall[c("b", "d", "ED50")]
          value
}

NLSLL.3 <- selfStart(NLSLL.3mean, NLSLL.3Init, parameters=c("b", "d", "ED50"))
