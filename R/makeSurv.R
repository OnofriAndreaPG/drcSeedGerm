
makeSurv2 <- function(counts, data) {
#This function organises a dataset to be submitted to survival analysis
#i.e. one row per each seed.

# dataset is a vector of counts

  anName <- deparse(substitute(counts))
  dfr <- data["anName" > 0,]
  frequency <- subset(dfr, select = anName)[,1]
  ## print(head(data))
  # data <- data[frequency > 0,]
  nr <- nrow(dfr)
  # print(frequency)
  # print(rep(seq_len(nr), frequency))
  df_surv <- dfr[rep(seq_len(nr), frequency), ]
  row.names(df_surv) <- 1:(length(df_surv[,1]))
  df_surv
}


makeSurv.old <- function(dataset, treat, seeds, moniTimes) {
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


