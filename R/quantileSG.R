# Function to get nonparametric quantiles for germination rate and pMax
# with bootstrap SEs
# Last edited: 12/01/19
##########################################################################
# quantileSG.new <- function(timeBef, timeAf, counts, probs, # nSeeds,
#                        se.fit = F, rate = F, type = 1){
#   # Quantiles for grouped data (Farooq, 2005)
#   # time <- c(2, 4, 6, 8, 10, 12); counts <- c(22, 12, 7, 6, 0, 0)
#   # time <- c(2, 4, 6, 8, 10, 12); counts <- c(0,0,0,0,25,0)
#   # probs <- c(50); nSeeds <- sum(counts); type = 1; i <- length(probs)
#   # germinations are uniformly distributed within the interval!
#   #
#   ret <- c(); g <- probs
#   nSeeds <- sum(counts)
#   for(i in 1:length(g)){
#     # For each percentile
#     dec <- (nSeeds + type - 1) * g[i]/100
#     if( dec > sum(counts)) { ret[i] <- Inf; next() }
#     pos <- head( which(cumsum(counts) >= dec), 1)
#     t1 <- timeBef[pos]
#     t2 <- timeAf[pos]
#     N1 <- ifelse(pos == 1, 0, cumsum(counts)[pos - 1])
#     N2 <- cumsum(counts)[pos]
#     ret[i] <- t1 +  (dec - N1) * abs( (t2 - t1)/(N2 - N1) )
#   }
#   names(ret) <- paste(probs, "%", sep="")
#   if(rate == F) {  return(ret)} else {
#     if(se.fit == F){return(1/ret)} else {
#       if(se.fit == T){
#         estimate <- 1/ret
#           se <- bootSGr(timeBef, timeAf, counts, probs, type)
#           return(list( estimate = estimate, se = se)) }
#       } } }

bootSGr <- function(timeBef, timeAf, counts, probs, type = 1){
temp <- data.frame()
for(i in 1:1000){
  start <- rep(timeBef, counts)
  end <- rep(timeAf, counts)
  bootSam <- sample(1:length(start), replace = T)
  # bySeed <- c( rep(time, counts), rep(10000, nSeeds - sum(counts)) )
  start <- start[bootSam]
  end <- end[bootSam]
  # resVec <- sample(bySeed, replace = T)
  resCount <- as.numeric( table(cut(end, breaks=c(0, sort(unique(timeAf)))) ) )
  values <- quantileSG(timeBef, timeAf, resCount, probs, type)
  #pMax <- sum(resCount)/nSeeds
  report <- c(1/values)
  temp <- rbind(temp, report)
}
names(temp) <- c(paste(probs, "%", sep = "") )
#apply(temp, 2, mean)
#apply(temp, 2, median)
bootSE <- apply(temp, 2, sd)
bootSE
}

quantileSG <- function(time, counts, probs, nSeeds,
                       se.fit = F, rate = F, type = 1){
  # Quantiles for grouped data (Farooq, 2005)
  # time <- c(2, 4, 6, 8, 10, 12); counts <- c(25, 12, 2, 6, 0, 0)
  # time <- c(2, 4, 6, 8, 10, 12); counts <- c(0,0,0,0,25,0)
  # probs <- c(10, 30, 50); nSeeds <- sum(counts); type = 1; i <- length(probs)
  # germinations are uniformly distributed within the interval!
  #
  ret <- c(); g <- probs
  for(i in 1:length(g)){
    # For each percentile
    dec <- (nSeeds + type - 1) * g[i]/100
    if( dec > sum(counts)) { ret[i] <- Inf; next() }
    pos <- head( which(cumsum(counts) >= dec), 1)
    t1 <- ifelse(pos == 1, 0, time[pos - 1])
    t2 <- time[pos]
    N1 <- ifelse(pos == 1, 0, cumsum(counts)[pos - 1])
    N2 <- cumsum(counts)[pos]
    ret[i] <- t1 +  (dec - N1) * abs( (t2 - t1)/(N2 - N1) )
  }
  names(ret) <- paste(probs, "%", sep="")
  if(rate == F) {  return(ret)} else {
    if(se.fit == F){return(1/ret)} else {
      if(se.fit == T){
        estimate <- 1/ret
          se <- bootSGr.old(time, counts, probs, nSeeds, type)
          return(list( estimate = estimate, se = se)) }
      } } }

bootSGr.old <- function(time, counts, probs, nSeeds, type=1){
temp <- data.frame()
for(i in 1:1000){
  bySeed <- c( rep(time, counts), rep(10000, nSeeds - sum(counts)) )
  # set.seed(1234)
  resVec <- sample(bySeed, replace = T)
  resCount <- as.numeric( table(cut(resVec, breaks=c(0, sort(unique(time)))) ) )
  values <- quantileSG.old(time, resCount, probs, nSeeds, type)
  #pMax <- sum(resCount)/nSeeds
  report <- c(1/values)
  temp <- rbind(temp, report)
}
names(temp) <- c(paste(probs, "%", sep = "") )
#apply(temp, 2, mean)
#apply(temp, 2, median)
bootSE <- apply(temp, 2, sd)
bootSE
}

bootSGpMax <- function(time, counts, nSeeds){
  temp <- c()
  for(i in 1:1000){
    bySeed <- c( rep(time, counts), rep(10000, nSeeds - sum(counts)) )
    resVec <- sample(bySeed, replace = T)
    resCount <- as.numeric( table(cut(resVec, breaks=c(0, sort(unique(time)) )) ) )
    pMax <- sum(resCount)/nSeeds
    temp <- c(temp, pMax)
  }
  bootSE <- sd(temp)
  bootSE
}

pMaxFin <- function(time, counts, nSeeds, se.fit = F){
  pMax <- sum(counts)/nSeeds
  if(se.fit == T) {
    pMaxSE <- bootSGpMax(time, counts, nSeeds)
    returnList <- list(pMaxFin = pMax, se = pMaxSE)
    }else{ returnList <- list (pMaxFin = pMax) }
  return(returnList)
  }

