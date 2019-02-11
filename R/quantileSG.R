#Function to get nonparametric quantiles for germination rate and pMax
#with bootstrap SEs
#Last edited: 12/01/19
#####################################################################
quantileSG <- function(time, counts, probs, nSeeds,
                       se.fit = F, rate = F, type=1){
  #time <- c(2, 4, 6, 8, 10, 12); counts <- c(25, 12, 2, 6, 0, 0)
  #probs <- 10; nSeeds <- 250; type = 1; i <- length(g)
  ret <- c(); g <- probs
  for(i in 1:length(g)){
    dec <- (nSeeds + type - 1) * g[i] /100
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
          se <- bootSGr(time, counts, probs, nSeeds, type)
          return(list( estimate = estimate, se = se)) }
      } } }

bootSGr <- function(time, counts, probs, nSeeds, type=1){
temp <- data.frame()
for(i in 1:1000){
  bySeed <- c( rep(time, counts), rep(10000, nSeeds - sum(counts)) )
  resVec <- sample(bySeed, replace = T)
  resCount <- as.numeric( table(cut(resVec, breaks=c(0, sort(unique(time)))) ) )
  values <- quantileSG(time, resCount, probs, nSeeds, type)
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

pMaxFin <- function(time, counts, nSeeds){
  pMax <- sum(counts)/nSeeds
  se <- bootSGpMax(time, counts, nSeeds)
  return(list (pMaxFin = pMax, se = se))
  }

