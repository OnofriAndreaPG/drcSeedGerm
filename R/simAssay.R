simAssay <- function(n, parms, numSeeds, schedule,
                     fct = "loglogistic"){
  b <- parms[1]
  d <- parms[2]
  e <- parms[3]
  nvals <- length(schedule) + 1
  dfres <- data.frame(Id = rep(NA, nvals*n),
                      timeBef = rep(NA, nvals*n),
                      timeAf = rep(NA, nvals*n),
                      counts = rep(NA, nvals*n))
  for(i in 1:n){
    # i <- 1
    #1: Generate the ungerminable seeds
    viable <- rbinom(1, numSeeds, d)
    nonViable <- numSeeds - viable

    #2: generate the germination time for each seed in the lot
    if(fct == "loglogistic"){
      lit <- survival::rsurvreg(viable, log(e), b, distribution='loglogistic')
      lit <- c(lit, rep(100000, times = numSeeds - sum(viable)))
      tg <- lit
    } else if(fct == "lognormal"){
      lit <- survival::rsurvreg(viable, log(e), b, distribution='lognormal')
      lit <- c(lit, rep(100000, times = numSeeds - sum(viable)))
      tg <- lit
    } else if(fct == "weibull"){
      lit <- survival::rsurvreg(viable, log(e), b, distribution='weibull')
      lit <- c(lit, rep(100000, times = numSeeds - sum(viable)))
      tg <- lit
    } else if(fct == "exponential"){
      lit <- survival::rsurvreg(viable, log(e), b, distribution='weibull')
      lit <- c(lit, rep(100000, times = numSeeds - sum(viable)))
      tg <- lit
    }

    # 3: generate germination intervals, according to schedule and
    #    generate the counts
    intervals <- cut(tg, breaks = c(0, schedule, Inf))
    counts <- table(intervals)
    nams <- names(counts)
    margins <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", nams) ),
        upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", nams) ))
    timeBef <- margins[ ,1]
    timeAf <- margins[ ,2]
    # timeBef <- timeBef[counts > 0]
    # timeAf <- timeAf[counts > 0]
    # counts <- counts[counts > 0]
    firstRow <- 1 + nvals * (i - 1)
    lastRow <- firstRow + nvals - 1
    dfres[(firstRow:lastRow), 1] <- i
    dfres[(firstRow:lastRow), 2] <- timeBef
    dfres[(firstRow:lastRow), 3] <- timeAf
    dfres[(firstRow:lastRow), 4] <- as.numeric(counts)
  }
  dfres <- dfres[counts > 0, ]
}
