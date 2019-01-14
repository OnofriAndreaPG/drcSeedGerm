#Several germination indices exist
rm(list=ls())
counts <- c(0, 25, 12, 2, 6, 0, 0)
timeAf <- c(0, 2, 4, 6, 8, 10, 12)
df <- data.frame(timeAf, counts)
rm(counts, timeAf)
df

#Germination times (only for germinated fraction)
library(germinationmetrics)
MeanGermTime(germ.counts = df$counts, intervals = df$timeAf)
sum(df$timeAf * df$counts)/sum(df$counts)
tmu <- sum(df$timeAf * df$counts)/sum(df$counts)

#Variance of germination time (only for germinated fraction)
VarGermTime(germ.counts = df$counts, intervals = df$timeAf)
sum((df$timeAf - tmu) ^ 2 * df$counts )/(sum(df$counts) - 1)

#Percentiles
t50(germ.counts = df$counts, intervals = df$timeAf, method = "coolbear")
quantileSG(df$timeAf, df$counts, probs = c(30, 50), nSeeds=250, type=1, 
           rate=F)
t50(germ.counts = df$counts, intervals = df$timeAf, method = "farooq")
quantileSG(df$timeAf, df$counts, probs = c(30, 50), nSeeds=250, type=1, 
           rate=T, se.fit = T)
#Germination extent
pMaxFin(df$timeAf, df$counts, nSeeds = 90)

#Germination rate
MeanGermRate(germ.counts = df$counts, intervals = df$timeAf)
1/tmu
CVG(germ.counts = df$counts, intervals = df$timeAf)
1/tmu * 100


GermRateRecip(germ.counts = df$counts, intervals = df$timeAf)
1/t50(germ.counts = df$counts, intervals = df$timeAf, method = "coolbear")

#Mixed
#Area under the curve
TimsonsIndex(germ.counts = df$counts, intervals = df$timeAf, total.seeds = 50)
sum(cumsum(df$counts) / 50 * 100) * 2 #Lunghezza intervallo
#sum(df$counts/50 * 100 * (15 - df$timeAf + 1))

PeakValue(germ.counts = df$counts, intervals = df$timeAf, total.seeds = 25)
max(df$nCum/25 * 100 / df$timeAf)

GermValue(germ.counts = df$counts, intervals = df$timeAf, total.seeds = 25)$Germ
max(df$nCum/25 * 100 / df$timeAf) * mean(df$counts/25 * 100)

GermSpeed(germ.counts = df$counts, intervals = df$timeAf)
sum(df$counts / df$timeAf)

#Uniformity
CUGerm(germ.counts = df$counts, intervals = df$timeAf)
varmu <- sum((df$timeAf - tmu) ^ 2 * df$counts )/sum(df$counts)
1/varmu

CVSEGermTime(germ.counts = df$counts, intervals = df$timeAf)
sd <- sqrt ( sum((df$timeAf - tmu) ^ 2 * df$counts )/(sum(df$counts) - 1) )
sd / tmu  #???

#Synchrony
GermUncertainty(germ.counts = df$counts, intervals = df$timeAf)
temp <- (df$counts/sum(df$counts) * log2(df$counts/sum(df$counts)) )[ is.na(df$counts/sum(df$counts) * log2(df$counts/sum(df$counts)) ) == F]
- sum( ( temp) )

GermSynchrony(germ.counts = df$counts, intervals = df$timeAf)
#?????

#The question is: are these indicators bioased?
#Get reasonable values for Monte Carlo simulation
library(drcSeedGerm)
library(drc)
library(lattice)
xyplot(propCum ~ timeAf|factor(Psi), data = rape)
data(rape)
head(rape)

rapeR <- subset(rape, Psi > -0.95 & Psi < -0.6)
mod <- drm(nSeeds ~ timeBef + timeAf, data = rapeR, fct=LL.3(), 
           curveid = Psi, type = "event")
summary(mod)
plot(mod, log="", xlim=c(1E-1, 15))

#Data generation: ##########################################################
#1: generate real unobservable germination times
rm(list=ls())
library(germinationmetrics)
library(drc)
library(actuar)
#quantileSG(df$timeAf, df$counts, c(10, 50, 90, 99), nSeeds=50)

GermInd <- function(nSeeds, timeLast, stepOss, mu, shape, pMax){
    germ <- rbinom(1, nSeeds, pMax)
    p <- rllogis(germ, shape = shape, scale = mu)
    GT <- c(p, rep(Inf, nSeeds - germ))
    
    #Generate the observed data
    obsT <- seq(1, timeLast, by=stepOss) #Observation schedule (by 5 o 2)
    counts <- table( cut(GT, breaks = c(0, obsT)) )
    df <- data.frame( timeAf = obsT, counts = as.numeric(counts), nCum = cumsum(counts))
    df2 <- data.frame(timeBef = c(0, df$timeAf), timeAf = c(df$timeAf, Inf), 
    counts = c(df$counts, nSeeds - sum(df$counts))) 
    
    #Calculate the statistics
    nMax <- max(df$nCum)
    #t501 <- t50(germ.counts = df$counts, intervals = df$timeAf, method = "coolbear")
    #t502 <- t50(germ.counts = df$counts, intervals = df$timeAf, method = "farooq")
    t501 <- quantileSG(df$timeAf, df$counts, 50, nSeeds=nMax, type=2)
    t502 <- quantileSG(df$timeAf, df$counts, 50, nSeeds=nMax)
    t503 <- quantileSG(df$timeAf, df$counts, 50, nSeeds=nSeeds)
    MGT <- MeanGermTime(germ.counts = df$counts, intervals = df$timeAf)
    CUG <- CUGerm(germ.counts = df$counts, intervals = df$timeAf)
    Pmax <- nMax/nSeeds
    
    df3 <- df2[df2$counts > 0 | is.finite(df2$timeAf) == F,]
    df4 <- df2[is.finite(df2$timeAf) == T,]
    df4$propCum <- cumsum(df4$counts)/nSeeds
    
    #summary(modNl, log="")
    # plot(modNl)
    #plot(mod, log="")
    
    
    # IMPORTANT: Exclude zeros and always include right censoring even if = 0
    #Try fitting
    mod1 <- try (drm( counts ~ timeBef + timeAf, data = df2, fct = LL.3(), 
      type="event", upperl = c(NA, 1, NA), lowerl = c(NA, 0, NA)))
    mod2 <- try (drm( counts ~ timeBef + timeAf, data = df3, fct = LL.3(), 
      type="event", upperl = c(NA, 1, NA), lowerl = c(NA, 0, NA)))
    mod3 <- try (drm( counts ~ timeBef + timeAf, data = df2, fct = LL.2(), 
                     type="event"))
    modNl <- try( drm( propCum ~ timeAf, data = df4, fct = LL.3(),
                       upperl = c(NA, 1, NA), lowerl = c(NA, 0, NA)) )
    modNl2 <- try( drm( propCum ~ timeAf, data = df4, fct = LL.2() ) )
    if(class(mod1) == "drc"){mod <- mod1
        t50f1 <- ED(mod, 50, display=F)[1]
        t50f2 <- ED(mod, 0.5, type = "absolute", display=F)[1]
        MGTf <- as.numeric( (coef(mod)[3] * pi / (- coef(mod)[1])) / (sin( pi/ (- coef(mod)[1]) )) )
        Pmaxf <- as.numeric( coef(mod)[2] )
        bPar <- as.numeric(coef(mod)[1])
        mu <- coef(mod)[3]; b <- - pi/bPar
        CUGf <- 1 / ( mu^2 * (2*b/sin(2*b) - (b^2)/(sin(b)^2)) )
        
        code <- 1
    }else{if(class(mod2) == "drc"){ mod <- mod2
        cat("Fitting LL3red")
        t50f1 <- ED(mod, 50, display=F)[1]
        t50f2 <- ED(mod, 0.5, type = "absolute", display=F)[1]
        MGTf <- as.numeric( (coef(mod)[3] * pi / (- coef(mod)[1])) / (sin( pi/ (- coef(mod)[1]) )) )
        Pmaxf <- as.numeric( coef(mod)[2] )
        mu <- coef(mod)[3]; b <- - pi/(coef(mod)[1])
        CUGf <- 1 / ( mu^2 * (2*b/sin(2*b) - (b^2)/(sin(b)^2)) )
        bPar <- as.numeric(coef(mod)[1])
        code <- 2
    } else {if(class(mod3) == "drc"){ mod <- mod3
        cat("Fitting LL2")
        #print(df2); stop()
        t50f1 <- ED(mod, 50, display=F)[1]
        t50f2 <- ED(mod, 0.5, type = "absolute", display=F)[1]
        MGTf <- as.numeric( (coef(mod)[2] * pi / (- coef(mod)[1])) / (sin( pi/ (- coef(mod)[1]) )) )
        Pmaxf <- 1
        mu <- coef(mod)[2]; b <- - pi/(coef(mod)[1])
        CUGf <- 1 / ( mu^2 * (2*b/sin(2*b) - (b^2)/(sin(b)^2)) )
        bPar <- as.numeric(coef(mod)[1])
        code <- 3
    } else {t50f1 = NA; t50f2 = NA; MGTf = NA; Pmaxf = NA}
    }}
    if(class(modNl) == "drc"){mod <- modNl
        t50n1 <- ED(mod, 50, display=F)[1]
        t50n2 <- ED(mod, 0.5, type = "absolute", display=F)[1]
        MGTn <- as.numeric( (coef(mod)[3] * pi / (- coef(mod)[1])) / (sin( pi/ (- coef(mod)[1]) )) )
        Pmaxn <- as.numeric( coef(mod)[2] )
        bParn <- as.numeric( coef(mod)[1] )
        coden <- 1
    } else {if(class(modNl2) == "drc"){ mod <- modNl2
        cat("Fitting LL2")
        #print(df2); stop()
        t50n1 <- ED(mod, 50, display=F)[1]
        t50n2 <- ED(mod, 0.5, type = "absolute", display=F)[1]
        MGTn <- as.numeric( (coef(mod)[2] * pi / (- coef(mod)[1])) / (sin( pi/ (- coef(mod)[1]) )) )
        Pmaxn <- 1
        bParn <- as.numeric( coef(mod)[1] )
        coden <- 2
    } else {t50n1 = NA; t50n2 = NA; MGTn = NA; Pmaxn = NA} }
    #result <- c(t501, t502, t503, MGT, CUG, Pmax, code, t50f1, t50f2, MGTf, CUGf, Pmaxf, bPar, t50n1, t50n2, MGTn, Pmaxn, bParn)
    result <-  list(t501, t502, t503, MGT, CUG, Pmax, code, t50f1, t50f2, MGTf, CUGf, Pmaxf, bPar, t50n1, t50n2, MGTn, Pmaxn, bParn)
    return(result)
}

############################################################
#First seed lot
mu <- 5; shape <- 15; pMax <- 0.9
curve( pMax * pllogis(x, shape = shape, scale = mu), xlim=c(0, 50 ), ylim=c(0, 1))
RMGT <- (mu * pi / shape) / (sin( pi/shape)) #5.045759
mode <- mu * ( (shape - 1)/( shape + 1) )^(1 / shape) #3.857565
b <- pi/shape
var <- mu^2 * (2*b/sin(2*b) - (b^2)/(sin(b)^2)) #8.4169
CUG <- 1/var
t50m <- qllogis(0.5, shape = shape, scale = mu) #4.47
t50 <- qllogis(0.5/pMax, shape = shape, scale = mu) #4.811284
RMGT; mode; var; CUG; t50m; t50
# [1] 5.344797
# [1] 4.61054
# [1] 4.465809
# [1] 0.2239236
# [1] 5
# [1] 5.228198

nSeeds <- 25; timeLast <- 15; stepOss <- 2
set.seed(1234)
result <- data.frame()
for(i in 1:1000){
  res <- GermInd(nSeeds, timeLast, stepOss, mu, shape, pMax)
  result <- rbind(result, res)
}
names(result) <- c("t501", "t502", "t503", "MGT", "CUG", "Pmax",
                   "code", "t50f1", "t50f2", "MGTf", "CUGf", "Pmaxf",
                   "bPar", "t50n1", "t50n2", "MGTn", "Pmaxn", "bParn")
#apply(result, 2, mean)
meanSim <- function(x) mean(x[is.na(x)==F & result$Pmaxf > 0])
nSim <- function(x) length(x[is.na(x)==F & result$Pmaxf > 0])
varSim <- function(x) var(x[is.na(x)==F & result$Pmaxf > 0])

apply(result, 2, meanSim)
apply(result, 2, nSim)
apply(result, 2, varSim)
(0.9*0.1)/25

table(result$code)
write.csv(result, file="A-25-15-1.csv")

nSeeds <- 25; timeLast <- 30; stepOss <- 4
set.seed(1234)
result <- data.frame()
for(i in 1:1000){
  res <- GermInd(nSeeds, timeLast, stepOss, mu, shape, pMax)
  result <- rbind(result, res)
}
names(result) <- c("t501", "t502", "t503", "MGT", "CUG", "Pmax",
                   "code", "t50f1", "t50f2", "MGTf", "CUGf", "Pmaxf",
                   "bPar", "t50n1", "t50n2", "MGTn", "Pmaxn", "bParn")
#apply(result, 2, mean)
meanSim <- function(x) mean(x[is.na(x)==F & result$Pmaxf > 0])
nSim <- function(x) length(x[is.na(x)==F & result$Pmaxf > 0])

apply(result, 2, meanSim)
apply(result, 2, nSim)
table(result$code)
write.csv(result, file="A-25-15-1.csv")



############################################################
nSeeds <- 200; timeLast <- 300; stepOss <- 4
mu <- 20.47; shape <- 2.1; pMax <- 0.51
set.seed(1234)
result <- data.frame()
for(i in 1:1000){
  res <- GermInd(nSeeds, timeLast, stepOss, mu, shape, pMax)
  result <- rbind(result, res)
}
names(result) <- c("t501", "t502", "t503", "MGT", "CUG", "Pmax",
                   "code", "t50f1", "t50f2", "MGTf", "CUGf", "Pmaxf",
                   "bPar", "t50n1", "t50n2", "MGTn", "Pmaxn", "bParn")
#apply(result, 2, mean)
#result <- read.csv("seconda.csv", header=T)

meanSim <- function(x) mean(x[is.na(x)==F & result$Pmaxf > 0 & result$Pmaxf <= 1])
meanSim <- function(x) mean(x[is.na(x)==F & result$Pmaxf > 0 & result$Pmaxf <= 1])
nSim <- function(x) length(x[is.na(x)==F & result$Pmaxf > 0 & result$Pmaxf <= 1])

apply(result, 2, meanSim)
apply(result, 2, nSim)

table(result$code)
result[is.na(result$t50n1)==T,]
write.csv(result, file="seconda.csv")

MGTll <- function(alpha, beta) { ( alpha * pi / beta ) / (sin( pi/ beta))}
str(result)
MGTll(result[1,]$t50f1, 0.001)

#Second seed lot
rm(list=ls())
mu <- 20.47; shape <- 2.1; pMax <- 0.51
curve( pMax * pllogis(x, shape = shape, scale = mu), xlim=c(0, 300 ), ylim=c(0, 1))

RMGT <- (mu * pi / shape) / (sin( pi/shape)) #30.70892
mode <- mu * ( (shape - 1)/( shape + 1) )^(1 / shape) #12.49818
b <- pi/shape
var <- mu^2 * (2*b/sin(2*b) - (b^2)/(sin(b)^2)) #7468.721
t50m <- qllogis(0.5, shape = shape, scale = mu) #20.47
t50 <- qllogis(0.5/pMax, shape = shape, scale = mu) #131.8716
RMGT; mode; var; t50m; t50

############################################################
#Third seed lot
nSeeds <- 25; timeLast <- 15; stepOss <- 1
mu <- 5; shape <- 5; pMax <- 0.95
set.seed(1234)
result <- data.frame()
for(i in 1:1000){
  res <- GermInd(nSeeds, timeLast, stepOss, mu, shape, pMax)
  result <- rbind(result, res)
}
names(result) <- c("t501", "t502", "t503", "MGT", "CUG", "Pmax",
                   "code", "t50f1", "t50f2", "MGTf", "CUGf", "Pmaxf",
                   "bPar", "t50n1", "t50n2", "MGTn", "Pmaxn", "bParn")
#apply(result, 2, mean)
meanSim <- function(x) mean(x[is.na(x)==F])
nSim <- function(x) length(x[is.na(x)==F])

apply(result, 2, meanSim)
apply(result, 2, nSim)

table(result$code)
#result[is.na(result$t50n1)==T,]
write.csv(result, file="terza.csv")

#Third seed lot
curve( pMax * pllogis(x, shape = shape, scale = mu), xlim=c(0, 100 ), ylim=c(0, 1))

RMGT <- (mu * pi / shape) / (sin( pi/shape)) #5.344797
mode <- mu * ( (shape - 1)/( shape + 1) )^(1 / shape) #4.61054
b <- pi/shape
var <- mu^2 * (2*b/sin(2*b) - (b^2)/(sin(b)^2)) #4.465809
1/4.465809
t50m <- qllogis(0.5, shape = shape, scale = mu) #5
t50 <- qllogis(0.5/pMax, shape = shape, scale = mu) #5.106478
RMGT; mode; var; t50m; t50

################################################################
#MonteCarlo package
nSeeds_g <- c(10, 25, 50, 100)
timeLast_g <- c(10, 15, 20, 40, 80)
stepOss_g <- c(0.5, 1, 2, 4)
mu_g <- c(5, 20); shape_g <- c(5, 15); pMax_g <- c(0.6, 0.9)
param_list <- list("nSeeds" = nSeeds_g,
                "timeLast" = timeLast_g,
                "stepOss" = stepOss_g,
                "mu" = mu_g, "shape" = shape_g,
                "pMax" = pMax_g)
set.seed(1234)
library(MonteCarlo)
MC_result <- MonteCarlo(func=GermInd, nrep=10, param_list = param_list)
summary(MC_result)  
mean(MC_result[[1]][[1]][1,1,1,1,1,1,1:1000])
mean(MC_result[[1]][[2]][1,1,1,1,1,1,1:1000])
mean(MC_result[[1]][[3]][1,1,1,1,1,1,1:1000])
mean(MC_result[[1]][[4]][1,1,1,1,1,1,1:1000])
mean(MC_result[[1]][[5]][1,1,1,1,1,1,1:1000])
mean(MC_result[[1]][[6]][1,1,1,1,1,1,1:1000])
mean(MC_result[[1]][[7]][1,1,1,1,1,1,1:1000])
mean(MC_result[[1]][[8]][1,1,1,1,1,1,1:1000])

#Variability of Pmax (1/10 of binomial!)
rm(list=ls())
library(drcSeedGerm)
data(barley)
rape <- barley
head(rape)
rapeR <- rape[is.finite(rape$timeAf)==F,]
head(rapeR)
rapeR$prop <- (25 - rapeR$count)/25
mod <- glm(cbind(50 - nSeeds, nSeeds) ~ Temp, data = rapeR, 
           family = binomial)
summary(mod)
library(plyr)
means <- ddply(rapeR, ~Stage, summarise, mean(prop), var(prop))
means
means[,2] * (1 - means[,2])/ means[,3]
