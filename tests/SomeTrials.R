rm(list=ls())
library(drcSeedGerm)
data(rape)
head(rape)
library(lattice)
library(aomisc)
xyplot(propCum ~ timeAf|factor(Psi), data = rape)
dataset <- subset(rape, Psi == -1.1)
mod3 <- drm(nSeeds ~ timeBef + timeAf, fct = LL.3(),
            data = dataset, type = "event")
summary(mod3)
mod3 <- drm(nSeeds ~ timeBef + timeAf, fct = LLsg.3(fixed = c(0, NA, 1)),
            data = dataset, type = "event")
summary(mod3)
plot(mod3)
dataset$timeAf
1/(1 + exp(-coef(mod3)[1]))
exp(coef(mod3)[1])




curve(expoLog.fun(x, 1, 0.5), xlim = c(0, 20000))

mod <- nls(propCum ~ expoLog.fun(timeAf, d, b), data = dataset,
           start = list(d = 1, b = 0.1))
summary(mod)
mod2 <- nls(propCum ~ NLS.negExp(timeAf, d, b), data = dataset,
           start = list(d = 1, b = 0.1))
summary(mod2)
mod3 <- drm(propCum ~ timeAf, fct = DRC.negExp(), data = dataset)
summary(mod3)
plot(mod3, xlim = c(0, 15), log = "")
plotnls(mod2, xlim = c(0, 15))

dataset <- subset(rape, Psi > -0.7)
mod <- drm(propCum ~ timeAf, fct = DRC.negExp(), curveid = factor(Psi),
           data = dataset)
summary(mod)

plot(mod, log = "")

thetaHT <- 9.68; Tb <- 3.43; Td <- 16.30; Psib50 <- -0.45; sigmaPsib <- 0.23
Kt <- 0.026

curve(HTTnormRF.fun(30, -0.5, x, thetaHT, Tb, Td, Psib50, Kt, sigmaPsib), xlim = c(0, 100),
      ylim = c(0,1))

curve(HTTnorm.fun(30, 0, x, thetaHT, Tb, Psib50, Kt, sigmaPsib), xlim = c(0, 100),
      ylim = c(0,1))
