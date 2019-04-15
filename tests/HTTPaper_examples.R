library(drcSeedGerm)
library(lmtest)
library(sandwich)


#HTE model fitting - Code Snippets 1-4
data(rape)
modHTE <- drm( nSeeds ~ timeBef + timeAf + Psi,
            data=rape, fct=HTE1(), type="event")
summary(modHTE)

# Model fitted: Hydro-time model with shifted exponential
# for Pmax and linear model for GR50
#
# Parameter estimates:
#
#                         Estimate Std. Error  t-value   p-value
# G:(Intercept)          0.9577943  0.0063663  150.448 < 2.2e-16 ***
# Psib:(Intercept)      -1.0397178  0.0047014 -221.152 < 2.2e-16 ***
# sigmaPsib:(Intercept)  0.1108836  0.0087593   12.659 < 2.2e-16 ***
# thetaH:(Intercept)     0.9060853  0.0301585   30.044 < 2.2e-16 ***
# b:(Intercept)          4.0272972  0.1960877   20.538 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

coeftest(modHTE, vcov=vcovCL, cluster=rape$Dish)

# t test of coefficients:
#
#                         Estimate Std. Error   t value  Pr(>|t|)
# G:(Intercept)          0.9577943  0.0080918  118.3661 < 2.2e-16 ***
# Psib:(Intercept)      -1.0397178  0.0047067 -220.9003 < 2.2e-16 ***
# sigmaPsib:(Intercept)  0.1108836  0.0121872    9.0983 < 2.2e-16 ***
# thetaH:(Intercept)     0.9060853  0.0410450   22.0754 < 2.2e-16 ***
# b:(Intercept)          4.0272972  0.1934579   20.8174 < 2.2e-16 ***

#Code snippet 4 ############################################
ED(modHTE, Psi=-1, respLev=c(50, 30, 10))

# Estimated effective doses
#
#         Estimate Std. Error
# e:1:50 0.0438345  0.0042176
# e:1:30 0.0540987  0.0029981
# e:1:10 0.0756414  0.0029027

ED(modHTE, Psi=-1, respLev=c(50, 30, 10),
     type="absolute")

# Estimated effective doses
#
#        Estimate Std. Error
# e:1:50 0.000000         NA
# e:1:30 0.000000         NA
# e:1:10 0.051297   0.006222


#Code snippet 5 ##############################
data(hordeum)
head(hordeum)
modHTTE <- drm(nSeeds ~ timeBef + timeAf + water + temp, data=hordeum,
  fct=HTTEM(), type="event",
  start=c(0.8,-2, 0.05, 3, 0.2, 2000, 0.5))

#Sandwich standard errors
coeftest(modHTTE, vcov=vcovCL, cluster=hordeum$Dish)
# t test of coefficients:
#
#                          Estimate  Std. Error  t value Pr(>|t|)
# G:(Intercept)          9.8820e-01  1.1576e-02  85.3692   <2e-16 ***
# Psib:(Intercept)      -2.9133e+00  3.4812e-02 -83.6874   <2e-16 ***
# kt:(Intercept)         7.4228e-02  1.2666e-03  58.6015   <2e-16 ***
# Tb:(Intercept)        -7.4525e-01  3.5254e-01  -2.1139   0.0346 *
# sigmaPsib:(Intercept)  5.5284e-01  2.8976e-02  19.0790   <2e-16 ***
# ThetaHT:(Intercept)    1.3091e+03  4.0638e+01  32.2130   <2e-16 ***
# b:(Intercept)          4.1650e+00  1.1332e-01  36.7548   <2e-16 ***


jackGroupSE(modHTTE, hordeum, cluster=hordeum$Dish) #Takes long!
#                            Estimate           SE    Robust SE
# G:(Intercept)            0.98819850  0.013094129  0.011940148
# Psib:(Intercept)        -2.91328602  0.035582631  0.030384971
# kt:(Intercept)           0.07422745  0.001365243  0.001276247
# Tb:(Intercept)          -0.74525205  0.469169426  0.168407632
# sigmaPsib:(Intercept)    0.55283604  0.032776129  0.026671776
# ThetaHT:(Intercept)   1309.05984297 50.155730920 23.381767805
# b:(Intercept)            4.16504599  0.092917626  0.104963192

# Code snippet 6
data(phalaris)


#Time-to-event model fit
modTTE <- drm(nSeeds ~ timeBef + timeAf + temp, fct=TTEM(),
  data=phalaris, type="event", upperl=c(1,NA,NA,NA,NA,NA))
summary(modTTE)

# Model fitted: Thermal-time model with shifted exponential for Pmax and Mesgaran model for GR50
#
# Parameter estimates:
#
#   Estimate Std. Error t-value   p-value
#   G:(Intercept)       1.0000e+00 1.1647e-01  8.5862 < 2.2e-16 ***
#   Tc:(Intercept)      4.8801e+01 2.2307e+00 21.8769 < 2.2e-16 ***
#   sigmaTc:(Intercept) 3.3369e+01 7.3972e+00  4.5110 6.452e-06 ***
#   Tb:(Intercept)      3.3513e+00 3.2648e-01 10.2650 < 2.2e-16 ***
#   ThetaT:(Intercept)  1.6556e+03 7.7549e+01 21.3485 < 2.2e-16 ***
#   b:(Intercept)       3.1570e+00 7.3112e-02 43.1805 < 2.2e-16 ***
#
coeftest(modTTE, vcov=vcovCL, cluster=phalaris$Dish)

# t test of coefficients:
#
#   Estimate Std. Error t value  Pr(>|t|)
#   G:(Intercept)       1.0000e+00 8.0768e-02 12.3811 < 2.2e-16 ***
#   Tc:(Intercept)      4.8801e+01 3.9325e+00 12.4097 < 2.2e-16 ***
#   sigmaTc:(Intercept) 3.3369e+01 5.5571e+00  6.0047 2.146e-09 ***
#   Tb:(Intercept)      3.3513e+00 5.2455e-01  6.3890 1.927e-10 ***
#   ThetaT:(Intercept)  1.6556e+03 1.2804e+02 12.9300 < 2.2e-16 ***
#   b:(Intercept)       3.1570e+00 8.1219e-02 38.8701 < 2.2e-16 ***
#


#Code snippet 7
data(barley)
barleyS <- subset(barley, Temp<40 )
barleyS$TempF <- factor(barleyS$Temp)

mod1 <- drm(nSeeds ~ timeBef + timeAf, fct=LL.3(),
      curveid=TempF,
      type="event",
      data=barleyS)
EDs <- 1/ED(mod1, 50, display=F)

EDs
#           Estimate Std. Error
# e:1:50  0.07267946   6.991125
# e:3:50  0.09969720   9.139807
# e:7:50  0.14589512  10.865188
# e:10:50 0.20769370   6.196740
# e:15:50 0.32547337  17.992770
# e:25:50 0.53596003  17.511659
# e:30:50 0.63365374  14.877071
# e:35:50 0.23827828   3.465546

modTTEM <- drm( nSeeds ~ timeBef + timeAf + Temp, data=barley,
               fct=TTEM(), type="event")
modTTERF <- drm( nSeeds ~ timeBef + timeAf + Temp, data=barley,
               fct=TTERF(), type="event")
modTTERFc <- drm(nSeeds ~ timeBef + timeAf + Temp, data=barley,
               fct=TTERFc(), type="event")
AIC(modTTEM, modTTERF, modTTERFc)

coeftest(modTTERFc, vcov=vcovCL, cluster=barley$Dish)

# t test of coefficients:
#
#                       Estimate Std. Error  t value  Pr(>|t|)
# G:(Intercept)        0.9344628  0.0066005 141.5746 < 2.2e-16 ***
# Tc:(Intercept)      35.8003756  0.0706928 506.4215 < 2.2e-16 ***
# sigmaTc:(Intercept)  3.4329529  0.1847095  18.5857 < 2.2e-16 ***
# Td:(Intercept)      33.7636739  0.2041060 165.4223 < 2.2e-16 ***
# Tb:(Intercept)      -3.4235076  0.3288439 -10.4107 < 2.2e-16 ***
# ThetaT:(Intercept)  62.8383924  3.8016598  16.5292 < 2.2e-16 ***
# b0:(Intercept)      19.5776657  6.6793792   2.9311  0.003474 **
# s:(Intercept)        0.0063950  0.0014294   4.4740 8.788e-06 ***

coeftest(modTTEM, vcov=vcovCL, cluster=barley$Dish)

# t test of coefficients:
#
#                      Estimate Std. Error t value  Pr(>|t|)
# G:(Intercept)        12.22588    2.38958  5.1163 3.899e-07 ***
# Tc:(Intercept)       82.63574   23.46074  3.5223 0.0004519 ***
# sigmaTc:(Intercept) 998.26416  319.63167  3.1232 0.0018533 **
# Tb:(Intercept)       -1.89278    0.38753 -4.8843 1.253e-06 ***
# ThetaT:(Intercept)   43.40393    4.37336  9.9246 < 2.2e-16 ***
# b:(Intercept)         5.62728    0.69158  8.1369 1.536e-15 ***

coeftest(modTTERF, vcov=vcovCL, cluster=barley$Dish)

# t test of coefficients:
#
#                       Estimate Std. Error t value  Pr(>|t|)
# G:(Intercept)        0.9353433  0.0066389 140.888 < 2.2e-16 ***
# Tc:(Intercept)      35.8076114  0.0708678 505.274 < 2.2e-16 ***
# sigmaTc:(Intercept)  3.4558548  0.1857632  18.604 < 2.2e-16 ***
# Td:(Intercept)      33.5851390  0.2207856 152.117 < 2.2e-16 ***
# Tb:(Intercept)      -2.9579836  0.2094144 -14.125 < 2.2e-16 ***
# ThetaT:(Intercept)  58.5184898  2.1283379  27.495 < 2.2e-16 ***
# b:(Intercept)        6.8949196  0.6802456  10.136 < 2.2e-16 ***
