thetaH <- 0.77466142
Psib50 <- -0.96329603
sigmaPsib <- 0.09244382
g <- 0.2072411
Psi <- -1
ret <- qnorm(g)

#Expression (type = absolute) ###############################################
expr <- expression( - (sigmaPsib * ret - Psi + Psib50 )/ thetaH )
eval(expr)
HTnorm.gra <- function(thetaH, Psib50, sigmaPsib, Psi, g) {
 GR <- - (sigmaPsib * qnorm(g) - Psi + Psib50 )/ thetaH
 GR <- ifelse(GR > 0, GR, 0)
 } 

#Derivata parziale thetaH ####################################
D(expr, "thetaH" )
eval(D(expr, "thetaH" ))
d1 <- (sigmaPsib * ret - Psi + Psib50)/thetaH^2

#Approximation (finite differences)
g <- 0.5
d1 <- HTnorm.gra(thetaH, Psib50, sigmaPsib, Psi, g)
d2 <- HTnorm.gra(thetaH + 10e-6, Psib50, sigmaPsib, Psi, g)
(d2 - d1)/10e-6

#Derivata parziale Psib ####################################
D(expr, "Psib" )
eval(D(expr, "Psib" ))
-(exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
    g)) + log(thetaH/((Psi - Psib)^2))) * ((1/b) * ((1/g) * (G *
    (exp(-(Psi - Psib) * (1/sigmaPsib)) * (1/sigmaPsib)))/((1/g) *
    (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) - g))) + thetaH *
    (2 * (Psi - Psib))/((Psi - Psib)^2)^2/(thetaH/((Psi - Psib)^2)))/exp(-(1/b) *
    log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
        g)) + log(thetaH/((Psi - Psib)^2)))^2)

#Approximation (finite differences)
d1 <- HTE2.gra(G, Psib, sigmaPsib, thetaH, b)
d2 <- HTE2.gra(G, Psib + 10e-6, sigmaPsib, thetaH, b)
(d2 - d1)/10e-6

#Derivata parziale sigmaPsib ####################################
D(expr, "sigmaPsib" )
eval(D(expr, "sigmaPsib" ))
-(exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
    g)) + log(thetaH/((Psi - Psib)^2))) * ((1/b) * ((1/g) * (G *
    (exp(-(Psi - Psib) * (1/sigmaPsib)) * ((Psi - Psib) * (1/sigmaPsib^2))))/((1/g) *
    (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) - g))))/exp(-(1/b) *
    log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
        g)) + log(thetaH/((Psi - Psib)^2)))^2)
#Approximation (finite differences)
d1 <- HTE2.gra(G, Psib, sigmaPsib, thetaH, b)
d2 <- HTE2.gra(G, Psib, sigmaPsib + 10e-6, thetaH, b)
(d2 - d1)/10e-6

#Derivata parziale sigmaPsib ####################################
parD <- "thetaH"
D(expr,  parD)
eval(D(expr, parD ))
-(exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
    g)) + log(thetaH/((Psi - Psib)^2))) * (1/((Psi - Psib)^2)/(thetaH/((Psi -
    Psib)^2)))/exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi -
    Psib) * (1/sigmaPsib))) - g)) + log(thetaH/((Psi - Psib)^2)))^2)
#Approximation (finite differences)
d1 <- HTE2.gra(G, Psib, sigmaPsib, thetaH, b)
d2 <- HTE2.gra(G, Psib, sigmaPsib, thetaH + 10e-6, b)
(d2 - d1)/10e-6

#Derivata parziale b ####################################
parD <- "b"
D(expr,  parD)
eval(D(expr, parD ))
-(exp(-(1/b) * log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
    g)) + log(thetaH/((Psi - Psib)^2))) * (1/b^2 * log((1/g) *
    (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) - g)))/exp(-(1/b) *
    log((1/g) * (G * (1 - exp(-(Psi - Psib) * (1/sigmaPsib))) -
        g)) + log(thetaH/((Psi - Psib)^2)))^2)
#Approximation (finite differences)
d1 <- HTE2.gra(G, Psib, sigmaPsib, thetaH, b)
d2 <- HTE2.gra(G, Psib, sigmaPsib, thetaH, b + 10e-6)
(d2 - d1)/10e-6

#Expression (type = relative) ###############################################
expr <- expression( 1/exp( - (1/b)*log( ((1 - g)/g) ) + log( thetaH/((Psi - Psib)^2))) )
eval(expr)
HTE2.grr <- function(G, Psib, sigmaPsib, thetaH, b)  1/exp( - (1/b)*log( ((1 - g)/g) ) + log( thetaH/((Psi - Psib)^2)))

#Derivata parziale G ####################################
D(expr, "G" )
eval(D(expr, "G" ))
0

#Approximation (finite differences)
d1 <- HTE2.grr(G, Psib, sigmaPsib, thetaH, b)
d2 <- HTE2.grr(G + 10e-6, Psib, sigmaPsib, thetaH, b)
(d2 - d1)/10e-6

#Derivata parziale Psib ####################################
D(expr, "Psib" )
eval(D(expr, "Psib" ))
-(exp(-(1/b) * log(((1 - g)/g)) + log(thetaH/((Psi - Psib)^2))) *
    (thetaH * (2 * (Psi - Psib))/((Psi - Psib)^2)^2/(thetaH/((Psi -
        Psib)^2)))/exp(-(1/b) * log(((1 - g)/g)) + log(thetaH/((Psi -
    Psib)^2)))^2)

#Approximation (finite differences)
d1 <- HTE2.grr(G, Psib, sigmaPsib, thetaH, b)
d2 <- HTE2.grr(G, Psib + 10e-6, sigmaPsib, thetaH, b)
(d2 - d1)/10e-6

#Derivata parziale sigmaPsib ####################################
D(expr, "sigmaPsib" )
eval(D(expr, "sigmaPsib" ))
0

#Approximation (finite differences)
d1 <- HTE2.grr(G, Psib, sigmaPsib, thetaH, b)
d2 <- HTE2.grr(G, Psib, sigmaPsib + 10e-6, thetaH, b)
(d2 - d1)/10e-6

#Derivata parziale sigmaPsib ####################################
parD <- "thetaH"
D(expr,  parD)
eval(D(expr, parD ))
-(exp(-(1/b) * log(((1 - g)/g)) + log(thetaH/((Psi - Psib)^2))) *
    (1/((Psi - Psib)^2)/(thetaH/((Psi - Psib)^2)))/exp(-(1/b) *
    log(((1 - g)/g)) + log(thetaH/((Psi - Psib)^2)))^2)
#Approximation (finite differences)
d1 <- HTE2.grr(G, Psib, sigmaPsib, thetaH, b)
d2 <- HTE2.grr(G, Psib, sigmaPsib, thetaH + 10e-6, b)
(d2 - d1)/10e-6

#Derivata parziale b ####################################
parD <- "b"
D(expr,  parD)
eval(D(expr, parD ))
-(exp(-(1/b) * log(((1 - g)/g)) + log(thetaH/((Psi - Psib)^2))) *
    (1/b^2 * log(((1 - g)/g)))/exp(-(1/b) * log(((1 - g)/g)) +
    log(thetaH/((Psi - Psib)^2)))^2)
#Approximation (finite differences)
d1 <- HTE2.grr(G, Psib, sigmaPsib, thetaH, b)
d2 <- HTE2.grr(G, Psib, sigmaPsib, thetaH, b + 10e-6)
(d2 - d1)/10e-6


