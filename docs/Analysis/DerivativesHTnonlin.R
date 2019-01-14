delta <- 1.0846335
mu <- -0.9363510
sigma <- 0.7817727
thetaH <- 0.7591181
g <- 0.2072411
Psi <- -1
time <- 20

#Expression ###############################################
HTW2.fun(time, Psi, thetaH, delta, mu, sigma)
expr <- expression( 1 - exp(-exp((log((Psi - thetaH/time + delta)/(mu + delta))/sigma)) ) )
eval(expr)

#Derivata parziale thetaH ####################################
D(expr, "thetaH" )
eval(D(expr, "thetaH" ))

#Approximation (finite differences)
d1.1 <- HTW2.fun(time, Psi, thetaH, delta, mu, sigma)
d1.2 <- HTW2.fun(time, Psi, thetaH + 10e-06, delta, mu, sigma)
d1 <- (d1.2 - d1.1)/10e-6
d1

#Derivata parziale Psib ####################################
D(expr, "delta" )
eval(D(expr, "delta" ))

#Approximation (finite differences)
d2.1 <- HTW2.fun(time, Psi, thetaH, delta, mu, sigma)
d2.2 <- HTW2.fun(time, Psi, thetaH, delta + 10e-06, mu, sigma)
d2 <- (d2.2 - d2.1)/10e-6
d2

