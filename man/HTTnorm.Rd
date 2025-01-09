\name{HTTnorm.M}
\alias{HTTnorm.M}
\alias{HTTnorm.M.fun}
\alias{HTTnorm.BS}
\alias{HTTnorm.BS.fun}

\title{
Hydro-thermal-time model with normal distribution of base water potential (Bradford, 2002)
}

\description{
This model relates the time-course of the proportion of germinated seeds to the water potential and temperature in the substrate and it is based on a normal distribution of base water potential within the seed lot. Two similar functions are available within the 'drcSeedGerm' package: the first one is 'HTTnorm.M()' that assumes that the base water potential decreases with temperature for any \eqn{T > T_b}. The equation is:

\deqn{ P(t, T, \Psi) = \Phi \left\{ \frac{\Psi - \left[ \frac{\theta_{HT}}{t (T - T_b)} \right] - \left[\Psi_b + K_t (T - T_b) \right] }{\sigma_{\Psi_b}} \right\} }

where \eqn{\Phi} is a gaussian cumulative distribution function for base water potential.

The second function is 'HTTNorm.BS()' and it assumes that the base water potential decreases with temperature only for any \eqn{T > T_o}, where \eqn{T_o} is the optimal temperature level. For this case, the element \eqn{K_t(T - T_b)} is modified as \eqn{ K_t[\max(T,T_o) - T_o]} and the optimal temperature \eqn{T_o} is included as an explicit parameter.

The 'HTTNorm.M.fun()' and 'HTTNorm.BS.fun()' are two generic functions, which can be used for plotting or other applications, while the 'HTTNorm.M()' and 'HTTNorm.BS()' functions are meant to be used for model fitting with the 'drmte()' function in the 'drcte()' package.
}

\details{
The detail of this time-to-event model and the meaning of parameters are described in Bradford (2002).
}

\usage{
HTTnorm.M()
HTTnorm.BS()
HTTnorm.M.fun(time, Psi, Temp, thetaHT, Tb, Psib50, Kt, sigmaPsib)
HTTnorm.BS.fun(time, Psi, Temp, thetaHT, Tb, To, Psib50, Kt, sigmaPsib)
}
\arguments{
The 'HTTnorm.M()' and 'HTTnorm.BS()' functions have no arguments. The 'HTTnorm.M.fun()' and the 'HTTnorm.BS.fun()' functions have the following arguments:
  \item{time}{time}
  \item{Psi}{water potential in the substrate}
  \item{Temp}{temperature}
  \item{thetaHT}{hydro-time parameter}
  \item{Tb}{base temperature}
  \item{To}{optimal temperature}
  \item{Psib50}{median base water potential}
  \item{Kt}{parameter measuring the effect of temperature on base water potential}
  \item{sigmaPsib}{standard deviation for the base water potential within the seed lot}
}
\value{
The 'HTTnorm.M.fun()' and 'HTTnorm.BS.fun()' functions return the proportion of germinated seeds, for any given values of time, water potential and temperature in the substrate. The 'HTTnorm.M()' and 'HTTnorm.BS()' functions return a list containing the nonlinear function, the self-starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}
\references{
Bradford, K.J., 2002. Applications of hydrothermal time to quantifying and modeling seed germination and dormancy. Weed Science 50, 248–260.
}
\author{
Andrea Onofri
}
\note{
This function is for use with the function 'drmte()'
}

\examples{
data(hordeum)
\donttest{
modHTTnorm.M <- drmte(nSeeds ~ timeBef + timeAf + water + temp,
                 data=hordeum,
                 fct = HTTnorm.M(),
  start=c(932,-2.5, -3, 0.07, 0.5))
summary(modHTTnorm.M)
}
}
\keyword{ seed germination }
\keyword{ hydro-thermal-time model }
\keyword{ time-to-event model }
