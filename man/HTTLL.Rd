\name{HTTLL}
\alias{HTTLL.M}
\alias{HTTLL.M.fun}
\alias{HTTLL.BS}
\alias{HTTLL.BS.fun}

\title{
Hydro-thermal-time model with log-logistic distribution of base water potential (Mesgaran et al., 2017)
}

\description{
This model relates the time-course of the proportion of germinated seeds to the water potential and temperature in the substrate and it is based on a log-logistic distribution of base water potential within the seed lot. Two similar functions are available within the 'drcSeedGerm' package: the first one is 'HTTLL.M()' that assumes that the base water potential decreases with temperature for any \eqn{T > T_b}. The equation is:

\deqn{ P(t, T, \Psi) = \frac{1}{1 + \exp \left\{ -  \frac{ \log \left[ \Psi - \left( \frac{\theta_{HT}}{t ( T - T_b )} \right) + \delta \right] - \log \left[ \Psi_b + K_t (T - T_b) + \delta \right] }{\sigma_{\Psi_b}} \right\} } }


The second function is 'HTTLL.BS()' and it assumes that the base water potential decreases with temperature only for any \eqn{T > T_o}, where \eqn{T_o} is the optimal temperature level. For this case, the element \eqn{K_t(T - T_b)} is modified as \eqn{ K_t[\max(T,T_o) - T_o]} and the optimal temperature \eqn{T_o} is included as an explicit parameter.

The 'HTTLL.M.fun()' and 'HTTLL.BS.fun()' are two generic functions, which can be used for plotting or other applications, while the 'HTTLL.M()' and 'HTTLL.BS()' functions are meant to be used for model fitting with the 'drmte()' function in the 'drcte()' package.
}

\details{
The detail of this time-to-event model and the meaning of parameters are described in Mesgaran et al. (2017).
}

\usage{
HTTLL.M()
HTTLL.BS()
HTTLL.M.fun(time, Psi, Temp, thetaHT, Tb, Psib50, Kt, delta, sigmaPsib)
HTTLL.BS.fun(time, Psi, Temp, thetaHT, Tb, To, Psib50, Kt, delta, sigmaPsib)
}
\arguments{
The 'HTTLL.M()' and 'HTTLL.BS()' functions have no arguments. The 'HTTLL.M.fun()' and the 'HTTLL.BS.fun()' functions have the following arguments:
  \item{time}{time}
  \item{Psi}{water potential in the substrate}
  \item{Temp}{temperature}
  \item{thetaHT}{hydro-time parameter}
  \item{Tb}{base temperature}
  \item{To}{optimal temperature}
  \item{Psib50}{median base water potential}
  \item{Kt}{parameter measuring the effect of temperature on base water potential}
  \item{delta}{shifting parameter to allow for a log-logistic distribution of negative values for base water potential}
  \item{sigmaPsib}{scale parameter for log-logistic distribution}
}
\value{
The 'HTTLL.M.fun()' and 'HTTLL.BS.fun()' functions return the proportion of germinated seeds, for any given values of time, water potential and temperature in the substrate. The 'HTTLL.M()' and 'HTTLL.BS()' functions return a list containing the nonlinear function, the self-starter function, the parameter names and other items which are internally used by the 'drmte()' function.
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
modHTTLL.M <- drmte(nSeeds ~ timeBef + timeAf + water + temp,
                 data=hordeum,
                 fct = HTTLL.M(),
  start=c(832,-2.5, -3, 0.07, 3, 0.5))
summary(modHTTLL.M)

modHTTLL.BS <- drmte(nSeeds ~ timeBef + timeAf + water + temp,
                 data=hordeum,
                 fct = HTTLL.BS(),
  start=c(932,-2.5, 15, -3, 0.07, 3, 0.5))
summary(modHTTLL.BS)
}
\keyword{ seed germination }
\keyword{ hydro-thermal-time model }
\keyword{ time-to-event model }
