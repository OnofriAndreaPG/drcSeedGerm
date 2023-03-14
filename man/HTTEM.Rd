\name{HTTEM}
\alias{HTTEM}
\alias{HTTEM.fun}

\title{
Hydro-thermal-time model with log-logistic distribution of germination time (Onofri et al., 2018)
}

\description{
This model relates the time-course of the proportion of germinated seeds to the water potential and temperature in the substrate. It is based on a truncated log-logistic distribution of germination time:


\deqn{P(t) = \frac{d}{1 + exp\left[ b (\log(x) - \log(e)\right]}}

where two of the three usual parameters ('d' and 'e') are expressed as functions of water potential (\eqn{\Psi}) and temperature (\eqn{T}). In the function 'HTTEM()', we implemented the following submodels: (1) for the parameter 'd', we implemented a shifted exponential function:

\deqn{d = G \, \left[ 1 - \exp \left( \frac{ \Psi - \Psi_b - k(T - T_b )}{\sigma_{\Psi_b}} \right) \right]}

while, (2) for the parameter 'e' we considered that its inverse corresponds to the median Germination Rate within the population (i.e. \eqn{1/e = GR_{50}}) and modelled this latter parameter as:

\deqn{ GR_{50} =  \frac{T - T_b }{\theta_{HT}} \left[\Psi - \Psi_b - k(T - T_b )\right] }

The 'HTTEM.fun()' is a generic function, which can be used for plotting or other applications, while the 'HTTEM()' function is meant to be used for model fitting with the 'drmte()' function in the 'drcte()' package.
}
\details{
The detail of this time-to-event model and the meaning of parameters are described in Onofri et al. (2018). See Table 2, where 'HTTEM()' is abbreviated as HTTE.
}
\usage{
HTTEM()
HTTEM.fun(time, Psi, Temp, G, Psib, kt, Tb, sigmaPsib, ThetaHT, b)
}
\arguments{
The 'HTTEM' function has no arguments. The 'HTTEM.fun()' has the following arguments:
  \item{time}{time}
  \item{Psi}{water potential in the substrate}
  \item{Temp}{temperature}
  \item{G}{maximum germination capability of the seed lot}
  \item{Psib}{base water potential}
  \item{kt}{parameter measuring the effect of temperature on base water potential}
  \item{Tb}{base temperature}
  \item{sigmaPsib}{standard deviation for the base water potential within the seed lot}
  \item{ThetaHT}{hydro-thermal-time parameter}
  \item{b}{shape parameter for the cumulative distribution function of germination time}
}

\value{
The 'HTTEM.fun()' function returns the proportion of germinated seeds, for any given values of time and water potential in the substrate. The 'HTTEM()' function returns a list containing the nonlinear function, the parameter names and other items which are internally used by the 'drmte()' function. At the moment, there is no self-starting function and starting parameters for fitting must be provided within the 'drcte' function.
}

\references{
Onofri, A., Benincasa, P., Mesgaran, M.B., Ritz, C., 2018. Hydrothermal-time-to-event models for seed germination. European Journal of Agronomy 101, 129–139.
\url{https://www.statforbiology.com/2023/stat_drcte_10-examplehtte/}
}
\author{
Andrea Onofri
}

\examples{
data(hordeum)
modHTTEM <- drmte(nSeeds ~ timeBef + timeAf + water + temp,
                 data=hordeum,
                 fct = HTTEM(),
  start=c(0.8,-2, 0.05, 3, 0.2, 2000, 0.5))
summary(modHTTEM)
}
\keyword{ seed germination }
\keyword{ hydro-thermal-time model }
\keyword{ time-to-event model }
