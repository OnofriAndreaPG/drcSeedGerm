\name{HTTEM}
\alias{HTTEM}
\alias{HTTEM.fun}

\title{
Hydro-thermal-time-to-event model for seed germination.
}
\description{
This model relates the time-course of the proportion of germinated seeds to the water potential and temperature in the substrate. It is based on four submodels: (1) a truncated log-logistic distribution of germination time, where the three usual parameters ('d', 'e' and 'b') are expressed as functions of water potential (water) and temperature (temp) in the substrate; (2) for the parameter 'd', the dependency on 'Psi'water' and 'temp' is described by using a shifted exponential distribution; (3) for the parameter 'e', it is assumed that its inverse (i.e. 1/e = GR50) is related to 'water' and 'temp' by using the function proposed in Mesgaran et al. (2017); (4) for the 'b' parameter, it is assumed that it is unaffected by the environmental covariates 'water' and 'temp'. The 'HTTEM.fun()' is a generic function, which can be used for plotting or other applications, while the HTTEM() function is meant to be used to fit this model by the 'drmte()' function in the 'drcte()' package. Self-starters are not provided for this model.
}
\usage{
HTTEM()
HTTEM.fun(time, Psi, Temp, G, Psib, kt, Tb, sigmaPsib, ThetaHT, b)
}
\details{
The detail of this time-to-event model is described in Onofri et al. (2018). See Table 2, where this model is abbreviated as HTTE.
}
\arguments{
  The HTTEM() function has no arguments. The HTTEM.fun() function has the following arguments:
  \item{time}{time}
  \item{Psi}{water potential in the substrate}
  \item{Temp}{temperature in the substrate}
  \item{G}{maximum germination capability of the seed lot}
  \item{Psib}{base water potential}
  \item{Kt}{regression parameter, expressing the effect of temperature on base water potential}
  \item{Tb}{base germination temperature}
  \item{sigmaPsib}{standard deviation for the base water potential within the seed lot}
  \item{thetaHT}{hydro-thermal-time parameter}
  \item{b}{shape parameter for the cumulative distribution function of germination time}
}
\value{
The 'HTTEM.fun()' functions returns the proportion of germinated seeds, for given values of time, water potential and temperature in the substrate, depending on model parameters. The HTTEM() function returns a list containing the nonlinear function, the self starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}
\references{
Onofri, A., Benincasa, P., Mesgaran, M.B., Ritz, C., 2018. Hydrothermal-time-to-event models for seed germination. European Journal of Agronomy 101, 129–139.
}
\author{
Andrea Onofri
}

\examples{
library(drcte)
data(hordeum)
modHTTE <- drmte(nSeeds ~ timeBef + timeAf + water + temp,
                 data=hordeum,
                 fct = HTTEM(),
  start=c(0.8,-2, 0.05, 3, 0.2, 2000, 0.5))
summary(modHTTE, units = Dish) # Cluster robust SEs
}
\keyword{ seed germination }
\keyword{ hydro-thermal-time model }
\keyword{ time-to-event model }
