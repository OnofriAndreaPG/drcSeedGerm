\name{TTEM}
\alias{TTEM}
\alias{TTEM.fun}

\title{
Thermal-time model for seed germination, relating the time-course of the proportion of germinated seeds to the environmental temperature
}
\description{
This model is based on a truncated log-logistic distribution of germination time, where the three parameters ('Pmax', 'GR50' and 'b') are expressed as functions of water potential (Psi). The 'TTEM.fun()' is a generic function, which can be used for plotting or other applications, while the TTEM() function is meant to be used to fit this model by the 'drmte()' function in the 'drcte()' package.
}
\usage{
TTEM()
TTEM.fun(time, Temp, G, Tc, sigmaTc, Tb, ThetaT, b)
}
\details{
The detail of these functions are described in Onofri et al. (2018).
}
\arguments{
  The TTEM() function has no arguments. The TTEM.fun() function has the following arguments:
  \item{time}{time}
  \item{Temp}{temperature}
  \item{G}{maximum germination capability of the seed lot}
  \item{Tc}{ceiling temperature}
  \item{sigmaTc}{standard deviation for the ceiling temperature within the seed lot}
  \item{Tb}{base temperature}
  \item{ThetaT}{Thermal-time parameter}
  \item{b}{shape parameter for the cumulative distribution function of germination time}
}
\value{
The 'TTEM.fun()' functions returns the proportion of germinated seeds, for given values of time and temperature. The TTEM() function returns a list containing the nonlinear function, the self starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}
\references{
Onofri, A., Benincasa, P., Mesgaran, M.B., Ritz, C., 2018. Hydrothermal-time-to-event models for seed germination. European Journal of Agronomy 101, 129–139.
\url{https://www.statforbiology.com/2020/stat_seedgermination_ht1step/}
}
\author{
Andrea Onofri
}

\examples{
data(barley)
barleyS <- subset(barley, Temp<40)
modTTEM <- drmte(nSeeds ~ timeBef + timeAf + Temp, data = barley,
               fct = TTEM())
summary(modTTEM)
}
\keyword{ seed germination }
\keyword{ thermal-time model }
\keyword{ time-to-event model }
