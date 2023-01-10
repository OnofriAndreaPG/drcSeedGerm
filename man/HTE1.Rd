\name{HTE1}
\alias{HTE1}
\alias{HTE1.fun}

\title{
Hydro-time-to-event model for seed germination.
}
\description{
This model relates the time-course of the proportion of germinated seeds to the water potential in the substrate. It is based on four submodels: (1) a truncated log-logistic distribution of germination time, where the three usual parameters ('d', 'e' and 'b') are expressed as functions of water potential (Psi); (2) for the parameter 'd', the dependency on Psi is described by using a shifted exponential distribution; (3) for the parameter 'e', it is assumed that its inverse (i.e. 1/e = GR50) is linearly related to Psi; (4) for the 'b' parameter, it is assumed that it is unaffected by Psi. The 'HTE1.fun()' is a generic function, which can be used for plotting or other applications, while the HTE1() function is meant to be used to fit this model by the 'drmte()' function in the 'drcte()' package.
}
\usage{
HTE1()
HTE1.fun(time, Psi, G, Psib, sigmaPsib, thetaH, b)
}
\details{
The detail of this time-to-event model is described in Onofri et al. (2018). See Table 2, where this model is abbreviated as HTE.
}
\arguments{
  The HTE1() function has no arguments. The HTE1.fun() function has the following arguments:
  \item{time}{time}
  \item{Psi}{water potential in the substrate}
  \item{G}{maximum germination capability of the seed lot}
  \item{Psib}{base water potential}
  \item{sigmaPsib}{standard deviation for the base water potential within the seed lot}
  \item{thetaH}{hydro-time constant}
  \item{b}{shape parameter for the cumulative distribution function of germination time}
}
\value{
The 'HTE1.fun()' functions returns the proportion of germinated seeds, for given values of time and water potential in the substrate, depending on model parameters G, Psib, sigmaPsib, thetaH and b. The HTE1() function returns a list containing the nonlinear function, the self starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}
\references{
Onofri, A., Benincasa, P., Mesgaran, M.B., Ritz, C., 2018. Hydrothermal-time-to-event models for seed germination. European Journal of Agronomy 101, 129–139.
\url{https://www.statforbiology.com/2020/stat_seedgermination_ht1step/}
}
\author{
Andrea Onofri
}

\examples{
library(drcte)
data(rape)
modHTE <- drmte(nSeeds ~ timeBef + timeAf + Psi,
                data = rape, fct = HTE1())
summary(modHTE)
}
\keyword{ seed germination }
\keyword{ hydro-time model }
\keyword{ time-to-event model }
