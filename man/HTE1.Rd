\name{HTE1}
\alias{HTE1}
\alias{HTE1.fun}

\title{
Hydrotime model with log-logistic distribution of germination time (Onofri et al., 2018)
}

\description{
This model relates the time-course of the proportion of germinated seeds to the water potential in the substrate. It is based on a truncated log-logistic distribution of germination time:


\deqn{P(t) = \frac{d}{1 + exp\left[ b (\log(x) - \log(e)\right]}}

where two of the three usual parameters ('d' and 'e') are expressed as functions of water potential (\eqn{\Psi}). In this function, the two submodels are: (1) for the parameter 'd', we used a shifted exponential function:

\deqn{d = G \left[ 1 - \exp \left( \frac{ \Psi - \Psi_b }{\sigma_{\Psi_b}} \right) \right]}

while, (2) for the parameter 'e' we considered that its inverse corresponds to the median Germination Rate within the population (i.e. \eqn{1/e = GR_{50}}) and modelled this latter parameter as:

\deqn{GR_{50} = \frac{\Psi - \Psi_{b}}{\theta_H} }


The 'HTE1.fun()' is a generic function, which can be used for plotting or other applications, while the 'HTE1()' function is meant to be used for model fitting with the 'drmte()' function in the 'drcte()' package.
}
\details{
The detail of this time-to-event model and the meaning of parameters are described in Onofri et al. (2018). See Table 2, where 'HTE1()' is abbreviated as HTE.
}
\usage{
HTE1(fixed = c(NA, NA, NA, NA, NA),
                  names = c("G", "Psib", "sigmaPsib", "thetaH", "b"))
HTE1.fun(time, Psi, G, Psib, sigmaPsib, thetaH, b)
}
\arguments{
These functions have the following arguments:
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed. At the moment, only the parameter G can be held fixed.}
 \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{time}{time}
  \item{Psi}{water potential in the substrate}
  \item{G}{maximum germination capability of the seed lot}
  \item{Psib}{base water potential}
  \item{sigmaPsib}{standard deviation for the base water potential within the seed lot}
  \item{thetaH}{hydro-time constant}
  \item{b}{shape parameter for the cumulative distribution function of germination time}
}
\value{
The 'HTE1.fun()' function returns the proportion of germinated seeds, for any given values of time and water potential in the substrate. The 'HTE1()' function returns a list containing the nonlinear function, the self-starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}

\references{
Onofri, A., Benincasa, P., Mesgaran, M.B., Ritz, C., 2018. Hydrothermal-time-to-event models for seed germination. European Journal of Agronomy 101, 129–139.
\url{https://www.statforbiology.com/2020/stat_seedgermination_ht1step/}
}
\author{
Andrea Onofri
}

\examples{
data(rape)
modHTE <- drmte( nSeeds ~ timeBef + timeAf + Psi,
            data=rape, fct=HTE1())
summary(modHTE)
}
\keyword{ seed germination }
\keyword{ hydro-time model }
\keyword{ time-to-event model }
