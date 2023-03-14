\name{HTE2}
\alias{HTE2}
\alias{HTE2.fun}

\title{
Hydrotime model with log-logistic distribution of germination time (Onofri et al., 2018)
}

\description{
This model relates the time-course of the proportion of germinated seeds to the water potential in the substrate. It is based on a truncated log-logistic distribution of germination time:


\deqn{P(t) = \frac{d}{1 + exp\left[ b (\log(x) - \log(e)\right]}}

where two of the three usual parameters ('d' and 'e') are expressed as functions of water potential (\eqn{\Psi}). In this function, the two submodels are: (1) for the parameter 'd', we used a shifted exponential function:

\deqn{d = G \left[ 1 - \exp \left( \frac{ \Psi - \Psi_b }{\sigma_{\Psi_b}} \right) \right]}

while, (2) for the parameter 'e' we considered that its inverse corresponds to the median Germination Rate within the population (i.e. \eqn{1/e = GR_{50}}) and modelled this latter parameter as:

\deqn{GR_{50} = \frac{\left[\Psi - \Psi_b \right]^2}{\theta_H} }

The difference with the 'HTE1()' function is that, in this case, the relationship between GR50 and water potential is not linear, but curvilinear (convex down).
The 'HTE2.fun()' is a generic function, which can be used for plotting or other applications, while the 'HTE2()' function is meant to be used for model fitting with the 'drmte()' function in the 'drcte()' package.
}

\details{
No more detail, at the moment.
}
\usage{
HTE2(fixed = c(NA, NA, NA, NA, NA),
                  names = c("G", "Psib", "sigmaPsib", "thetaH", "b"))
HTE2.fun(time, Psi, G, Psib, sigmaPsib, thetaH, b)
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
The 'HTE2.fun()' function returns the proportion of germinated seeds, for any given values of time and water potential in the substrate, depending on model parameters G, Psib, sigmaPsib, thetaH and b. The 'HTE2()' function returns a list containing the nonlinear function, the self starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}

\references{
Onofri, A., Benincasa, P., Mesgaran, M.B., Ritz, C., 2018. Hydrothermal-time-to-event models for seed germination. European Journal of Agronomy 101, 129–139.
\url{https://www.statforbiology.com/2020/stat_seedgermination_ht1step/}
}
\author{
Andrea Onofri
}

\examples{
# Fitting model
data(rape)
modHTE2 <- drmte( nSeeds ~ timeBef + timeAf + Psi,
            data=rape, fct=HTE2())
summary(modHTE2)
}
\keyword{ seed germination }
\keyword{ hydro-time model }
\keyword{ time-to-event model }
