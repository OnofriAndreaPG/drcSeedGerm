\name{HTnorm}
\alias{HTnorm}
\alias{HTnorm.fun}

\title{
Hydrotime model with normal distribution of base water potential (Bradford, 2002)
}

\description{
This model relates the time-course of the proportion of germinated seeds to the water potential in the substrate. It is based on a normal distribution of base water potential within the seed lot. The equation is:

\deqn{ P(t) = \Phi \left\{ \frac{\Psi - (\theta_H / t) -\Psi_b }{\sigma_{\Psi_b}} \right\} }

where \eqn{\Phi} is a gaussian cumulative distribution function for base water potential. In contrast to other hydrotime models (e.g., 'HTE1()', 'HTE2()', 'HTE3()', describing the distribution of germination time), this models describes the distribution of base water potential within the population, while the distribution of germination times is indirectly modelled, but it is not, in itself, gaussian (you see that 't' is at the denominator).
The 'HTNorm.fun()' is a generic function, which can be used for plotting or other applications, while the 'HTNorm()' function is meant to be used for model fitting with the 'drmte()' function in the 'drcte()' package.
}

\details{
The detail of this time-to-event model and the meaning of parameters are described in Bradford (2002).
}

\usage{
HTnorm()
HTnorm.fun(time, Psi, ThetaH, Psib50, sigmaPsib)
}
\arguments{
The 'HTNorm()' function has no arguments. The 'HTNorm.fun()' function has the following arguments:
  \item{time}{time}
  \item{Psi}{water potential in the substrate}
  \item{ThetaH}{hydro-time constant}
  \item{Psib50}{median base water potential}
  \item{sigmaPsib}{standard deviation for the base water potential within the seed lot}
}
\value{
The 'HTNorm.fun()' function returns the proportion of germinated seeds, for any given values of time and water potential in the substrate. The 'HTNorm()' function returns a list containing the nonlinear function, the self-starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}
\references{
Bradford, K.J., 2002. Applications of hydrothermal time to quantifying and modeling seed germination and dormancy. Weed Science 50, 248–260.
}
\author{
Andrea Onofri
}
\note{
This function is for use with the R function 'drmte()'
}

\examples{
data(rape)
# non linear regression (wrong)
modg <- drm( propCum ~ timeAf + Psi, fct=HTnorm(), data=rape)
# time-to-event regression (right)
modg <- drmte(nSeeds ~ timeBef + timeAf + Psi, fct=HTnorm(), data=rape)
summary(modg)
}
\keyword{ seed germination }
\keyword{ hydro-time model }
\keyword{ time-to-event model }
