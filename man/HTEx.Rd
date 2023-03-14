\name{HTex}
\alias{HTex}
\alias{HTex.fun}

\title{
Hydrotime model with Type II Extreme Value distribution of base water potential (Mesgaran et al., 2013)
}

\description{
This model relates the time-course of the proportion of germinated seeds to the water potential in the substrate. It is based on Type II Extreme Value distribution of base water potential within the seed lot. The equation is:

\deqn{ P(t) = 1 - \exp \left\{ - \exp \left[ \frac{\Psi - (\theta _H / t ) - \mu }{\sigma } \right] \right\} }

In contrast to other hydrotime models (e.g., 'HTE1()', 'HTE2()', 'HTE3()', describing the distribution of germination time), this models describes the distribution of base water potential within the population, while the distribution of germination times is indirectly modelled, but it is not, in itself, gaussian (you see that 't' is at the denominator).
The 'HTex.fun()' is a generic function, which can be used for plotting or other applications, while the 'HTex()' function is meant to be used for model fitting with the 'drmte()' function in the 'drcte()' package.
}

\details{
The detail of this time-to-event model and the meaning of parameters are described in Mesgaran et al. (2013).
}

\usage{
HTex()
HTex.fun(time, Psi, thetaH, mu, sigma)
}
\arguments{
The 'HTex()' function has no arguments. The 'HTex.fun()' function has the following arguments:
  \item{time}{time}
  \item{Psi}{water potential in the substrate}
  \item{thetaH}{hydro-time constant}
  \item{mu}{location parameter for the Type II Extreme Value distribution}
  \item{sigma}{scale parameter for the Type II Extreme Value distribution}
}
\value{
The 'HTex.fun()' function returns the proportion of germinated seeds, for any given values of time and water potential in the substrate. The 'HTex()' function returns a list containing the nonlinear function, the self-starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}
\references{
Mesgaran, M.B., Mashhadi, H.R., Alizadeh, H., Hunt, J., Young, K.R., Cousens, R.D., 2013. Importance of distribution function selection for hydrothermal time models of seed germination. Weed Research 53, 89–101.
}
\author{
Andrea Onofri
}
\note{
This function is for use with the function 'drmte()'
}

\examples{
data(rape)
modg <- drm( propCum ~ timeAf + Psi, fct=HTex(), data=rape)
summary(modg)
}
\keyword{ seed germination }
\keyword{ hydro-time model }
\keyword{ time-to-event model }
