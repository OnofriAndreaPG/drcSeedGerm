\name{Pmax_models}
\alias{PmaxPsi1}
\alias{PmaxPsi1.fun}
\alias{PmaxT1}
\alias{PmaxT1.fun}
\title{Effect of environmental water potential or temperature on the germination capability of a seed lot
}
\description{
These models are used to describe the germination capability of a seed lot, depending on the environmental water potential or temperature.}

\usage{
PmaxPsi1(fixed = c(NA, NA, NA), names = c("G", "Psib", "sigma"))
PmaxT1(fixed = c(NA, NA, NA), names = c("G", "Tc", "sigmaTc"))
PmaxPsi1.fun(Psi, G, Psib, sigma)
PmaxT1.fun(Temp, G, Tc, sigmaTc)
}
\details{
The R functions 'PmaxPsi1()' and 'PmaxT1()' are meant to be used for model fitting with the 'drm()' function, within the 'drc' package. On the other hand, 'PmaxPsi1.fun()' and 'PmaxT1.fun()' are general purpose functions, to bu used for plotting or other applications.}
\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. Defaults to NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). The default is reasonable.}
  \item{Psi}{Water potential level}
  \item{Temp}{Temperature level}
  \item{G}{Maximum proportion of germinated seeds}
  \item{Psib}{Base water potential level}
  \item{sigma}{Variability of base water potential within the seed lot}
  \item{Tc}{Ceiling temperature level}
  \item{sigmaTc}{Variability of ceiling temperature level within the population}
}
\value{
'PmaxPsi1()' and 'PmaxT1()' return a list containing the nonlinear function, the self starter function and the parameter names, that are internally used for model fitting. 'PmaxPsi1.fun()' and 'PmaxT1.fun()' return the maximum proportion of germinated seeds for any given level of temperature or water potential in the substrate.
}
\references{
\url{https://www.statforbiology.com/2021/stat_drcte_12_HTT2step}
}
\author{
Andrea Onofri
}

\examples{
library(drcte)
# Pmax vs Psi (shifted exponential)
Psi <- seq(-2.2, 0, by = 0.2)
Pmax <- c(0, 0, 0.076, 0.413, 0.514, 0.643, 0.712,
          0.832, 0.865, 0.849, 0.89, 0.90)
mod <- drm(Pmax ~ Psi, fct = PmaxPsi1())
summary(mod)

# Pmax vs Psi (shifted exponential, with asymptote)
Psi <- seq(-2.2, 0, by = 0.2)
Pmax <- c(0, 0, 0.076, 0.413, 0.514, 0.643, 0.712,
          0.832, 0.865, 0.849, 0.89, 0.90)
mod <- drm(Pmax ~ Psi, fct = PmaxPsi1(fixed = c(1, NA, NA)))
summary(mod)

# Pmax vs temperature
Tval <- c(0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5,
          20, 22.5, 25, 27.5, 30, 32.5, 35)
Pmax2 <- c(0.79, 0.81, 0.807, 0.776, 0.83,
           0.73, 0.744, 0.73, 0.828, 0.818,
           0.805, 0.706, 0.41, 0.002, 0)
mod2 <- drm(Pmax2 ~ Tval, fct = PmaxT1())
summary(mod2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{seed germination}
\keyword{capability of germination}
