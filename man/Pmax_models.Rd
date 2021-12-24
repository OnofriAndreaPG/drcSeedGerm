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
}
%- maybe also 'usage' for other objects documented here.
\details{
This R function is meant to be used with the 'drm()' function, within the 'drc' package.}
\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. Defaults to NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). The default is reasonable.}
}
\value{
The value returned is a list containing the nonlinear function, the self starter function and the parameter names.
}
\references{
\url{https://www.statforbiology.com/2021/stat_seedgermination_htt2step/#one-threshold-exponential}
}
\author{
Andrea Onofri
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
library(drcte)
# Pmax vs Psi (shifted exponential)
Psi <- seq(-2.2, 0, by = 0.2)
Pmax <- c(0, 0, 0.076, 0.413, 0.514, 0.643, 0.712,
          0.832, 0.865, 0.849, 0.89, 0.90)
mod <- drm(Pmax ~ Psi, fct = PmaxPsi1())
summary(mod)

#Pmax vs Psi (shifted exponential, with asymptote)
Psi <- seq(-2.2, 0, by = 0.2)
Pmax <- c(0, 0, 0.076, 0.413, 0.514, 0.643, 0.712,
          0.832, 0.865, 0.849, 0.89, 0.90)
mod <- drm(Pmax ~ Psi, fct = PmaxPsi1(fixed = c(NA, NA, 1)))
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
