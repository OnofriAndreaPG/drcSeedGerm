\name{GR_models}
\alias{GRT.M}
\alias{GRT.M.fun}
\title{Modelling the effect of temperature on the germination rate of a seed}
\description{
These models are used to describe the germination rate of a seed, depending on the environmental temperature.}

\usage{
GRT.M()
GRT.M.fun(predictor, Tc, Tb, ThetaT)
}
\details{
The 'GRT.M.fun()' is a generic R function, while the GRT.M() function is meant to be used with the 'drm()' function, within the 'drc' package.}
\arguments{
  \item{predictor}{Temperature variable}
  \item{Tc}{ceiling temperature}
  \item{Tb}{base temperature}
  \item{ThetaT}{Hydro-time parameter}
}
\value{
The 'GRT.M.fun()' functions returns a vector of responses, for given values of temperature, Tc, Tb and ThetaH. The GRT.M() function returns a list containing the nonlinear function, the self starter function, the parameter names and other items which are internally used by the 'drm()' function.
}
\references{
\url{https://www.statforbiology.com/2021/stat_seedgermination_htt2step/#polynomial-model}
}
\author{
Andrea Onofri
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
library(drcte)
Tval <- c(2, 5, 10, 15, 20, 25, 30, 35, 40)
GR <- c(0, 0, 0.209, 0.435, 0.759, 0.821, 0.417, 0.145, 0)
modM <- drm(GR ~ Tval, fct = GRT.M())
plot(modM, log="", xlim = c(0, 40), ylim=c(0,1.2),
     legendPos = c(5, 1.0), xlab = "Temperature (°C)")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{seed germination}
\keyword{germination rate}
