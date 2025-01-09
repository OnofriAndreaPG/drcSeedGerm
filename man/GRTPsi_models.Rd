\name{GR_TPsi_models}
\alias{GRTPsi.BS}
\alias{GRTPsi.BS.fun}
\alias{GRTPsi.M}
\alias{GRTPsi.M.fun}

\title{Hydro-Thermal-time models based on the germination rate}
\description{
These models are used to describe the germination rate of a seed, depending on the environmental temperature and water potential.}

\usage{
GRTPsi.M()
GRTPsi.M.fun(Temp, Psi, k, Tb, ThetaHT, Psib)
GRTPsi.BS()
GRTPsi.BS.fun(Temp, Psi, k, Tb, To, ThetaHT, Psib)
}
\details{
The 'GRT.M.fun()' is a generic R function, while the GRT.M() function is meant to be used with the 'drm()' function, within the 'drc' package.}
\arguments{
  \item{Temp}{Temperature variable}
  \item{Psi}{Water potential variable}
  \item{Psib}{Base water potential variable}
  \item{k}{regression parameter}
  \item{Tb}{base temperature}
  \item{To}{optimal temperature}
  \item{ThetaHT}{Hydro-thermal-time parameter}
}
\value{
The 'GRT.M.fun()' functions returns a vector of responses, for given values of temperature, Tc, Tb and ThetaH. The GRT.M() function returns a list containing the nonlinear function, the self starter function, the parameter names and other items which are internally used by the 'drm()' function.
}
\references{
\url{https://www.statforbiology.com/2023/stat_drcte_12-htt2step/}
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
