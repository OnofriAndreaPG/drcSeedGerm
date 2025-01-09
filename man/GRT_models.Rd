\name{GR_T_models}
\alias{GRT.GH}
\alias{GRT.GH.fun}
\alias{GRT.GH2}
\alias{GRT.GH2.fun}
\alias{GRT.YL}
\alias{GRT.YL.fun}
\alias{GRT.BS}
\alias{GRT.BS.fun}
\alias{GRT.BSb}
\alias{GRT.BSb.fun}
\alias{GRT.Ex}
\alias{GRT.Ex.fun}
\alias{GRT.Exb}
\alias{GRT.Exb.fun}
\alias{GRT.M}
\alias{GRT.M.fun}
\alias{GRT.Mb}
\alias{GRT.Mb.fun}
\alias{GRT.RF}
\alias{GRT.RF.fun}
\alias{GRT.RFb}
\alias{GRT.RFb.fun}

\title{Thermal-time models based on the germination rate}
\description{
These models are used to describe the germination rate of a seed, depending on the environmental temperature.}

\usage{
GRT.GH()
GRT.GH2()
GRT.YL()
GRT.BS()
GRT.BSb()
GRT.Ex()
GRT.Exb()
GRT.M()
GRT.Mb()
GRT.RF()
GRT.RFb()
GRT.GH.fun(Temp, Tb, ThetaT)
GRT.GH2.fun(Temp, Tb, beta)
GRT.YL.fun
GRT.BS.fun(Temp, k, Tb, To, ThetaT)
GRT.BSb.fun(Temp, Tc, Tb, To, ThetaT)
GRT.Ex.fun(Temp, k, Tb, Tc, ThetaT)
GRT.Exb.fun(Temp, k, Tb, Tc, ThetaT)
GRT.M.fun(Temp, k, Tb, ThetaT)
GRT.Mb.fun(Temp, Tb, Tc, ThetaT)
GRT.RF.fun(Temp, k, Tb, Td, ThetaT)
GRT.RFb.fun(Temp, Tc, Tb, Td, ThetaT)
}
\details{
All these functions are named according to this rule: 'GRT' (Germination Rate Temperature), followd by the 'function name' (e.g., BS, RF, M, Ex, YL). The R functions 'GR.funName().fun' are generic R function, that are meant to be used for general purposes, such as plotting or predicting. The corresponding 'GR.funName()' (without the '.fun' ending) are meant to be used for fitting with the 'drm()' function, within the 'drc' package.}
\arguments{
  The 'GR.funName()' functions have no arguments. The general purpose 'GR.funName.fun()' functions have some of the following arguments (depending on function):
  \item{Temp}{Temperature variable}
  \item{k}{regression parameter}
  \item{beta}{regression parameter}
  \item{Tc}{ceiling temperature}
  \item{Tb}{base temperature}
  \item{To}{optimal temperature}
  \item{Td}{close-to-optimal temperature (Rowse-Fintch-Savage equation)}
  \item{ThetaT}{Hydro-time parameter}
}
\value{
The 'GRT.funName.fun()' functions return a vector of responses, for given values of temperature and parameters. The 'GRT.funName()' functions return a list containing the nonlinear function, the self starter function, the parameter names and other items which are internally used by the 'drm()' function.
}
\references{
\url{https://www.statforbiology.com/2023/stat_drcte_12-htt2step/}
}
\author{
Andrea Onofri
}

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
