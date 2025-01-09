\name{Thermal-Time_models}
\alias{TTEM}
\alias{TTEM.fun}
\alias{TTERF}
\alias{TTERF.fun}
\alias{TTERFc}
\alias{TTERFc.fun}

\title{
Thermal-time models with log-logistic distribution of germination time (Onofri et al., 2018)
}
\description{
These models relates the time-course of germinations to the environmental temperature and they are based on a truncated log-logistic distribution of germination time:

\deqn{P(t) = \frac{d}{1 + exp\left[ b (\log(x) - \log(e)\right]}}

where the usual parameters ('b', 'd' and 'e') are expressed as functions of wtemperature (\eqn{T}). In the function 'TTEM()', we implemented the following submodels: (1) for the parameter 'd', we implemented the following submodels:

\deqn{d = G \, \left[ 1 - \exp \left( - \frac{ T_c - T }{\sigma_{T_c}} \right) \right]}

\deqn{ 1/[e(T)] = GR_{50}(T) = \frac{T - T_b }{\theta_T} \left[1 - \frac{T - T_b}{T_c - T_b}\right] }

while 'b' was regarded as constant and independent from temperature.

In the 'TTERF()' function, the last submodel was modified, according to Rowse and Finch-Savage (2003):

\deqn{
1/[e(T)] = GR_{50}(T) = \left\{ {\begin{array}{ll}
\frac{T - T_b}{\theta_T} & \textrm{if} \,\,\, T_b < T < T_d \\
\frac{T - T_b}{\theta_T} \left[ 1 - \frac{T - T_d}{T_c - T_d}  \right] & \textrm{if} \,\,\, T_d < T < T_c \\
0 & \textrm{if} \,\,\, T \leq T_b \,\,\, or \,\,\, T \geq T_c
\end{array}} \right.
}

In the In the 'TTERFc()' function a further submodel was introduced, to model the effect of temperature in the shape parameter 'b':

\deqn{
\sigma(T) = \frac{1}{b} = \frac{1}{b_0} + s (T - T_b)
}

The 'TTEM.fun()', 'TTERF.fun()' and 'TTERFc.fun()' are generic functions, which can be used for plotting or other applications, while the 'TTEM()', 'TTERF()' and 'TTERFc()' functions are meant to be used for model fitting with the 'drmte()' function in the 'drcte()' package.
}
\usage{
TTEM()
TTEM.fun(time, Temp, G, Tc, sigmaTc, Tb, ThetaT, b)
TTERF()
TTERF.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b)
TTERFc()
TTERFc.fun(time, Temp, G, Tc, sigmaTc, Td, Tb, ThetaT, b0, s)
}
\details{
The detail of these functions are described in Onofri et al. (2018).
}
\arguments{
  The TTEM(), TTERF() and TTERFc() functions have no arguments. The TTEM.fun(), TTERF.fun() and TTERFc.fun() functions have the following arguments:
  \item{time}{time}
  \item{Temp}{temperature}
  \item{G}{maximum germination capability of the seed lot}
  \item{Tc}{ceiling temperature}
  \item{sigmaTc}{standard deviation for the ceiling temperature within the seed lot}
  \item{Tb}{base temperature}
  \item{Td}{close to optimal temperature}
  \item{ThetaT}{Thermal-time parameter}
  \item{b}{shape parameter for the cumulative distribution function of germination time}
  \item{b0}{shape parameter for the cumulative distribution function of germination time at T = Tb}
  \item{s}{parameter for the effect of temperature on the scale parameter for the log-logistic distribution of germination times}
}
\value{
The 'TTEM.fun()', 'TTERF.fun()' and 'TTERFc.fun()' functions return the proportion of germinated seeds, for given values of time and temperature. The 'TTEM()', 'TTERF()' and 'TTERFc()' functions return a list containing the nonlinear function, the self starter function, the parameter names and other items which are internally used by the 'drmte()' function.
}
\references{
Onofri, A., Benincasa, P., Mesgaran, M.B., Ritz, C., 2018. Hydrothermal-time-to-event models for seed germination. European Journal of Agronomy 101, 129–139.
Rowse, H.R., Finch-Savage, W.E., 2003. Hydrothermal threshold models can describe the germination response of carrot (Daucus carota) and onion (Allium cepa) seed populations across both sub- and supra-optimal temperatures. New Phytologist 158, 101–108.
\url{https://www.statforbiology.com/2023/stat_drcte_11-examplette/}
}
\author{
Andrea Onofri
}

\examples{
data(barley)
modTTE <- drmte(nSeeds ~ timeBef + timeAf + Temp, data = barley,
               fct = TTERF())
summary(modTTE)
}
\keyword{ seed germination }
\keyword{ thermal-time model }
\keyword{ time-to-event model }
