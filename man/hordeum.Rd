\name{hordeum}
\alias{hordeum}
\docType{data}
\title{
Germination of Hordeum spontaneum at different temperatures and water potentials
}
\description{
This dataset was obtained from a germination assay with four replicated Petri dishes with 20 seeds, tested at six different water potential levels (0, -0.3, -0.6, -0.9, -1.2 and -1.5 MPa). Osmotic potentials were produced using variable amount of polyethylene glycol (PEG, molecular weight 8000) adjusted for the temperature level. Petri dishes were incubated at six constant temperature levels (8, 12, 16, 20, 24 and 28 °C), under a photoperiod of 12 h. Germinated seeds (radicle protrusion > 3 mm) were counted and removed daily for 20 days.
}
\usage{data("hordeum")}
\format{
  A data frame with 3024 observations on the following 8 variables.
  \describe{
    \item{\code{temp}}{a numeric vector: temperature level}
    \item{\code{water}}{a numeric vector: water potential level}
    \item{\code{Dish}}{a numeric vector: code for Petri dishes}
    \item{\code{timeBef}}{a numeric vector: beginning of scoring interval}
    \item{\code{timeAf}}{a numeric vector: end of scoring interval}
    \item{\code{nViable}}{a numeric vector: number of viable seeds at the beginning of assay, in each dish}
    \item{\code{nSeeds}}{a numeric vector: number of germinated seeds, between timeBef and timeAf}
    \item{\code{nCum}}{a numeric vector: cumulative number of germinated seeds at timeAf}
  }
}
\details{
This dataset was analysed in the time-to-event framework in Onofri et al (2018). See Example 2.
}
\source{
Mesgaran, MB, A Onofri, HR Mashhadi, RD Cousens (2017) Water availability shifts the optimal temperatures for seed germination: A modelling approach. Ecological Modelling 351:87–95
}
\references{
Mesgaran, MB, A Onofri, HR Mashhadi, RD Cousens (2017) Water availability shifts the optimal temperatures for seed germination: A modelling approach. Ecological Modelling 351:87–95

Onofri, A, P Benincasa, MB Mesgaran, C Ritz (2018) Hydrothermal-time-to-event models for seed germination. European Journal of Agronomy 101:129–139

}
\examples{
# Fitting a hydrotime model
gmod <- drm(propCum ~ timeAf + Psi, fct=HT(), data=rape)
summary(gmod)
jackGroupSE(gmod, rape, rape$Dish)
}
\keyword{seed germination, hydro-time models}
