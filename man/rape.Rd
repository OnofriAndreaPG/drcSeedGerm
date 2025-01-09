\name{rape}
\alias{rape}
\docType{data}
\title{
Germination data from an assay of rapeseed at decreasing water potential levels
}
\description{
This files describes the germination of rapeseed (cv. Excalibur) at different water potential levels in the substrate.
}
\usage{data("rape")}
\format{
  A data frame with 294 observations on the following 7 variables.
  \describe{
    \item{\code{Psi}}{a numeric vector: waterpotential level}
    \item{\code{Dish}}{a numeric vector: code for Petri dishes}
    \item{\code{timeBef}}{a numeric vector: beginning of scoring interval}
    \item{\code{timeAf}}{a numeric vector: end of scoring interval}
    \item{\code{nSeeds}}{a numeric vector: number of germinated seeds, between timeBef and timeAf}
    \item{\code{nCum}}{a numeric vector: cumulative number of germinated seeds at timeAf}
    \item{\code{propCum}}{a numeric vector: cumulative proportion of germinated seeds at timeAf}
  }
}
\source{
Pace, R., Benincasa, P., Ghanem, M.E., Quinet, M., Lutts, S., 2012. Germination of untreated and primed seeds in rapeseed (Brassica napus var. oleifera Del.) under salinity and low matric potential. Experimental Agriculture 48, 238–251.
}
\references{
Pace, R., Benincasa, P., Ghanem, M.E., Quinet, M., Lutts, S., 2012. Germination of untreated and primed seeds in rapeseed (Brassica napus var. oleifera Del.) under salinity and low matric potential. Experimental Agriculture 48, 238–251.}
\examples{
#Fitting a hydrotime model
library(drcte)
data(rape)
gmod <- drmte(nSeeds ~ timeBef + timeAf + Psi,
                 data = rape, fct = HTE1())
summary(gmod)
}
\keyword{datasets}
\keyword{seed germination}
\keyword{hydrotime models}
