\name{rape2G}
\alias{rape2G}
\docType{data}
\title{
Germination data from an assay of rapeseed at decreasing water potential levels
}
\description{
This files describes the germination of rapeseed (cv. Excalibur and Toccata) at different water potential levels in the substrate.
}
\usage{data("rape")}
\format{
  A data frame with 294 observations on the following 8 variables.
  \describe{
    \item{\code{Psi}}{a factor: rape genotype}
    \item{\code{Psi}}{a numeric vector: waterpotential level}
    \item{\code{Dish}}{a numeric vector: code for Petri dishes}
    \item{\code{timeBef}}{a numeric vector: beginning of scoring interval}
    \item{\code{timeAf}}{a numeric vector: end of scoring interval}
    \item{\code{nSeeds}}{a numeric vector: number of germinated seeds, between timeBef and timeAf}
    \item{\code{nCum}}{a numeric vector: cumulative number of germinated seeds at timeAf}
    \item{\code{propCum}}{a numeric vector: cumulative proportion of germinated seeds at timeAf}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Pace, R., Benincasa, P., Ghanem, M.E., Quinet, M., Lutts, S., 2012. Germination of untreated and primed seeds in rapeseed (Brassica napus var. oleifera Del.) under salinity and low matric potential. Experimental Agriculture 48, 238–251.
}
\references{
Pace, R., Benincasa, P., Ghanem, M.E., Quinet, M., Lutts, S., 2012. Germination of untreated and primed seeds in rapeseed (Brassica napus var. oleifera Del.) under salinity and low matric potential. Experimental Agriculture 48, 238–251.}
\examples{
#Fitting a hydrotime-to-event model
modHTE <- drm( nSeeds ~ timeBef + timeAf + Psi,
            data=rape, fct=HTE1(), type="event", 
            curveid = CV)
summary(modHTE)
jackGroupSE(modHTE, data = rape, cluster = rape$Dish)
}
\keyword{seed germination, hydro-time models}
