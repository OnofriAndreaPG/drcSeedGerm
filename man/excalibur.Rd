\name{excalibur}
\alias{excalibur}
\docType{data}
\title{
Relationship between germination rate and water potential in oilseed rape (var. Excalibur)
}
\description{
This files describes the relationship between germination rate and water potential in the substrate for seeds of oilseed rape (var. Excalibur). Three germination percentiles are considered for germination rate, i.e. GR10, GR30 and GR50.
}
\usage{data("excalibur")}
\format{
  A data frame with 27 observations on the following 5 variables.
  \describe{
    \item{\code{Perc}}{a numeric vector: the germination percentile}
    \item{\code{Psi}}{a numeric vector: water potential in the substrate (in MPa)}
    \item{\code{Tg}}{a numeric vector: germination time in days}
    \item{\code{SE}}{a numeric vector: standard errors for germination times (in days)}
    \item{\code{GR}}{a numeric vector: germination rates in d^{-1}}
  }
}
\source{
Pace, R., Benincasa, P., Ghanem, M.E., Quinet, M., Lutts, S., 2012. GERMINATION OF UNTREATED AND PRIMED SEEDS IN RAPESEED (BRASSICA NAPUS VAR OLEIFERA DEL.) UNDER SALINITY AND LOW MATRIC POTENTIAL. Experimental Agriculture 48, 238–251.

}
\references{
Pace, R., Benincasa, P., Ghanem, M.E., Quinet, M., Lutts, S., 2012. Germination of untreated and primed seeds in rapeseed (Brassica napus var. oleifera Del.) under salinity and low matric potential. Experimental Agriculture 48, 238–251.
}
\examples{
data(excalibur)
head(excalibur)
library(drc)

modGR1 <- drm(GR ~ Psi, fct=GRPsiLin(), data=excalibur, curveid=Perc)
summary(modGR1)
plot(modGR1, log="", legendPos=c(-0.9, 1))

modGR2 <- drm(GR ~ Psi, fct=GRPsiPol2(), data=excalibur, curveid=Perc)
summary(modGR2)
plot(modGR2, log="", legendPos=c(-0.9, 1))
}
\keyword{datasets, seed germination, hydrotime models}
