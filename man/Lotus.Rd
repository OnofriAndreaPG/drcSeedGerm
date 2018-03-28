\name{lotus}
\alias{lotus}
\docType{data}
\title{
Germination assay with Lotus ornithopodioides
}
\description{
Seeds of L. ornithopodioides were harvested in three stages after anthesis and the time course of their germination was evaluated in Petri dishes. The dataset is organised as required by time-to-event analyses with the drc package.
}
\usage{data("lotus")}
\format{
  A data frame with 192 observations on the following 7 variables.
  \describe{
    \item{\code{Stage}}{a character vector coding for the maturation stage}
    \item{\code{Dish}}{a numeric vector coding for the Petri Dish}
    \item{\code{timeBef}}{a numeric vector: beginning of scoring interval}
    \item{\code{timeAf}}{a numeric vector: end of scoring interval}
    \item{\code{count}}{a numeric vector: count of germinated seeds in each scoring interval}
    \item{\code{nCum}}{a numeric vector: cumulative count of germinated seeds}
    \item{\code{propCum}}{a numeric vector: cumulative proportions of germinated seeds}
  }
}
\details{
There were 25 seeds per Petri dish. Seeds of Lotus ornithopodioides were harvested at the following stages (A) Green fleshy pods with brilliant gree sedds (27 Days After Anthesis; DAA), (B) green pergamenous pods with brown soft seeds (35 DAA) and (C) green pergamenous pods with brown seeds and moderately hard (41 DAA).

}
\source{
Data are kindly provided by Fabio Gresta (University of Reggio Calabria, Italy).}
\references{
Gresta, F., Avola, G., Onofri, A., Anastasi, U., Cristaudo, A., 2011. When Does Hard Coat Impose Dormancy in Legume Seeds? Lotus and Scorpiurus Case Study. Crop Science 51, 1739–1747.
}
\examples{
data(lotus)
library(drc)
#A constraint was put on the second curve (Stage B), in order to avoid an estimate
#higher than 1.
mod <- drm(count ~ timeBef + timeAf, data=lotus, fct=LL.3(),
           curveid=Stage, type="event", upperl = c(NA, 1, NA))
summary(mod)
}
\keyword{time-to-event analysis}
