\name{festuca}
\alias{festuca}
\docType{data}
\title{
Relationship between germination rate and water potential in Festuca arundinacea L.
}
\description{
This files describes the relationship between germination rate and Festuca arundinacea. Three germination percentiles are considered for germination rate, i.e. GR10, GR30 and GR50.
}
\usage{data("festuca")}
\format{
  A data frame with 36 observations on the following 3 variables.
  \describe{
    \item{\code{g}}{a numeric vector: the germination percentile}
    \item{\code{Psi}}{a numeric vector: water potential in the substrate (in MPa)}
    \item{\code{GR}}{a numeric vector: germination rates in 1/d}
  }
}
\source{
TEI F, BENINCASA P and CIRICIOFOLO E (2001) Effetto del potenziale idrico e della temperatura sulla germinazione di alcune specie graminacee da tappeto erboso. In: Atti XXXIV Convegno della Società Italiana di Agronomia, Pisa, Italy, 200–201.
}
\references{
TEI F, BENINCASA P and CIRICIOFOLO E (2001) Effetto del potenziale idrico e della temperatura sulla germinazione di alcune specie graminacee da tappeto erboso. In: Atti XXXIV Convegno della Società Italiana di Agronomia, Pisa, Italy, 200–201.
}
\examples{
data(festuca)
modGR1 <- drm(GR ~ Psi, fct=GRPsiLin(),
            data=festuca, curveid=g)
summary(modGR1)
}
\keyword{datasets}
\keyword{seed germination}
\keyword{hydrotime models}
