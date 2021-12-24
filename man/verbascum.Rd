\name{verbascum}
\alias{verbascum}
\docType{data}
\title{
Germination time-course for three species of the Verbascum genus
}
\description{
This dataset relates to a germination assay with three species of Verbascum, four replicated Petri dishes and 25 seeds per Petri dish (12 dishes in all). Inspections were made daily for 15 days.
}
\usage{data("verbascum")}
\format{
  A data frame with 192 observations on the following 6 variables.
  \describe{
    \item{\code{Dish}}{a numeric vector with the coding for Petri dishes}
    \item{\code{Species}}{a factor with the name of the species}
    \item{\code{timeBef}}{a numeric vector with the start date for the each assessment interval}
    \item{\code{timeAf}}{a numeric vector with the end date for the each assessment interval. It contains the value Inf (Infinity), that indicates the final interval of time, after the final assessment time}
    \item{\code{nSeeds}}{a numeric vector with the number of germinated seeds in each interval. Relating to final interval of time, from the last assessment date to Infinity, 'nSeeds' contains the final number of ungerminated seeds}
    \item{\code{nCum}}{a numeric vector with the cumulative number of germinated seeds in each interval of time}
   }
}
\references{
Catara, S., Cristaudo, A., Gualtieri, A., Galesi, R., Impelluso, C., Onofri, A., 2016. Threshold temperatures for seed germination in nine species of Verbascum (Scrophulariaceae). Seed Science Research 26, 30–46.
}
\examples{
data(verbascum)
## maybe str(verbascum) ; plot(verbascum) ...
}
\keyword{datasets}
