\name{barley}
\alias{barley}
\docType{data}
\title{
A series of germination assays with barley
}
\description{
The germination of barley was assayed at 9 temperature levels, on three replicated Petri dishes at each temperature and 50 seeds per Petri dish. Inspections were made in several times after the beginning of the assay.
}
\usage{data("barley")}
\format{
  A data frame with 810 observations on the following 7 variables.
  \describe{
    \item{\code{Dish}}{a numeric vector with the coding for Petri dishes}
    \item{\code{Temp}}{a numeric vector with the temperature level}
    \item{\code{timeBef}}{a numeric vector, with the start time for each inspection interval}
    \item{\code{timeAf}}{a numeric vector, with the end time for each inspection interval}
    \item{\code{nSeeds}}{a numeric vector, with the number of germinated seeds at each inspection interval}
    \item{\code{nCum}}{a numeric vector, with the cumulative number of germinated seeds at each assessment time}
    \item{\code{propCum}}{a numeric vector, the cumulative proportion of germinated seeds at each assessment time}
  }
}
\details{
The variable 'timeAf' contains the value 'Inf' (Infinity), that corresponds to the seeds which did not germinate during the assay, for which we the germination time might be comprised from the last assessment time to infinity.
}
\source{
no reference yet
}
\references{
Onofri, A., Benincasa, P., Mesgaran, M.B., Ritz, C., 2018. Hydrothermal-time-to-event models for seed germination. European Journal of Agronomy 101, 129–139.}
\author{
Andrea Onofri
}

\examples{
data(barley)
## maybe str(barley) ; plot(barley) ...
}
\keyword{datasets}
