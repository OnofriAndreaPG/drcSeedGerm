\name{phalaris}
\alias{phalaris}
\docType{data}
\title{
A series of germination assays with Phalaris minor
}
\description{
The germination of Phalaris minor was assayed at 6 mater potential levels and 6 temperature levels, on four replicated Petri dishes at each temperature and water potential. Twenty-five seeds per Petri dish were used. Inspections were made daily for twenty days.
}
\usage{data("phalaris")}
\format{
  A data frame with 3024 observations on the following 9 variables.
  \describe{
    \item{\code{temp}}{a numeric vector with the temperature level}
    \item{\code{water}}{a numeric vector with the water potential level}
    \item{\code{Dish}}{a numeric vector with the coding for Petri dishes}
    \item{\code{timeBef}}{a numeric vector, with the start time for each inspection interval}
    \item{\code{timeAf}}{a numeric vector, with the end time for each inspection interval}
    \item{\code{nViable}}{a numeric vector, with the number of viable seeds per dish, at the beginning of the assay}
    \item{\code{nSeeds}}{a numeric vector, with the number of germinated seeds during each inspection interval}
    \item{\code{nCum}}{a numeric vector, with the cumulative number of germinated seeds at each assessment time}
    \item{\code{propCum}}{a numeric vector, the cumulative proportion of germinated seeds at each assessment time}
  }
}
\details{
The variable 'timeAf' contains the value 'Inf' (Infinity), that corresponds to the seeds which did not germinate during the assay, for which the germination times might be comprised from the last assessment time to infinity.
}
\source{
Mesgaran, MB, A Onofri, HR Mashhadi, RD Cousens (2017) Water availability shifts the optimal temperatures for seed germination: A modelling approach. Ecological Modelling 351, 87–95
}
\references{
Mesgaran, MB, A Onofri, HR Mashhadi, RD Cousens (2017) Water availability shifts the optimal temperatures for seed germination: A modelling approach. Ecological Modelling 351, 87–95
}
\author{
Mohsen Mesgaran
}

\examples{
data(phalaris)
}
\keyword{datasets}
