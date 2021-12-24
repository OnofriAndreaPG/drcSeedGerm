\name{alfalfaSG}
\alias{alfalfaSG}
\docType{data}
\title{
Field-book for a germination assay with alfalfa
}
\description{
The germination of alfalfa was assayed at 7 temperature levels, on three replicated Petri dishes per temperature and 100 seeds per Petri dish. Inspections were made in several times after the beginning of the assay.
}
\usage{data("alfalfaSG")}
\format{
  A data frame with 21 observations on the following variables.
  \describe{
    \item{\code{Dish}}{a numeric vector with the coding for Petri dishes}
    \item{\code{Temp}}{a numeric vector with the temperature level}
    \item{\code{Dish}}{a numeric vector: code of Petri dish}
    \item{\code{nViable}}{a numeric vector: number of viable seeds per dish}
    \item{\samp{1}}{a numeric vector: count at day 1}
    \item{\samp{2}}{a numeric vector: count at day 2}
    \item{\samp{3}}{a numeric vector: count at day 3}
    \item{\samp{4}}{a numeric vector: count at day 4}
    \item{\samp{5}}{a numeric vector: count at day 5}
    \item{\samp{6}}{a numeric vector: count at day 6}
    \item{\samp{7}}{a numeric vector: count at day 7}
    \item{\samp{8}}{a numeric vector: count at day 8}
    \item{\samp{9}}{a numeric vector: count at day 9}
    \item{\samp{10}}{a numeric vector: count at day 10}
    \item{\samp{11}}{a numeric vector: count at day 11}
    \item{\samp{12}}{a numeric vector: count at day 12}
    \item{\samp{13}}{a numeric vector: count at day 13}
    \item{\samp{14}}{a numeric vector: count at day 14}
    \item{\samp{15}}{a numeric vector: count at day 15}
    \item{\samp{16}}{a numeric vector: count at day 16}
    \item{\samp{17}}{a numeric vector: count at day 17}
    \item{\samp{18}}{a numeric vector: count at day 18}
    \item{\samp{19}}{a numeric vector: count at day 19}
    \item{\samp{20}}{a numeric vector: count at day 20}
    \item{\samp{21}}{a numeric vector: count at day 21}
    \item{\samp{22}}{a numeric vector: count at day 22}
    \item{\samp{23}}{a numeric vector: count at day 23}
    \item{\samp{24}}{a numeric vector: count at day 24}
    \item{\samp{27}}{a numeric vector: count at day 27}
    \item{\samp{28}}{a numeric vector: count at day 28}
    \item{\samp{29}}{a numeric vector: count at day 29}
    \item{\samp{31}}{a numeric vector: count at day 31}
    \item{\samp{34}}{a numeric vector: count at day 34}
  }
}
\details{
Every line of data represents a Petri dish. There were 100 seeds per Petri dish. The columns represent the characteristics of each Petri dish. The columns from 4 to 32 represent the number of germinated seeds counted at each assessment time.
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
# Transform to LONG format
data(alfalfaSG)
## maybe str(barley) ; plot(barley) ...
}
\keyword{datasets}
