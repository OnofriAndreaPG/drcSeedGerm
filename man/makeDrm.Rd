\name{makeDrm}
\alias{makeDrm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reshape a seed germination datasets for time-to-event analyses
}
\description{
MakeDrm: this function reshapes a common field book for germination assays into the form that is required for time-to-event analyses with drm() in the drc package and drmte() in the drcte package. The common field book has one row per each Petri dish and the counts of germinated seeds at each assessment time are listed in different columns.
}
\usage{
makeDrm(counts, treat, nViable, moniTimes)
}

\arguments{
\item{counts}{
Dataframe listing the counts of germinated seeds in each Petri dish (rows) at each assessment time (columns).
}
\item{treat}{
Dataframe listing for each dish the levels of each treatment. It has as many columns as there are treatments in the experiment.
}
\item{nViable}{
A vector with the number of viable seeds per dish, at the beginning of the assay
}
\item{moniTimes}{
Vector of monitoring times. It needs to be of the same length as the number of columns in 'counts'.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a dataframe
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andrea Onofri
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# makeDrm
library(drcSeedGerm)
data(lotusOr)
counts <- lotusOr[,3:length(lotusOr[1,])]
treat <- data.frame(tratt=lotusOr[,1])
nViable <- rep(25,12)
moniTimes <- c(1:15)
datasetDrc <- makeDrm(counts=counts, treat=treat, nViable=nViable, moniTimes)
head(datasetDrc, 30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
