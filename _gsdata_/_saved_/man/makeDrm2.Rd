\name{makeDrm2}
\alias{makeDrm2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reshape a seed germination datasets to obey the requirements of drm()
}
\description{
MakeDrm2: this function reshapes a datasets organised as necessary for nonlinear regression into the kind of dataset required by the function drm() in the drc package. It works with either the counts of germinated seeds at each monitoring time or the cumulative counts at each monitoring time.
}
\usage{
makeDrm2(counts, treat, nViable, moniTimes, Dish,  cumulative=T)
}

\arguments{
\item{counts}{
Vector listing the counts of germinated seeds in each Petri dish at each assessment time.
}
\item{treat}{
Dataframe listing, for each row of data, the corresponding level of experimental factors (one factor per column)
}
\item{nViable}{
A vector listing the number of viable seeds, at the beginning of the assay. This number is the same for all observations belonging to the same dish.
}
\item{moniTimes}{
Vector of monitoring times.
}
\item{Dish}{
Vector of codes for dishes.
}
\item{cumulative}{
Logical: True if counts are cumulative, False if they are not.
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
#makeDrm2
data(lotusCum)
moniTime <- lotusCum$Time
count <- lotusCum$nCum
nViable <- rep(25, length(lotusCum[,1]))
Dish <- as.factor(lotusCum$Dish)
treatGroups <- lotusCum[,1]
dataset_sd <- makeDrm2(count, treatGroups, nViable, moniTime, Dish)
head(dataset_sd)
count <- lotusCum$nSeeds
dataset_sd <- makeDrm2(count, treatGroups, nViable, moniTime, Dish,  cumulative = FALSE)
head(dataset_sd)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
