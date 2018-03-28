\name{HTW2}
\alias{HTW2}
\alias{HTW2.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hydrotime model, with Weibull type II distribution of base water potential
}
\description{
This equation corresponds to the hydrotime model for seed germination, which describes the time-course for the proportion of germinated seeds, based on a Weibull type II distribution of base water potential within the seed lot.
}
\usage{
HTW2()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
This function has no arguments}

\details{
g =  f(Psi)}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Watt, 2006
}
\author{
Andrea Onofri
}
\note{
This function is for use with the function drm
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(rape)
modg <- drm( propCum ~ timeAf + Psi, fct=HTW2(), data=rape)
summary(modg)
jackGroupSE(modg, rape, rape$Dish)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
