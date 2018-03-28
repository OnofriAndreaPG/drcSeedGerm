\name{HTEx}
\alias{HTEx}
\alias{HTEx.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hydrotime model with Type II Extreme distribution of base water potential
}
\description{
This equation corresponds to the hydrotime model for seed germination, which describes the time-course for the proportion of germinated seeds, based on a type II extreme distribution of base water potential within the seed lot.}
\usage{
HTEx()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  This function has no arguments
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
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
modg <- drm( propCum ~ timeAf + Psi, fct=HTEx(), data=rape)
summary(modg)
jackGroupSE(modg, rape, rape$Dish)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ seed germination }% use one of  RShowDoc("KEYWORDS")
\keyword{ hydrotime models }% __ONLY ONE__ keyword per line
