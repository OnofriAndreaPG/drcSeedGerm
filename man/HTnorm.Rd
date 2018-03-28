\name{HTNorm}
\alias{HTnorm}
\alias{HTnorm.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hydrotime model with normal distribution of base water potential(Bradford, 2002)
}
\description{
This equation corresponds to the hydrotime model for seed germination, which describes the time-course for the proportion of germinated seeds, based on a log-logistic distribution of base water potential within the seed lot.}
\usage{
HT()
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
Bradford, K.J., 2002. Applications of hydrothermal time to quantifying and modeling seed germination and dormancy. Weed Science 50, 248–260.
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
modg <- drm( propCum ~ timeAf + Psi, fct=HTnorm(), data=rape)
summary(modg)
jackGroupSE(modg, rape, rape$Dish)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ seed germination }% use one of  RShowDoc("KEYWORDS")
\keyword{ hydrotime models }% __ONLY ONE__ keyword per line
