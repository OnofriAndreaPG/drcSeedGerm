\name{HTG}
\alias{HTG}
\alias{HTG.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hydrotime model with logistic distribution of base water potential (Mesgaran et al., 2013)
}
\description{
This equation corresponds to the hydrotime model for seed germination, which describes the time-course for the proportion of germinated seeds, based on a logistic distribution of base water potential within the seed lot.}
\usage{
HTG()
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
Mesgaran, M.B., Mashhadi, H.R., Alizadeh, H., Hunt, J., Young, K.R., Cousens, R.D., 2013. Importance of distribution function selection for hydrothermal time models of seed germination. Weed Research 53, 89–101.
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
modg <- drm( propCum ~ timeAf + Psi, fct=HTG(), data=rape)
summary(modg)
jackGroupSE(modg, rape, rape$Dish)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ seed germination }% use one of  RShowDoc("KEYWORDS")
\keyword{ hydrotime models }% __ONLY ONE__ keyword per line
