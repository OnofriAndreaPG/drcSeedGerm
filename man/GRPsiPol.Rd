\name{GRPsiPol}
\alias{GRPsiPol}
\alias{GRPsiPol.fun}
\alias{GRPsiPol2}
\alias{GRPsiPol2.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polynomial hydrotime model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GRPsiPol()
}
%- maybe also 'usage' for other objects documented here.
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(excalibur)
head(excalibur)
library(drc)

modGR1 <- drm(GR ~ Psi, fct=GRPsiLin(), data=excalibur, curveid=Perc)
summary(modGR1)
plot(modGR1, log="", legendPos=c(-0.9, 1))

modGR2 <- drm(GR ~ Psi, fct=GRPsiPol2(), data=excalibur, curveid=Perc)
summary(modGR2)
plot(modGR2, log="", legendPos=c(-0.9, 1))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
