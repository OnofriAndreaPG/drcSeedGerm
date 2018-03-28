\name{PmaxPsi1}
\alias{PmaxPsi1}
\alias{PmaxPsi1.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shifted exponential distribution of base osmotic potential
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PmaxPsi1()
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function ()
{
    Pmax1.fct <- function(x, parm) {
        Pmax <- PmaxPsi1.fun(x, parm[, 3], parm[, 1], parm[,
            2])
        return(Pmax)
    }
    Pmax1.names <- c("Psib", "sigma", "G")
    Pmax1.text <- "Shifted exponential distribution of base osmotic potential"
    Pmax1.ss <- function(data) {
        data <- subset(data, data[, 2] != 0)
        x <- data[, 1]
        y <- data[, 2]
        G <- max(y) * 1.05
        pseudoY <- -log((G - y)/G)
        pseudoX <- x
        coefs <- coef(lm(pseudoY ~ pseudoX))
        a <- coefs[1]
        b <- coefs[2]
        sigma <- 1/b
        Psib <- -a * sigma
        return(c(Psib, sigma, G))
    }
    Pmax1 <- list(fct = Pmax1.fct, ssfct = Pmax1.ss, names = Pmax1.names,
        text = Pmax1.text)
    class(Pmax1) <- "drcMean"
    invisible(Pmax1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
