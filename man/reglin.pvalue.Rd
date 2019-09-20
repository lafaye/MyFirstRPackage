\name{reglin.pvalue}
\alias{reglin.pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
P-value for the Simple Linear Regression Model
}
\description{
Computation of p-values for the simple linear regression model using a randomisation approach.
}
\usage{
reglin.pvalue(x, y, M = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
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
function (x, y, M = 10000) 
{
    tobs <- reglin(x, y)$test.statistic
    tM <- rep(NA, M)
    for (m in 1:M) {
        tM[m] <- reglin(sample(x), y)$test.statistic
    }
    pvalue <- mean(abs(tM) > abs(tobs))
    return(pvalue)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
