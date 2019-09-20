\name{reglin}
\alias{reglin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation in Simple Linear Regression
}
\description{
Estimation of all parameters in the the simple linear regression model.
}
\usage{
reglin(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
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
function (x, y) 
{
    n <- length(x)
    xbar <- mean(x)
    ybar <- mean(y)
    x2bar <- mean(x^2)
    y2bar <- mean(y^2)
    xybar <- mean(x * y)
    beta1hat <- (xybar - xbar * ybar)/(x2bar - xbar^2)
    beta0hat <- ybar - beta1hat * xbar
    sigmahat <- sqrt(n * (y2bar + beta0hat^2 + beta1hat^2 * x2bar - 
        2 * beta0hat * ybar - 2 * beta1hat * xybar + 2 * beta0hat * 
        beta1hat * xbar)/(n - 2))
    sigmahat.beta1hat <- sigmahat/sqrt(n * (x2bar - xbar^2))
    test.statistic <- beta1hat/sigmahat.beta1hat
    return(list(beta0hat = beta0hat, beta1hat = beta1hat, sigmahat = sigmahat, 
        sigmahat.beta1hat = sigmahat.beta1hat, test.statistic = test.statistic))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
