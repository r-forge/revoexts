\name{factanal.rxCovCor}
\alias{factanal.rxCovCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Factor Analysis for rxCovCor object
}
\description{
A wrapper function around factanal function for rxCovCor object.
}
\usage{
factanal.rxCovCor(x, factors = 2, start = NULL, rotation = "varimax", control = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
rxCovCor object.
}
  \item{factors}{
The number of factors to be fitted.
}
  \item{start}{
NULL or a matrix of starting values, each column giving an initial set of uniquenesses.
}
  \item{rotation}{
character. "none" or the name of a function to be used to rotate the factors: it will be called with first argument the loadings matrix, and should return a list with component loadings giving the rotated loadings, or just the rotated loadings.
}
  \item{control}{
A list of control values, 
\itemize{
    \item{ nstart }{
        The number of starting values to be tried if start = NULL. Default 1. 
    }
    \item { trace } {
        logical. Output tracing information? Default FALSE. 
    }
    \item {lower}{
        The lower bound for uniquenesses during optimization. Should be > 0. Default 0.005. 
    }
    \item{ opt}{
        A list of control values to be passed to optim's control argument. 
    }
    \item{ rotate }{
        A list of additional arguments for the rotation function. }
    }
}
  \item{\dots}{
Other options to be passed to factanal function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{fit }{This is a fitted factanal object.}
\item{scores }{Not used.}}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Masanao Yajima \email{yajima@stat.ucla.edu}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link[stats]{factanal}}
}

\examples{
# Covariance matrix calculation
# define formula
form   <- ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width
# calculate covariance matrix
allCov <- rxCovCor( form, data = iris, type = "Cov" )

fit <- factanal.rxCovCor( allCov, factors = 1 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
