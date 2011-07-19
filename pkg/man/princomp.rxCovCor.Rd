\name{princomp.rxCovCor}
\alias{princomp.rxCovCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Principal Component Analysis for rxCovCor object
}
\description{
A wrapper function around princomp function for rxCovCor object.

}
\usage{
princomp.rxCovCor(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
rxCovCor object.
}
  \item{\dots}{
  Other options to be passed to princomp function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{fit }{This is a fitted princomp object.}
\item{scores }{Not used.}
}
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
    \code{\link[stats]{princomp}}
}
\examples{
# Covariance matrix calculation
# define formula
form   <- ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width
# calculate covariance matrix
allCov <- rxCovCor( form, data = iris, type = "Cov" )

# using princomp directly from stats package
fit.1 <- princomp( covmat = allCov$CovCor )
# using princomp for rxCovCor object
fit.2 <- princomp( allCov )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
