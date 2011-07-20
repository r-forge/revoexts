\name{rsrLoadingPlot}
\alias{rsrLoadingPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image plot of loading matrix
}
\description{
Wrapper on image plot for loading matrix.
}
\usage{
rsrLoadingPlot(loading, main = "Loading matrix", ylabels, xlabels, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loading}{
Loading matrix
}
  \item{main}{
Title of the plot
}
  \item{ylabels}{
Labels for the rows of the loading matrix.
}
  \item{xlabels}{
Labels for the columns of the loading matrix.
}
  \item{\dots}{
Other options to be passed to image function}
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
Masanao Yajima \email{yajima@stat.ucla.edu}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link[stats]{image}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }