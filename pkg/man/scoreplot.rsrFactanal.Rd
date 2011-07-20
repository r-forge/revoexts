\name{scoreplot.rsrFactanal}
\alias{scoreplot.rsrFactanal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoreplot.rsrFactanal(x, choices = 1L:2L, label = NULL, label.nm, alpha = 1, 
                      numBreaksX = 100, numBreaksY = 100, main = NULL, xlab, ylab, xlim, ylim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
rsrFactanal Object
}
  \item{choices}{
Choices of two factors to plot.
}
  \item{label}{
Labels used for plots.
}
  \item{label.nm}{
%%     ~~Describe \code{label.nm} here~~
}
  \item{alpha}{
Alpha level for the points in the plot
}
  \item{numBreaksX}{
Bin breaks for x.
}
  \item{numBreaksY}{
Bin breaks for y.
}
  \item{main}{
Main title.
}
  \item{xlab}{
Title of x axis.
}
  \item{ylab}{
Titile of y axis
}
  \item{xlim}{
Limit of x axis.
}
  \item{ylim}{
Limit of y axis.
}
  \item{\dots}{
Other options to be passed to rsrGridMap option.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Masanao Yajima \email{yajima@stat.ucla.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link[stats]{biplot}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }