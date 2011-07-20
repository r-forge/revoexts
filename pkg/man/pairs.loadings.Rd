\name{pairs.loadings}
\alias{pairs.loadings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix scatter plot for loading matrix
}
\description{
Produces bivariate scatter plot matrix for loading matrix with spceified columns.
}
\usage{
pairs.loadings(loading, include = 1:3, col = "black", labels = NULL, compnames = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loading}{
%%     ~~Describe \code{loading} here~~
}
  \item{include}{
Columns of loading matrix to include in the plot
}
  \item{col}{
color of each labels.  If only one color is specified it will be used for all labels.  Otherwise the color must have length of rows of the loading matrix.
}
  \item{labels}{
Label to be used instead of row names of the loading matrix.
}
  \item{compnames}{
Names to be used for the diagonal of the plot matrix instead of the column names of the loading matrix.
}
  \item{\dots}{
Other options to be passed to the plot function.
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
%}\references{
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
    \code{\link[stats]{pairs}}
}
\examples{
pc.fit <- princomp( iris[,1:4] )
pairs( pc.fit$loadings, include=1:3 ) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }