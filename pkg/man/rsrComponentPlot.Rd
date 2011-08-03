\name{rsrComponentPlot}
\alias{rsrComponentPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Principal Component Plot
}
\description{
Plots the principal components of the fitted princomp object
}
\usage{
rsrComponentPlot(x, componentName = "Comp.1", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  rsrPrincomp object
}
  \item{componentName}{
Name of a component, it should be "Comp." followed by a number.
}
  \item{\dots}{
  Other options to be passed to rxLinePlot
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
    \code{\link[RevoExts]{rsrFactorPlot}}
}
\examples{
# conver the dataset to xdf file
# rxDataFrameToXdf( iris, "inData.xdf", overwrite = TRUE )

# XDF
dir     <- system.file( "data", package = "RevoExts" )
inFile  <- file.path( dir, "inData.xdf" )
outFile <- file.path( dir, "outData.xdf" )

# dataset
form    <- ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width

# Fitting Principal Component Analsys Model using stats::princomp
pc.fit <- rsrPrincomp( form, inFile = inFile, outFile=outFile, score = TRUE, overwrite = TRUE )

rsrComponentPlot( pc.fit, componentName = "Comp.1", type = "p" )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }