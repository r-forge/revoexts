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
rsrLoadingPlot(x, main = "Loading matrix", ylabels, xlabels, mar = c(0, 5, 4, 0), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
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
  \item{mar}{
Size of margin see help on par for detail.
}

  \item{\dots}{
Other options to be passed to image function}
}
\author{
Masanao Yajima \email{yajima@stat.ucla.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link[graphics]{image}}
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

library( RColorBrewer )
rsrLoadingPlot( loadings(pc.fit) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }