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
\author{
Masanao Yajima \email{yajima@stat.ucla.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link[stats]{image}}
}
\examples{
# XDF
inFile  <- "inData.xdf"
outFile <- "outData.xdf"

# dataset
dataset <- iris[ ,1:4]
form    <- ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width

# conver the dataset to xdf file
rxDataFrameToXdf( dataset, inFile, overwrite = TRUE )

# Fitting Principal Component Analsys Model using stats::princomp
pc.fit <- rsrPrincomp( form, inFile = inFile, outFile=outFile, score = TRUE, overwrite = TRUE )

rsrLoadingPlot( loadings(pc.fit) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }