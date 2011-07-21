\name{biplot.rsrPrincomp}
\alias{biplot.rsrPrincomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binned biplot for rsrPrincomp object
}
\description{
biplot for rsrPrincomp object.  Since the data may be massive, the data are binned.
}
\usage{
biplot.rsrPrincomp(x, choices = 1L:2L, label = NULL, label.nm = NULL, numBreaksX = 100, numBreaksY = 100, expand = 1, scale = 1, pc.biplot = FALSE, arrow.len = 0.1, alpha = 1, main = NULL, xlab, ylab, xlim, ylim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  rsrPrincomp object.
}
  \item{choices}{
choices of components to plot.
}
  \item{label}{
XDF file where the labels are stored.
}
  \item{label.nm}{
Name of the variable for the label. 
}  \item{numBreaksX}{
Bin breaks for x.
}
  \item{numBreaksY}{
Bin breaks for y.
}
  \item{expand}{
An expansion factor to apply when plotting the second set of points relative to the first. This can be used to tweak the scaling of the two sets to a physically comparable scale.
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{pc.biplot}{
%%     ~~Describe \code{pc.biplot} here~~
}
  \item{arrow.len}{
The length of the arrow heads on the axes plotted in var.axes is true. The arrow head can be suppressed by arrow.len = 0.
}
\item{alpha}{
Alpha level for the points in the plot
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
    \code{\link[stats]{biplot}}
}
\examples{
# XDF
inFile  <- "inData.xdf"
outFile <- "outData.xdf"

# dataset
dataset <- iris[ ,1:4]
form    <- ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width

# conver the dataset to xdf file
rxDataFrameToXdf( iris, inFile, overwrite = TRUE )

# Fitting Principal Component Analsys Model using stats::princomp
pc.fit <- rsrPrincomp( form, inFile = inFile, outFile=outFile, score = TRUE, overwrite = TRUE )

biplot( pc.fit, label=inFile, label.nm="Species" )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }