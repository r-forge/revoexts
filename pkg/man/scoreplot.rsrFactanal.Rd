\name{scoreplot.rsrFactanal}
\alias{scoreplot.rsrFactanal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binned plot for rsrFactanal scores
}
\description{
Binned plot for rsrFactanal scores.  Since the data may be massive, the data are binned.
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
XDF file where the labels are stored.
}
  \item{label.nm}{
Name of the variable for the label. 
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
# Example
v1 <- c(1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,4,5,6) 
v2 <- c(1,2,1,1,1,1,2,1,2,1,3,4,3,3,3,4,6,5) 
v3 <- c(3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,5,4,6) 
v4 <- c(3,3,4,3,3,1,1,2,1,1,1,1,2,1,1,5,6,4) 
v5 <- c(1,1,1,1,1,3,3,3,3,3,1,1,1,1,1,6,4,5) 
v6 <- c(1,1,1,2,1,3,3,3,4,3,1,1,1,2,1,6,5,4) 
label <-c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4)

dataset <- data.frame( cbind( v1, v2, v3, v4, v5, v6 ) )
form   <- ~ v1 + v2 + v3 + v4 +v5 +v6

# conver the dataset to xdf file
rxDataFrameToXdf( dataset, "Data.xdf", overwrite = TRUE )

# XDF
dir     <- system.file( "data", package = "RevoExts" )
inFile  <- file.path( dir, "Data.xdf" )
outFile <- file.path( dir, "outData.xdf" )

# Principal Component Analsys ( stats::princomp )
fa.fit  <- rsrFactanal( form, inFile = inFile, outFile = outFile, 
                        factors = 2, overwrite = TRUE, scores = "regression" )

scoreplot( fa.fit, choices = 1:2, label = inFile, label.nm = "label", numBreaksX = 10, numBreaksY = 10,  )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }