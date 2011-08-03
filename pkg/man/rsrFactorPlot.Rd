\name{rsrFactorPlot}
\alias{rsrFactorPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot for factor scores
}
\description{
Plots factor scores using rxLinePlot.
}
\usage{
rsrFactorPlot(x, factorName = "Factor1", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
rsrFactanal Object
}
  \item{factorName}{
Name of a factor.  It should be "Factor" followed by a number with no space.
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link[RevoExts]{rsrComponentPlot}}
}
\examples{
# Example
v1 <- c(1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,4,5,6) 
v2 <- c(1,2,1,1,1,1,2,1,2,1,3,4,3,3,3,4,6,5) 
v3 <- c(3,3,3,3,3,1,1,1,1,1,1,1,1,1,1,5,4,6) 
v4 <- c(3,3,4,3,3,1,1,2,1,1,1,1,2,1,1,5,6,4) 
v5 <- c(1,1,1,1,1,3,3,3,3,3,1,1,1,1,1,6,4,5) 
v6 <- c(1,1,1,2,1,3,3,3,4,3,1,1,1,2,1,6,5,4) 

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
rsrFactorPlot( fa.fit, factorName = "Factor1", type = "p" )
                        
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graph }