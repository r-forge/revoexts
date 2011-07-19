\name{rsrFactanal}
\alias{rsrFactanal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Factor analysis for xdf object
}
\description{
A wrapper function around factanal function for xdf file object.
}
\usage{
rsrFactanal(formula, factors, inFile, outFile = NULL, append = "none", overwrite = FALSE, scores = c("none", "regression", "Bartlett"), start = NULL, rotation = "varimax", blocksPerRead = rxGetOption("blocksPerRead"), reportProgress = rxGetOption("reportProgress"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Formula object
}
  \item{factors}{
 The number of factors to be fitted. 
}
  \item{inFile}{
XDF file where the data is stored.
}
  \item{outFile}{
XDF file where calculated factor score will be stored.
}
  \item{append}{
Not used.
}
  \item{overwrite}{
  Logical with default as FALSE.  If you are using preexisting XDF file this option need to be set to TRUE.
}
  \item{scores}{
  Type of scores to produce, if any. The default is none, "regression" gives Thompson's scores, "Bartlett" given Bartlett's weighted least-squares scores. Partial matching allows these names to be abbreviated.  See help on factanal for detail.
}
  \item{start}{
NULL or a matrix of starting values, each column giving an initial set of uniquenesses.
}
  \item{rotation}{
character. "none" or the name of a function to be used to rotate the factors: it will be called with first argument the loadings matrix, and should return a list with component loadings giving the rotated loadings, or just the rotated loadings.
}
  \item{blocksPerRead}{
number of blocks to read for each chunk of data read from the data source.
}
  \item{reportProgress}{
integer value with options: 
\begin{itemize}
	\item 0: no progress is reported. 
	\item 1: the number of processed rows is printed and updated. 
	\item 2: rows processed and timings are reported. 
	\item 3: rows processed and all timings are reported. 
\end{itemize}
}
  \item{\dots}{
Other options to be passed to factanal
}
}
\details{
This function is a wrapper around factanal function.  See help on factanal for detail.
}
\value{
If it is a LIST, use
\item{fit }{This is a fitted factanal object.}
\item{scores }{Name of the XDF file where the factor score is stored.}
}
\references{

}
\author{
Masanao Yajima \email{yajima@stat.ucla.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link[stats]{factanal}}
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

# XDF file names 
inFile  <- "inData.xdf"
outFile <- "outData.xdf"

# dataset
rxDataFrameToXdf( dataset,  inFile, overwrite = TRUE )

# Principal Component Analsys ( stats::princomp )
fa.fit  <- rsrFactanal( form, inFile = inFile, outFile = outFile, 
                        factors = 2, overwrite = TRUE, scores = "regression" )
}
\keyword{ manip }
