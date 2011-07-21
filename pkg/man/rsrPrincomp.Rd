\name{rsrPrincomp}
\alias{rsrPrincomp}
\alias{print.rsrPrincomp}
\alias{summary.rsrPrincomp}
\alias{plot.rsrPrincomp}
\alias{scores.rsrPrincomp}
\alias{loadings.rsrPrincomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Principal Component Analysis for xdf object
}
\description{
A wrapper function around princomp function for xdf file object.
}
\usage{
rsrPrincomp(formula, inFile, score = FALSE, outFile = NULL, append = "none", overwrite = FALSE, blocksPerRead = rxGetOption("blocksPerRead"), reportProgress = rxGetOption("reportProgress"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
 formula, as described in rxFormula, with all the terms on the right-hand side of the ~ separated by + operators. Each term may be a single variable, a transformed variable, or the interaction of (transformed) variables separated by the : operator. e.g. ~ x1 + log(x2) + x3 : x4 
}
  \item{inFile}{
XDF file where the data is stored.
}
  \item{score}{
Logical if score should be calculated or not, default is FALSE.
}
  \item{outFile}{
XDF file where calculated factor score will be stored.
}  \item{append}{
Not used.
}
  \item{overwrite}{
  Logical with default as FALSE.  If you are using preexisting XDF file this option need to be set to TRUE.
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

}  \item{\dots}{
Other options to be passed to princomp function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{fit }{This is a fitted princomp object.}
\item{scores }{Name of the XDF file where the component score is stored.}
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
    \code{\link[stats]{princomp}}
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
