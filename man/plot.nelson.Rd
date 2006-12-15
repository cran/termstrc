\name{plot.nelson}
\alias{plot.nelson}

\title{S3 Plot Method }
\description{
 S3 plot method for an object of the class \code{nelson}
}
\usage{
\method{plot}{nelson}(x,matrange,pdf=FALSE,...)
}
\arguments{
  \item{x}{ object of the class \code{nelson}}
  \item{matrange}{maturity range for the plot, e.g. c(2,10)}
  \item{pdf}{ If \code{TRUE} the graphics are plotted in a pdf file,
    which is saved in the current working directory}
  \item{\dots}{other arguments}
}

\details{
A plot of the estimated zero-coupon yield curve including the calculated yields for each specified group of bonds is offered. 
In addition all zero-coupon yield curves are plotted in one graphic and afterwards the spread curves are plotted. 
The occurrence of a partially dashed line in a graphic indicates that the the curve has been extrapolated. 
}
 
%\value{ }
%\references{ }
\author{Robert Ferstl, Josef Hayden}

%\note{ }

\seealso{ \code{\link{print.nelson}}, \code{\link{summary.nelson}}}
%\examples{}
\keyword{  }

