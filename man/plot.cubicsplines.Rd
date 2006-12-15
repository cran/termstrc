\name{plot.cubicsplines}
\alias{plot.cubicsplines}

\title{S3 Plot Method for Cubic Splines }
\description{
 S3 plot method for an object of the class \code{cubicsplines}
}
\usage{
\method{plot}{cubicsplines}(x,matrange,pdf=FALSE,...)
}
\arguments{
  \item{x}{object of the class \code{cubicsplines}}
  \item{matrange}{maturity range for the plot, e.g. c(2,10)}
  \item{pdf}{if \code{TRUE} the graphics are plotted in a pdf file, which is saved in the current working directory}	
  \item{\dots}{other arguments}
}

\details{
A plot of the estimated zero-coupon yield curve including the calculated yields for each specified group of bonds is offered. 
The occurrence of dashed vertical lines in a graphic represents the knot
 points used for the cubic splines estimation. 
In addition all zero-coupon yield curves are plotted in one graphic and afterwards the spread curves are plotted. 
}
 
%\value{ }
%\references{ }
\author{Robert Ferstl, Josef Hayden}

%\note{ }

\seealso{ \code{\link{print.cubicsplines}}, \code{\link{summary.cubicsplines}}}
%\examples{}
\keyword{  }

