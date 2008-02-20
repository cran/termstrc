\name{plot.s_curves}
\alias{plot.s_curves}

\title{S3 Plot Method}
\description{
 S3 plot method for an object of the class \code{"s_curves"}. 
}
\usage{
\method{plot}{s_curves}(x, xlim = c(range(mapply(function(i) range(x[[i]][, 1]), seq(x)))),
	ylim = c(range(mapply(function(i) range(x[[i]][, 2]), 
	seq(x)))) * 10000, expoints = NULL, 
	xlab = "Maturity (years)", ylab = "Spread (basis points)",
	lwd = 2, lty = 1, main = "Spread curves", ...)
}

\arguments{
  \item{x}{object of the class \code{"s_curves"}.}
  \item{ylim}{the y limits of the plot, for details see \code{\link{plot.default}}.}
  \item{xlim}{the x limits of the plot, for details see \code{\link{plot.default}}.}
  \item{lty}{the line type, for details see \code{\link{par}}.}
  \item{lwd}{the line width, for details see \code{\link{par}}.}
  \item{expoints}{extrapolation points (default: \code{NULL}).}
  \item{ylab}{a label for the y axis, for details see \code{\link{plot.default}}.}
  \item{xlab}{a label for the x axis, for details see \code{\link{plot.default}}.}
  \item{main}{a main title for the plot, for details see \code{\link{title}}.}
  \item{\dots}{other graphical parameters, see \code{\link{par}}.}
}
\details{
The spread curves (the difference of zero-coupon yield curves) are plotted, if at least two groups of bonds
were specified.
}
%\value{}
%\references{  }
%\author{ Robert Ferstl, Josef Hayden }
%\note{ }
\seealso{\code{\link{plot.df_curves}}, \code{\link{plot.fwr_curves}}, \code{\link{plot.spot_curves}}}
%\examples{}

\keyword{  }
