\name{plot.spot_curves}
\alias{plot.spot_curves}
\title{ S3 Plot Method }
\description{
S3 plot method for an object of the class \code{"spot_curves"}.
}
\usage{
\method{plot}{spot_curves}(x, multiple = FALSE,
	ylim = c(range(mapply(function(i) range(x[[i]][, 2]),
	seq(x)))) * 100, xlim = c(), type = "l", lty = 1, 
	lwd = 2, expoints = NULL, ylab = "Zero-coupon yields (percent)",
	xlab = "Maturity (years)", main = "Zero-coupon yield curves", ...)
}

\arguments{
  \item{x}{object of the class \code{"spot_curves"}.}
  \item{multiple}{if \code{TRUE} all zero-coupon yield curves are plotted together (default: \code{FALSE}).}
  \item{ylim}{the y limits of the plot, for details see \code{\link{plot.default}}.}
  \item{xlim}{the x limits of the plot, for details see \code{\link{plot.default}}.}
  \item{type}{1-character string giving the type of plot desired, for details see \code{\link{plot.default}}.}
  \item{lty}{the line type, for details see \code{\link{par}}.}
  \item{lwd}{the line width, for details see \code{\link{par} }.}
  \item{expoints}{extrapolation points (default: \code{NULL}.)}
  \item{ylab}{a label for the y axis, for details see \code{\link{plot.default}}.}
  \item{xlab}{a label for the x axis, for details see \code{\link{plot.default}}.}
  \item{main}{a main title for the plot, for details see \code{\link{title}}.}
  \item{\dots}{other graphical parameters, see \code{\link{par}}.}
}
%\details{}
%\value{}
  
%\references{  }
%\author{ Robert Ferstl, Josef Hayden }
%\note{ }
\seealso{\code{\link{plot.df_curves}}, \code{\link{plot.fwr_curves}}, \code{\link{plot.s_curves}}}
%\examples{}

\keyword{  }
