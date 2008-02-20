\name{plot.cubicsplines}
\alias{plot.cubicsplines}

\title{S3 Plot Method for Cubic Splines }
\description{
 S3 plot method for an object of the class \code{"cubicsplines"}.
}

\usage{
\method{plot}{cubicsplines}(x, matrange = c(min(mapply(function(i) min(x$y[[i]][, 1]), 
	seq(x$n_group))), max(mapply(function(i) max(x$y[[i]][, 1]), 
	seq(x$n_group)))), multiple = FALSE, expoints= NULL,
	ctype = "spot", lwd=2, lty=1, type = "l", 
	errors = "price", inset = c(0.8, 0.1), ask=TRUE, ...)
}

\arguments{
 \item{x}{object of the class \code{"cubicsplines"}.}
  \item{matrange}{ maturity range for the plot, e.g. c(2,10).}
  \item{multiple}{if \code{TRUE} all curves are plotted together (default: \code{FALSE}).}
  \item{expoints}{extrapolation points (default: \code{NULL}).}
  \item{ctype}{ parameter setting for the desired curve type, \code{"spot"} (\code{"forward"}, \code{"discount"}, \code{"spread"}) for the spot rate (forward rate, discount factor, spread) curves. Use \code{"none"} if no curve plot is desired.  }
 \item{errors}{Specify the type of the error plot. If \code{"price"} (\code{"yield"}) the price (yield) errors will be plot. Use \code{"none"} if no error plot is desired. }
\item{lwd}{the line width, for details see \code{\link{par}}.}
\item{lty}{the line type, for details see \code{\link{par}}.}
  \item{type}{1-character string giving the type of plot desired, for details see \code{\link{plot.default}}.}
  \item{inset}{ inset distance(s) from the margins as a fraction of the plot region, for details see \code{\link{legend}}.}
  \item{ask}{if \code{TRUE} (and the \R session is interactive) the user is asked for input, before a new figure is drawn, see \code{\link{par}} for details.}
  \item{\dots}{other graphical parameters, see \code{\link{par}}. }
}
\details{
Depending on the choice of the curve type (\code{"spot"}, \code{"forward"}, \code{"discount"}, \code{"spread"}) the corresponding
curves will be plot. Either separately or together (\code{mulitple = TRUE}). If the curves are plotted separately also the knot points
used for the estimation of the cubic splines and the yield-to-maturities will be plotted. In addition, with a zero-coupon yield curve plot the 95 \% confidence interval of the curve will be plot. To ease the analysis of the goodness of the estimation, serveral error plots for the yield and price error are offered. }  

%\value{}
%\references{  }
%\author{Robert Ferstl, Josef Hayden }
%\note{}
\seealso{ \code{\link{plot.df_curves}},  \code{\link{plot.error}},  \code{\link{plot.fwr_curves}},  \code{\link{plot.ir_curve}},  \code{\link{plot.s_curves}},  \code{\link{plot.spot_curves}},  \code{\link{plot.cubicsplines}} }
%\examples{}

\keyword{  }

