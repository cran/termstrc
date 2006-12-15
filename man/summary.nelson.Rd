\name{summary.nelson}
\alias{summary.nelson}

\title{ S3 Summary Method }
\description{
  S3 summary method for objects of the class \code{nelson}
}
\usage{
\method{summary}{nelson}(object,...)
}

\arguments{
  \item{object}{object of the class \code{nelson}}
  \item{\dots}{other arguments}
}
\details{
  The summary method for an object of the class \code{nelson} prints the solution of
 the goodness of fit tests (RMSE,AABSE) of the optimisation. Moreover a convergence information of the used
 optimiser (nlminb) is printed. 
}
%\value{}
%\references{ }
\author{Robert Ferstl, Josef Hayden }
%\note{}
\seealso{ \code{\link{nlminb}}, \code{\link{plot.nelson}}, \code{\link{print.nelson}}, \code{\link{rmse}}, \code{\link{aabse}}}
\examples{}

\keyword{ }

