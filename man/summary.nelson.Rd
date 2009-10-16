\name{summary.nelson}
\alias{summary.nelson}

\title{ S3 Summary Method }
\description{
  S3 summary method for objects of the class \command{"nelson"}.
}
\usage{
\method{summary}{nelson}(object,...)
}

\arguments{
  \item{object}{object of the class \command{"nelson"}.}
  \item{\dots}{other arguments.}
}
\details{
  The summary method for an object of the class \command{"nelson"} prints the solution of
 the goodness of fit statistics (RMSE,AABSE) of the optimisation. Moreover a convergence information of the used
 optimiser (\command{nlminb}) is printed. 
}
%\value{}
%\references{ }
%\author{Robert Ferstl, Josef Hayden }
%\note{}
\seealso{ \command{\link{nlminb}}, \command{\link{plot.nelson}}, \command{\link{print.nelson}}, \command{\link{rmse}}, \command{\link{aabse}}}
%\examples{}


