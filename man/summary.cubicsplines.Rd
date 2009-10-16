\name{summary.cubicsplines}
\alias{summary.cubicsplines}

\title{S3 Summary Method for Cubicsplines }
\description{
  S3 summary method for objects of the class \command{"cubicsplines"}.
}
\usage{
\method{summary}{cubicsplines}(object,...)
}

\arguments{
  \item{object}{object of the class \command{"cubicsplines"}.}
  \item{\dots}{other arguments.}
}
\details{
 The summary method for an object of the class \command{"cubicsplines"} prints goodness of fit statistics (RMSE, AABSE) of the optimisation. Additionally summary statistics of the regression analysis of the parameters are printed. 
}
%\value{}
%\references{ }
%\author{Robert Ferstl, Josef Hayden }
%\note{}
\seealso{\command{\link{plot.cubicsplines}}, \command{\link{print.cubicsplines}}, \command{\link{rmse}}, \command{\link{aabse}}, \command{\link{summary.lm}}}
%\examples{}

