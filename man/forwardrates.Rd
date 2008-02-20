\name{forwardrates}
\alias{forwardrates}

\title{Forward Rate Calculation }
\description{Calculates forward rates according to the Nelson/Siegel or Svensson approach.}
\usage{
forwardrates(method, beta, m)
}

\arguments{
  \item{method}{method used for the forward rate calculation (\code{"Nelson/Siegel"} or \code{"Svensson"}).}
  \item{beta}{parameters corresponding to the chosen method.}
  \item{m}{one maturity or a vector of maturities.}
}
%\details{}
\value{The function returns a vector with the calculated forward rates.}
  
%\references{  }

%\note{ }
\seealso{\code{\link{fwr_ns}}, \code{\link{fwr_sv}}}
%\examples{}

\keyword{  }