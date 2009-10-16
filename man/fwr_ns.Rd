\name{fwr_ns}
\alias{fwr_ns}

\title{Forward Rate Calculation according to Nelson/Siegel }
\description{Calculate forward rates according to the Nelson/Siegel (1987) approach.}
\usage{
fwr_ns(beta, m)
}

\arguments{
  \item{beta}{parameter vector  \eqn{{\bm{\beta}} = \left(\beta_0,\beta_1,\beta_2,\tau_1\right)}. }
  \item{m}{maturity or maturity vector.}
}
\details{The forward rate for a maturity \eqn{m} is calculated using the following relationsship:

\deqn{f(m,\bm{\beta}) = \beta_0+\beta_1\exp\left(-\frac{m}{\tau_1}\right)+\beta_2\left[\left(\frac{m}{\tau_1}\right)\exp\left(-\frac{m}{\tau_1}\right)\right].}
}
\value{The function returns the calculated forward rate (vector).}


\references{  
Charles R. Nelson and Andrew F. Siegel (1987): 
Parsimonious Modeling of Yield Curves.
\emph{The Journal of Business}, \bold{60(4)}:473--489.}

%\author{ Robert Ferstl, Josef Hayden }
%\note{ }
\seealso{\code{\link{fwr_sv}}, \code{\link{forwardrates}}}
%\examples{}


