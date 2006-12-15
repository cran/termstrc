\name{loss_function}
\alias{loss_function}

\title{Loss Function used for the Term Structure Estimation }
\description{
  The loss function defines the objective function used for the optimisation. 
  
  
}
\usage{
loss_function(p, phat, omega, weights)
}

\arguments{
  \item{p}{observed prices, yields}
  \item{phat}{estimated prices, yields}
  \item{omega}{weights }
  \item{weights}{if \code{"none"} the squared deviation of \code{p} and
    \code{phat} is not weighted, use \code{"duration"} for a duration weighted optimization}
}
%\details{
% For details please refer to the vignette of the package or to reference.
%}
%\value{}
\references{David Bolder and David Streliski (1999): Yield Curve Modelling at the Bank of Canada.\emph{Technical Report No 84 Bank of Canada}}
\author{Robert Ferstl, Josef Hayden}
%\note{ }
%\seealso{ }
%\examples{}
\keyword{}
