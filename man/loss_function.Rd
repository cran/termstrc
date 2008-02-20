\name{loss_function}
\alias{loss_function}

\title{Loss Function used for the Term Structure Estimation }
\description{
  The loss function defines the objective function used for the optimisation. In case of term structure estimation
 the objective function is either the minimsation of the (weighted) squared price or the yield errors.
}
\usage{
loss_function(p, phat, omega, weights)
}

\arguments{
  \item{p}{vector of observed prices (yields).}
  \item{phat}{vector of estimated prices (yields).}
  \item{omega}{weights vector. }
  \item{weights}{if \code{"none"} the squared deviation of \code{p} and
    \code{phat} is not weighted, use \code{"duration"} for a duration weighted optimisation.}
}
%\details{
% 
%}
%\value{}
%\author{Robert Ferstl, Josef Hayden}
%\note{ }
%\seealso{ }
%\examples{}
\keyword{}
