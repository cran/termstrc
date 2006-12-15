\name{rmse}
\alias{rmse}
\title{Root Mean Squared Error}
\description{Root mean squared error}

\usage{
rmse(actual, estimated)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{actual}{vector, consisting of the observed values}
  \item{estimated}{vector, consisting of the estimated values}
}
\details{
  Calculation of the RMSE according to the formula:
\deqn{RMSE=\sqrt{\frac{1}{n}\sum_{i=1}^{n}{(e_i-\overline{e_i})^2}}}

}
\references{David Bolder and David Streliski (1999): Yield Curve Modelling at the Bank of Canada.\emph{Technical Report No 84 Bank of Canada}.}

\author{Robert Ferstl, Josef Hayden}

\seealso{\code{\link{aabse}}}

\keyword{}
