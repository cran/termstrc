\name{aabse}
\alias{aabse}

\title{Average Absolute Mean Error}

\description{Average absolute mean error}
\usage{
aabse(actual, estimated)
}

\arguments{
  \item{actual}{vector, consisting of the observed values}
  \item{estimated}{vector, consisting of the estimated values}
}

}
\details{
  Calculation of the AABSE according to the formula:
\deqn{AABSE=\frac{1}{n}\sum_{i=1}^{n}{|e_i-\overline{e_i}|}}
}
%\value{}
\references{David Bolder and David Streliski (1999): Yield Curve Modelling at the Bank of Canada.\emph{Technical Report No 84 Bank of Canada}}

\author{Robert Ferstl, Josef Hayden}
\seealso{\code{\link{rmse}}}

\keyword{}
