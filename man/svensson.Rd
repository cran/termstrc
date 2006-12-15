\name{svensson}
\alias{svensson}

\title{Spot Rate Function according to Svensson }
\description{This function calculates the spot-rates for certain maturity dates and a parameter vector
  according to Svensson.}
\usage{
svensson(beta, m)
}

\arguments{
  \item{beta}{ a vector of parameters \eqn{(\beta_0, \beta_1,\beta_2,\tau_1, \beta_3, \tau_2)} }
  \item{m}{ one maturity (or a vector of maturities) }
}
%\details{}
\value{
returns a vector consisting of the calculated spot rates 
}
\references{
Lars E.O. Svensson (1994):
Estimating and interpreting forward interest rates: Sweden 1992 -1994. 
\emph{Technical Reports 4871, National Bureau of Economic Research}. }

\author{Robert Ferstl, Josef Hayden }

\examples{
 svensson(c(0.07,0.3,0.05,0.1,0.08,0.2),1:30)
}
\keyword{}