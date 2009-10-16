\name{spotrates}
\alias{spotrates}

\title{Function for the Calculation of the Spot Rates }
\description{
 The function calculates the spot rates for the chosen approach, a provided maturity vector and
 parameter set.
}

\usage{
spotrates(method, beta, m)
}



\arguments{
  \item{method}{\command{"Nelson/Siegel"} or \command{"Svensson"}.}
  \item{beta}{ parameter set \eqn{\bm{\beta}}.}
  \item{m}{maturity or a vector of maturities.}
}
\details{The spot rates according to Nelson/Siegel are defined as:
\deqn{
  s(m,\bm{\beta}) = \beta_0 + \beta_1\frac{1-\exp(-\frac{m}{\tau_1})}{\frac{m}{\tau_1}} + \beta_2\left(\frac{1-\exp(-\frac{m}{\tau_1})}{\frac{m}{\tau_1}} - \exp(-\frac{m}{\tau_1})\right).}

Svensson defines the spot rate function as follows:
\deqn{
s(m,\bm{\beta}) = \beta_0 + \beta_1\frac{1-\exp(-\frac{m}{\tau_1})}{\frac{m}{\tau_1}} + \beta_2\left(\frac{1-\exp(-\frac{m}{\tau_1})}{\frac{m}{\tau_1}} - \exp(-\frac{m}{\tau_1})\right) + \beta_3\left(\frac{1-\exp(-\frac{m}{\tau_2})}{\frac{m}{\tau_2}} - \exp(-\frac{m}{\tau_2})\right)
}

}

\value{Returns a vector with the calculated spot rates.}

\references{
Charles R. Nelson and Andrew F. Siegel (1987): 
Parsimonious Modeling of Yield Curves.
\emph{The Journal of Business}, \bold{60(4)}:473--489.


Lars E.O. Svensson (1994): 
Estimating and Interpreting Forward Interest Rates: Sweden 1992 -1994. 
\emph{Technical Reports 4871, National Bureau of Economic Research}.  }
%\author{Robert Ferstl, Josef Hayden }

%\note{  }
%\seealso{ }

%\examples{}

