\name{impl_fwr}
\alias{impl_fwr}

\title{Implied Forward Rate Calculation}
\description{Calculates the implied forward rates from given spot rates.
}
\usage{
impl_fwr(m, s)
}


\arguments{
  \item{m}{maturity vector.}
  \item{s}{spot rate vector.}
}
\details{Implied forward rates can be calculated using the following relationship:

 \deqn{f(t',T) = \frac{s(m_T)m_T - s(m_{t'})m_{t'}}{m_T-m_{t'}},}

whereas \eqn{s(m_T), s(m_{t'})} is the spot rate for a maturity \eqn{m_T,m_{t'}} respectively.}
\value{
The function returns the calculated forward rate vector.
}
%\references{  }
%\author{ Robert Ferstl, Josef Hayden }
%\note{ }
%\seealso{}
%\examples{}

\keyword{  }

