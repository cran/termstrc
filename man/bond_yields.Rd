\name{bond_yields}
\alias{bond_yields}

\title{Bond Yield Calculation}

\description{Function for the calculation of bond yields.}

\usage{
bond_yields(cashflows, m, searchint = c(-1, 1), tol = 1e-10)
}

\arguments{
  \item{cashflows}{ matrix with the bonds cashflows.}
  \item{m}{maturity matrix.}
  \item{searchint}{search interval for root finding.}
  \item{tol}{desired accuracy for function \command{uniroot}.}
}
%\details{}

\value{
The function returns a matrix with the bond yields and the associated maturities. 
}

%\references{}

%\author{Robert Ferstl, Josef Hayden }
%\note{}

\seealso{ \command{\link{uniroot}}}
%\examples{}

