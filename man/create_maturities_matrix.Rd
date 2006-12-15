\name{create_maturities_matrix}
\alias{create_maturities_matrix}

\title{Maturity Matrix Creation }
\description{
  Creates a matrix of maturities for a specified group of bonds. The number of rows is the number of cashflows for the bond with the longest maturity.

}

\usage{
create_maturities_matrix(group, include_price = FALSE)
}

\arguments{
  \item{group}{group name}
  \item{include_price}{if TRUE the dirty price is included}
}


\value{matrix
}

%\references{ ~put references to the literature/web site here ~ }
\author{Robert Ferstl, Josef Hayden}

%\note{ ~~further notes~~ 

\seealso{\code{\link{create_cashflows_matrix}}}

%\examples{}
\keyword{ }
