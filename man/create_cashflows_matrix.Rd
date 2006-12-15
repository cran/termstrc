\name{create_cashflows_matrix}
\alias{create_cashflows_matrix}

\title{Cashflows Matrix Creation}
\description{
  Creates a matrix of cashflows for a specified group of  bonds. The number of rows is the number of cashflows for the bond with the longest maturity.

}
\usage{
create_cashflows_matrix(group, include_price = FALSE)
}

\arguments{
  \item{group}{group name}
  \item{include_price}{if \code{TRUE} the dirty price is included }
}


\value{matrix}

\author{Robert Ferstl, Josef Hayden}
%\note{ }
\seealso{\code{\link{create_maturities_matrix}} }
%\examples{}
\keyword{ }

