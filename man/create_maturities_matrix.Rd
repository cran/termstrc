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
  \item{group}{name of the group of bonds (default: \code{FALSE}).}
  \item{include_price}{if \code{TRUE} the dirty price is included.}
}


\value{The maturity matrix is returned. 
}

%\references{  }
%\author{Robert Ferstl, Josef Hayden}

%\note{ }

\seealso{\code{\link{create_cashflows_matrix}}}

%\examples{}
\keyword{ }
