\name{bond_prices}
\alias{bond_prices}

\title{Bond Price Calculation }
\description{
  Function for the calculation of bond prices according to the chosen approach 
  (Nelson and Siegel or Svensson) based on the cashflows and maturities of the bonds. 
  
}
\usage{
bond_prices(method = "Nelson/Siegel", beta, m, cf)
}

\arguments{
  \item{method}{ defines the desired method,  \code{"Nelson/Siegel"} for the Nelson/Siegel approach
  		or \code{"Svensson"} for the Svensson approach}
  \item{beta}{parameter vector, is linked to the chosen approach}
  \item{m}{maturities matrix, consists
  of the maturity dates which are appended to the cashflows of the bonds   }
  \item{cf}{ cashflows matrix}
 }
%\details{}
\value{
  
  Returns a list with:
  
  \item{spot_rates}{spot rates}
  \item{discount_factors}{discount factors}
  \item{bond_prices}{bond prices}
}

\references{David Bolder and David Streliski (1999): Yield Curve Modelling at the Bank of Canada.\emph{Technical Report No 84 Bank of Canada}}

\author{ Robert Ferstl, Josef Hayden }
\note{  


}
\seealso{ \code{\link{svensson}}, \code{\link{nelson_siegel}}}

%\examples{}

\keyword{  }