\name{splines_estim}
\alias{splines_estim}

\title{Term Structure and Credit Spread Estimation with Cubic Splines Method}
\description{
  Term structure and credit spread estimation with cubic splines method
}
\usage{
splines_estim(group,
              bonddata,
              matrange = "all")
}

\arguments{
  \item{group}{vector defining the group of bonds used for the estimation,\newline
               e.g. \code{c("GERMANY","AUSTRIA")}}
  \item{bonddata}{a dataset of bonds in list format}
  \item{matrange}{use \code{"all"} for no restrictions, or restrict the
    maturity range used for the estimation with \code{c(lower,upper)} }
}
\details{
  \item{\code{group}}{The first element of the vector will be used as the
   reference country for the credit spread estimation. The group can be either a vector of groups or a scalar.}
\item{\code{bonddata}}{The package is tailored to the included data sets. Therefore the structure and the naming convention of other used data sets has to be identical.  
 Use the function \code{str()} to explore the structure of the provided datasets.}
}
\value{
  The function \code{splines_estim} returns a list with the following elements or sub-lists:

  \item{group}{groups used from data set}
  \item{matrange}{maturity range}
  \item{n_group}{the number of groups used for the optimisation}
  \item{zcy_curves}{values for plotting the estimated zero-coupon yield curves}
  \item{scurves}{values for plotting the spread curves}
  \item{cf}{cashflows matrix for all specified groups}
  \item{m}{maturity matrix for all specified groups}
  \item{duration}{duration, weighted duration and duration based weights}
  \item{p}{dirty prices}
  \item{phat}{estimated bond prices}
  \item{y}{bond yields}
  \item{yhat}{theoretical bond yields calculated with the estimated
    bond prices \code{phat}}
  \item{alpha}{OLS coefficients of cubic splines estimation}
}
\references{
David Bolder and David Streliski (1999):
Yield Curve Modelling at the Bank of Canada.
\emph{Technical Report No 84 Bank of Canada} 


J.Huston McCulloch (1971):
Measuring the Term Structure of Interest Rates.
\emph{The Journal of Business, \bold{44} 19--31.} 


J. Huston McCulloch (1975):
The Tax-Adjusted Yield Curve.
\emph{The Journal of Finance, \bold{30} 811--830.} 


Sanjay K. Nawalkha and Gloria M. Soto and Natalia K. Beliaeva (2005):
Interest Rate Risk Modeling : The Fixed Income Valuation Course
\emph{Wiley Finance,60--67}

}

\author{Robert Ferstl, Josef Hayden}

\seealso{for another estimation method see \code{\link{nelson_estim}}}

\examples{
demo(euro02)
}

\keyword{}

