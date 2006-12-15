\name{nelson_estim}
\alias{nelson_estim}
\title{ Term Structure and Credit Spread Estimation with Nelson/Siegel and Svensson Method}
\description{
 Term structure and credit spread estimation with Nelson/Siegel and Svensson method
}
\usage{
nelson_estim(group,
             bonddata,
             matrange = "all",
             method = "Nelson/Siegel",
             fit = "prices",
             weights = "none",
             startparam, control = list(eval.max = 1000))
}
\arguments{
  \item{group}{vector defining the group of bonds used for the estimation,\newline
               e.g. \code{c("GERMANY","AUSTRIA")}}
  \item{bonddata}{a dataset of bonds in list format}
  \item{matrange}{use \code{"all"} for no restrictions, or restrict the
    maturity range used for the estimation with \code{c(lower,upper)} }
  \item{method}{\code{"Nelson/Siegel"} or \code{"Svensson"} }
  \item{fit}{use \code{"prices"} (\code{"yields"}) for minimising the squared price (yield) error}
  \item{weights}{If a weighted minimisation of the squared price deviation is required, use \code{"duration"}, otherwise \code{"none"} }
  \item{startparam}{matrix of start parameters, for the Nelson/Siegel (Svensson) method 4 (6) parameters for each each group are required (one row per group).}
  \item{control}{list of control parameters for the function \code{nlminb}}
}
\details{
 \item{\code{group}}{The first element of the vector will be used as the
   reference country for the credit spread estimation. The group can be either a vector of groups or a scalar.}
\item{\code{bonddata}}{The package is tailored to the included data sets. Therefore the structure and the naming convention of other used data sets has to be identical.  
 Use the function \code{str()} to explore the structure of the provided datasets.}
 \item{\code{startparam}}{For the estimation a matrix of start parameters must be provided. The estimation is very sensitive to the start parameters, due to multiple
 local minima of the objective functions. Therefore a great effort
 should be invested in finding a good start parameter set. }

 \item{\code{control}}{Please refer to the documentation of the function \code{mlimnb}.}
}
\value{
The function \code{nelson_estim} returns a list with the following elements or sub-lists:

  \item{group}{groups used from data set}
  \item{matrange}{includes the chosen maturity range}
  \item{method}{includes the chosen estimation method}
  \item{fit}{includes the chosen objective function}
  \item{weights}{type of weighting used in optimisation}
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
  \item{opt_result}{the optimal parameter vector for the specified groups according to the chosen estimation approach}
   
}
\references{ 
Charles R. Nelson and Andrew F. Siegel (1987): 
Parsimonious modeling of yield curves.
\emph{The Journal of Business}, \bold{60(4)}:473--489.

Lars E.O. Svensson (1994): 
Estimating and interpreting forward interest rates: Sweden 1992 -1994. 
\emph{Technical Reports 4871, National Bureau of Economic Research}. 
 
David Bolder and David Streliski (1999):
Yield Curve Modelling at the Bank of Canada.
\emph{Technical Report No 84 Bank of Canada} 

 }
\author{Robert Ferstl, Josef Hayden}
\seealso{for another estimation method see \code{\link{splines_estim}}}
\examples{ 
demo(euro01)
}

\keyword{}

