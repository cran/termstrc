\name{splines_estim}
\alias{splines_estim}

\title{Term Structure and Credit Spread Estimation with Cubic Splines}
\description{
  Term structure and credit spread estimation using an cubic splines approach of McCulloch(1975). }

\usage{
splines_estim(group,
              bonddata,
              matrange = "all")
}

\arguments{
  \item{group}{vector defining the group of bonds used for the estimation,\newline
               e.g. \code{c("GERMANY","AUSTRIA")}.}
  \item{bonddata}{a data set of bonds in a list format.}
  \item{matrange}{use \code{"all"} for no restrictions, or restrict the
    maturity range used for the estimation with \code{c(lower,upper)}. }
}
\details{
  \item{\code{group}}{	The first element of the vector will be used as the reference country for the spread curve calculation. The group can be either a vector of groups or a scalar.}

\item{\code{bonddata}}{	The package is tailored to a certain structure. Therefore the structure and 
			    	       	the naming convention of other used data sets has to be identical to the provided ones.  
 					Use the function \code{str()} to explore the structure of the provided data sets.}
}

\value{
  The function \code{splines_estim} returns an object of the class \code{"cubicsplines"}. The object
contains the following items (mainly lists):

  \item{group}{vector consisting of the bond groups used for the estimation.}
  \item{matrange}{either \code{"none"} or a vector with the maturity range.}
  \item{n_group}{the number of groups used for the optimisation.}
  \item{knotpoints}{one list for each group with the used knot points for the cubic splines estimation.}
  \item{spot}{one list for each group with an object of the class \code{"spot_curves"}.}
  \item{spread}{one list for each group with  an object of the class \code{"s_curves"}.}
 \item{discount}{one list for each group with  an object of the class \code{"df_curves"}.} 
\item{forward}{one list for each group with an object of the class \code{"fwr_curves"}.} 
 \item{cf}{one list for each group with the cashflows matrix.}
  \item{m}{one list for each group with the maturity matrix.}
  \item{p}{one list for each group with the dirty prices.}
  \item{phat}{one list for each group with the estimated bond prices.}
  \item{perrors}{one list for each group with an object of the class \code{"error"} containing the price errors of the bonds and the associated maturities.}
  \item{y}{one list for each group with the bond yields and the associated maturities.}
  \item{yhat}{one list for each group with the theoretical bond yields calculated with the estimated bond prices \code{phat}.}
  \item{yerrors}{ one list for each group with an object of the class \code{"error"} containing the yield errors of the bonds and the associated maturities.}
  \item{alpha}{OLS coefficients of cubic splines estimation.}
  \item{regout}{one list for each group with an object of the class \code{"lm"}.}
}
\references{
J.Huston McCulloch (1971):
Measuring the Term Structure of Interest Rates.
\emph{The Journal of Business, \bold{44} 19--31.} 

J. Huston McCulloch (1975):
The Tax-Adjusted Yield Curve.
\emph{The Journal of Finance, \bold{30} 811--830.} 
}

%\author{Robert Ferstl, Josef Hayden}

\note{ For objects of the class \code{"spot_curves"}, 
\code{"s_curves"}, \code{"df_curves"}, \code{"fwr_curves"}, \code{"error"} appropriate plot methods are offered. For objects of the list item \code{regout} standard \code{lm} methods apply. For objects of the class \code{"cubicsplines"} print, summary and plot methods are available. Another term structure estimation method is provided by the function \code{\link{nelson_estim}}.} 

\seealso{\code{\link{print.cubicsplines}}, \code{\link{summary.cubicsplines}},  \code{\link{plot.cubicsplines}},  \code{\link{nelson_estim}},  \code{\link{plot.spot_curves}},  \code{\link{plot.s_curves}},  \code{\link{plot.df_curves}},  \code{\link{plot.fwr_curves}},  \code{\link{plot.error}},  \code{\link{summary.lm}},  \code{\link{print.lm}},  \code{\link{plot.lm}}.}

\examples{
# load data set
data(eurobonds)

# define countries, for which the estimation 
# of the zero-coupon yield curves will be carried out
group <- c("GERMANY", "AUSTRIA", "ITALY")

# define data set
bonddata <- eurobonds

# set maturtiy range
matrange <- c(0, 19)  

# perform estimation
x <- splines_estim(group, bonddata, matrange)

# print the obtained parameters of the estimation
print(x)

# goodness of fit measures
summary(x)

# plot the zero-coupon yield curve for each country
plot(x,errors="none")

# plot all zero-coupon yield curves together
plot(x,multiple=TRUE,errors="none")

# spread curve splot
plot(x,ctype="spread",errors="none")

# price error plot for all countries
plot(x,ctype="none")
}

\keyword{}

