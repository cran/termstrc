\name{nelson_estim}
\alias{nelson_estim}
\title{ Term Structure and Credit Spread Estimation with the Nelson/Siegel, Svensson method}
\description{
 Term structure and credit spread estimation according to Nelson and Siegel (1987), and Svensson (1994).
}
\usage{
nelson_estim(group,
             bonddata,
             matrange = "all",
             method = "Nelson/Siegel",
             fit = "prices",
             weights = "duration",
             startparam, control = list(eval.max = 1000, iter.max= 500))
}
\arguments{
  \item{group}{vector defining the group of bonds used for the estimation,\newline
               e.g. \code{c("GERMANY","AUSTRIA")}.}
  \item{bonddata}{a data set of bonds in a list format.}
  \item{matrange}{use \code{"all"} for no restrictions, or restrict the
    maturity range used for the estimation with \code{c(lower,upper)}. }
  \item{method}{\code{"Nelson/Siegel"} or \code{"Svensson"}. }
  \item{fit}{use \code{"prices"} (\code{"yields"}) for minimising the squared price (yield) error.}
  \item{weights}{If a weighted minimisation of the squared price deviation is required, use \code{"duration"}, otherwise \code{"none"}. }
  \item{startparam}{a matrix of start parameters, for the Nelson/Siegel (Svensson) method 4 (6) parameters for each each group are required (one row per group).}
  \item{control}{list of control parameters for the function \code{nlminb}.}
}
\details{
 \item{\code{group}}{	The first element of the vector will be used as the
  			    	  	reference country for the spread curve calculation.
					 \code{group} can be either a vector of bond groups or a scalar.}
\item{\code{bonddata}}{	The package is tailored to a certain structure. Therefore the structure
			    	        and the naming convention of other used data sets has to be identical to the provided one.  
 					Use the function \code{str()} to explore the structure of the provided data sets.}

 \item{\code{startparam}}{For the estimation a matrix of start parameters must be provided.
 			      	  	 The estimation is very sensitive to the start parameters, due to multiple
 					 local minima of the objective functions. Therefore a great effort
 					 should be invested in finding a good start parameter set.}

 \item{\code{control}}{	Please refer to the documentation of the function \code{mlimnb}.}
}
\value{
 The function \code{nelson_estim} returns an object of the class \code{"nelson"}. The object
contains the following items (mainly lists):

  \item{group}{vector consisting of the bond groups used for the estimation.}
  \item{matrange}{either \code{"none"} or a vector with the maturity range.}
  \item{method}{the chosen estimation method (\code{"Nelson/Siegel"} or \code{"Svensson"}) .}
  \item{fit}{the chosen objective function (\code{"prices"}, or \code{"yields"}).}
  \item{weights}{type of weighting used in the optimisation (\code{"none"} or \code{"duration"}).}
  \item{n_group}{the number of groups used for the optimisation.}
 \item{spot}{one list for each group with an object of the class \code{"spot_curves"}.}
  \item{spread}{one list for each group with  an object of the class \code{"s_curves"}.}
  \item{forward}{one list for each group with an object of the class \code{"fwr_curves"}.} 
 \item{discount}{one list for each group with  an object of the class \code{"df_curves"}.} 
 \item{expoints}{one list for each group with the extrapolation points.}
 \item{cf}{one list for each group with the cashflows matrix.}
  \item{m}{one list for each group with the maturity matrix.}
  \item{p}{one list for each group with the dirty prices.}
  \item{phat}{one list for each group with the estimated bond prices.}
  \item{perrors}{one list for each group with an object of the class \code{"error"} containing the price errors of the bonds and the associated maturities.}
  \item{y}{one list for each group with the bond yields and the associated maturities.}
  \item{yhat}{one list for each group with the theoretical bond yields calculated with the estimated bond prices \code{phat}.}
  \item{yerrors}{ one list for each group with an object of the class \code{"error"} containing the yield errors of the bonds and the associated maturities.}
  \item{opt_result}{one liste for each group with the optimal parameter vector obtained from the optimisation and further information concerning the optimisation.}
   
}
\references{ 
Charles R. Nelson and Andrew F. Siegel (1987): 
Parsimonious Modeling of Yield Curves.
\emph{The Journal of Business}, \bold{60(4)}:473--489.

Lars E.O. Svensson (1994): 
Estimating and Interpreting Forward Interest Rates: Sweden 1992 -1994. 
\emph{Technical Reports 4871, National Bureau of Economic Research}. 
 }
%\author{Robert Ferstl, Josef Hayden}


\note{An error message concerning the function \code{uniroot} is in general due to wrongly specified start parameters.
 
For objects of the class \code{"spot_curves"}, \code{"s_curves"}, \code{"df_curves"}, \code{"fwr_curves"}, \code{"error"} appropriate plot methods are offered.  For objects of the class \code{"nelson"} print, summary and plot methods are available. Another term structure estimation method is provided by the function \code{\link{splines_estim}}.} 

\seealso{\code{\link{print.nelson}}, \code{\link{summary.nelson}},  \code{\link{plot.nelson}},  \code{\link{splines_estim}},  \code{\link{plot.spot_curves}},  \code{\link{plot.s_curves}},  \code{\link{plot.df_curves}},  \code{\link{plot.fwr_curves}},  \code{\link{plot.error}},  \code{\link{uniroot}}.}


\examples{ 
# load data set
data(eurobonds)

# define countries, for which the estimation 
# of the zero-coupon yield curves will be carried out
group <- c("GERMANY", "AUSTRIA", "ITALY")

# define data set
bonddata <- eurobonds

# set maturity range
matrange <- c(2,12)

# define the used method
method <- "Nelson/Siegel"

# the weighted squared price errors will be minimized
fit <- "prices"
weights <- "duration"

# naive start parameters
b <- matrix(rep(c(0,0,0, 1),3),nrow=3,byrow=TRUE)
rownames(b) <- group
colnames(b) <- c("beta0","beta1","beta2","tau1")

# perform estimation
x <- nelson_estim(group, bonddata, matrange, 
                  method, fit, weights, startparam=b)

# prints the obtained parameters of the estimation
print(x)

# goodness of fit measures
summary(x)

# plots the zero-coupon yield curves for each country 
plot(x,errors="none")

# plots all zero-coupon yield curves togther
plot(x,multiple=TRUE, errors="none")

# spread curve splot
plot(x,ctype="spread",errors="none")

# price error plot for all countries
plot(x,ctype="none")

}

\keyword{}

