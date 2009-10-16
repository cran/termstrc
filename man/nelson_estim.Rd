\name{nelson_estim}
\alias{nelson_estim}
\title{Zero-coupon Yield Curve Estimation with the Nelson/Siegel, Svensson Method}
\description{
 Zero-coupon yield curve estimation with the parametric Nelson and Siegel (1987) method and the Svensson (1994) extension.
}
\usage{
nelson_estim(group,
             bonddata,
             matrange = "all",
             method = "Nelson/Siegel",
             fit = "prices",
             weights = "duration",
             startparam,
             control = list(eval.max = 1000, iter.max= 500))
}
\arguments{
  \item{group}{vector defining the group of bonds used for the estimation,
               e.g. \command{c("GERMANY","AUSTRIA")}.}
  \item{bonddata}{a data set of bonds in list format.}
  \item{matrange}{use \command{"all"} for no restrictions, or restrict the
    maturity range used for the estimation with \command{c(lower,upper)}. }
  \item{method}{\command{"Nelson/Siegel"} or \command{"Svensson"}. }
  \item{fit}{use \command{"prices"} (\command{"yields"}) for minimizing the squared pricing (yield) errors.}
  \item{weights}{weighting of the errors used in the optimization (\command{"none"} or \command{"duration"}). }
  \item{startparam}{matrix of start parameters, for the Nelson/Siegel (Svensson) method 4 (6) parameters for each each group are required (one row per group).}
  \item{control}{list of control parameters for the function \command{nlminb()}.}
}
\details{
\itemize{
 \item{\command{group}}{The first element of the vector will be used as the
  			    	  	reference country for the spread curve calculation.
					 \command{group} can be either a vector of bond groups or a scalar.}
\item{\command{bonddata}}{The package is designed to work with a certain list data structure. For more information use the function \command{str()} to explore the structure of the example data sets.}

 \item{\command{startparam}}{A matrix of start parameters must be provided for the estimation. Due to the possible existence of multiple
 					 local minima, a great effort should be invested in finding a good start parameter set. For rolling estimation procedures, the previous optimal parameters should be used as inputs for the next iteration.}

 \item{\command{control}}{Please refer to the documentation of the function \command{\link{nlminb}}.}
}}
\value{
 The function \command{nelson_estim} returns an object of the class \command{"nelson"}. The object
contains the following items (mainly lists):
  \item{group}{group of bonds (e.g. countries) used for the estimation.}
  \item{matrange}{\command{"none"} or a vector with the maturity range.}
  \item{method}{estimation method (\command{"Nelson/Siegel"} or \command{"Svensson"}).}
  \item{fit}{objective function (\command{"prices"}, or \command{"yields"}).}
  \item{weights}{weighting of the errors used in the optimization (\command{"none"} or \command{"duration"}).}
  \item{n_group}{length of object \command{group}, i.e. the number of countries.}
 \item{spot}{zero-coupon yield curves as object of the class \command{"spot_curves"}.}
  \item{spread}{spread curves as object of the class \command{"s_curves"}.}
  \item{forward}{forward curves as object of the class \command{"fwr_curves"}.} 
 \item{discount}{discount curves as object of the class \command{"df_curves"}.} 
 \item{expoints}{extrapolation points for Nelson/Siegel method.}
 \item{cf}{cashflow matrices.}
  \item{m}{maturity matrices.}
  \item{p}{dirty prices.}
  \item{phat}{estimated bond prices.}
  \item{perrors}{pricing errors and maturities as object of the class \command{"error"}.}
  \item{y}{bond yields.}
  \item{yhat}{one list for each group with the theoretical bond yields calculated with the estimated bond prices \command{phat}.}
  \item{yerrors}{yield errors and maturities as object of the class \command{"error"}.}
  \item{opt_result}{optimization results from \command{nlminb}, e.g. optimal parameters, convergence info.}
   
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


\note{An error message concerning the function \command{uniroot()} is in general caused by wrongly specified start parameters.
 
For objects of the class \command{"spot_curves"}, \command{"s_curves"}, \command{"df_curves"}, \command{"fwr_curves"}, \command{"error"} appropriate plot methods are offered.  For objects of the class \command{"nelson"} print, summary and plot methods are available. Another term structure estimation method is provided by the function \command{\link{splines_estim}}.} 

\seealso{\command{\link{print.nelson}}, \command{\link{summary.nelson}},  \command{\link{plot.nelson}},  \command{\link{splines_estim}},  \command{\link{plot.spot_curves}},  \command{\link{plot.s_curves}},  \command{\link{plot.df_curves}},  \command{\link{plot.fwr_curves}},  \command{\link{plot.error}},  \command{\link{uniroot}}.}


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


