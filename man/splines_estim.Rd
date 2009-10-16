\name{splines_estim}
\alias{splines_estim}

\title{Discount Curve Estimation with McCulloch Cubic Splines}
\description{
  Discount curve estimation with the cubic splines approach by McCulloch (1975). }

\usage{
splines_estim(group,
              bonddata,
              matrange = "all")
}

\arguments{
  \item{group}{vector defining the group of bonds used for the estimation, e.g. \command{c("GERMANY","AUSTRIA")}.}
  \item{bonddata}{a data set of bonds in list format.}
  \item{matrange}{use \command{"all"} for no restrictions, or restrict the
    maturity range used for the estimation with \command{c(lower,upper)}. }
}
\details{\itemize{
 \item{\command{group}}{The first element of the vector will be used as the
  			    	  	reference country for the spread curve calculation.
					 \command{group} can be either a vector of bond groups or a scalar.}
\item{\command{bonddata}}{The package is designed to work with a certain list data structure. For more information use the function \command{str()} to explore the structure of the example data sets.}
}}

\value{
  The function \command{splines_estim} returns an object of the class \command{"cubicsplines"}. The object
contains the following items (mainly lists):
  \item{group}{group of bonds (e.g. countries) used for the estimation.}
  \item{matrange}{\command{"none"} or a vector with the maturity range.}
  \item{n_group}{length of object \command{group}, i.e. the number of countries.}
  \item{knotpoints}{selected knot points for the cubic splines estimation.}
 \item{spot}{zero-coupon yield curves as object of the class \command{"spot_curves"}.}
  \item{spread}{spread curves as object of the class \command{"s_curves"}.}
  \item{forward}{forward curves as object of the class \command{"fwr_curves"}.} 
 \item{discount}{discount curves as object of the class \command{"df_curves"}.} 
 \item{cf}{cashflow matrices.}
  \item{m}{maturity matrices.}
  \item{p}{dirty prices.}
  \item{phat}{estimated bond prices.}
  \item{perrors}{pricing errors and maturities as object of the class \command{"error"}.}
  \item{y}{bond yields.}
  \item{yhat}{one list for each group with the theoretical bond yields calculated with the estimated bond prices \command{phat}.}
  \item{yerrors}{yield errors and maturities as object of the class \command{"error"}.}
  \item{alpha}{OLS coefficients of cubic splines estimation.}
  \item{regout}{OLS estimation results as object of the class \command{"lm"}.}
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

\note{ For objects of the class \command{"spot_curves"}, 
\command{"s_curves"}, \command{"df_curves"}, \command{"fwr_curves"}, \command{"error"} appropriate plot methods are offered. For objects of the list item \command{regout} standard \command{lm} methods apply. For objects of the class \command{"cubicsplines"} print, summary and plot methods are available. Another term structure estimation method is provided by the function \command{\link{nelson_estim}}.} 

\seealso{\command{\link{print.cubicsplines}}, \command{\link{summary.cubicsplines}},  \command{\link{plot.cubicsplines}},  \command{\link{nelson_estim}},  \command{\link{plot.spot_curves}},  \command{\link{plot.s_curves}},  \command{\link{plot.df_curves}},  \command{\link{plot.fwr_curves}},  \command{\link{plot.error}},  \command{\link{summary.lm}},  \command{\link{print.lm}},  \command{\link{plot.lm}}.}

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



