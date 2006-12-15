\name{corpbonds}
\alias{corpbonds}
\docType{data}
\title{Corporate Bonds}
\description{Corporate bonds}
\usage{data(corpbonds)}

\details{
The data set eurobonds consists of bonds of the rating classes AAA ,AA+,
AA, AA-, A+, A, A-, BBB+, BBB, BBB-
}

\note{
If you use your own dataset, make sure that the structure is identical to the provided data sets.
Use the function \code{str()} to explore the data set. 
 
Every element of the list except \code{$RATING, $NAME} is required for the estimation. The elements 
\code{$COUPONRATE, $PRICE} are not necessarily required. The calculation of the accrued interest 
is not implemented and has therefore be provided by the user. For the provided datasets the list element
\code{$ACCRUED} includes the accrued interests. 
}

\seealso{ \code{\link{eurobonds}}}

\examples{
data(corpbonds)
str(corpbonds) 


# You may use this to generate an empty data set,
# which can be filled with bond data:

ISIN <- vector()
MATURITYDATE <- vector()
STARTDATE <- vector()
COUPONRATE <- vector()
PRICE <- vector()
ACCRUED <- vector()

CFISIN <- vector()
CF <- vector()
DATE <- vector()

CASHFLOWS <- list(CFISIN,CF,DATE)
names(CASHFLOWS) <- c("ISIN","CF","DATE")

TODAY <- vector()

mycountry1 <- list(ISIN,MATURITYDATE,STARTDATE,
                   COUPONRATE,PRICE,ACCRUED,CASHFLOWS,TODAY)
mycountry2 <- list(ISIN,MATURITYDATE,STARTDATE,
                   COUPONRATE,PRICE,ACCRUED,CASHFLOWS,TODAY)

names(mycountry1) <- c("ISIN","MATURITYDATE","STARTDATE","COUPONRATE",
                       "PRICE","ACCRUED","CASHFLOWS","TODAY")
names(mycountry2) <- c("ISIN","MATURITYDATE","STARTDATE","COUPONRATE",
                       "PRICE","ACCRUED","CASHFLOWS","TODAY")

mybonds <- list(mycountry1,mycountry2)

names(mybonds) <- c("mycountry1","mycountry2")

}
\keyword{datasets}
