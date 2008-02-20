###################################################################
#            cashflows matrix for a given group                   #
###################################################################

create_cashflows_matrix <- function(group,include_price=FALSE) {
  
  n_of_cf <- summary(factor(group$CASHFLOWS$ISIN,levels=group$ISIN))
  n_of_bonds <- length(n_of_cf)
  max_cf <- max(n_of_cf)
  pos_cf <- c(0,cumsum(n_of_cf))
  
  # the number of rows of the matrix is the number of 
  # cashflows of the bond with the maximum of cashflows
  # all missing elements of the matrix are filled up with zeros 
  	
  CASHFLOWMATRIX <-
    mapply(function(i) c(group$CASHFLOWS$CF[(pos_cf[i]+1):pos_cf[i+1]],
                         rep(0,max_cf-n_of_cf[i])),
           1:n_of_bonds)
  
  if (include_price == TRUE) {CASHFLOWMATRIX <- rbind(-(group[["PRICE"]]
      +group[["ACCRUED"]]),CASHFLOWMATRIX)}              
  colnames(CASHFLOWMATRIX) <- group$ISIN
  CASHFLOWMATRIX
}

###################################################################
#              maturities matrix for a given group                #
###################################################################

create_maturities_matrix <-
  function(group,include_price=FALSE) {

  n_of_cf <- summary(factor(group$CASHFLOWS$ISIN,levels=group$ISIN))
  n_of_bonds <- length(n_of_cf)
  max_cf <- max(n_of_cf)
  pos_cf <- c(0,cumsum(n_of_cf))
  year_diff <- as.numeric(difftime(as.Date(group$CASHFLOW$DATE),
  				as.Date(group$TODAY),units="days"))/365
 
  # the number of rows of the matrix is the number of 
  # maturity dates of the bond with the longest maturity
  # all missing elements of the matrix are filled up with zeros 
  
  MATURITYMATRIX <-
     mapply(function(i) c(year_diff[(pos_cf[i]+1):pos_cf[i+1]],
                     rep(0,max_cf-n_of_cf[i])),
            1:n_of_bonds)
  
  if (include_price == TRUE) {MATURITYMATRIX <- rbind(rep(0,n_of_bonds),
                           MATURITYMATRIX)}  
  colnames(MATURITYMATRIX) <- group$ISIN
  MATURITYMATRIX
}

###################################################################
#              bond yields for a given group                      #
###################################################################

bond_yields <- function(cashflows, m, searchint=c(-1,1), tol=1e-10) {

  # convert input data to matrices if necessary
  if (!is.matrix(cashflows))
    cashflows <- as.matrix(cashflows)
  if (!is.matrix(m))
    m <- as.matrix(m)

  # create empty bond yields matrix in appropriate size
  bondyields <- matrix(0, nrow=ncol(cashflows), ncol=2)                                                                                                                     

  # put maturity of every bond into first column of bond yields matrix
  bondyields[,1] <- apply(m, 2, max)

  # traverse list of bonds
  for (i in seq_len(ncol(cashflows))) {
                                                  
    # present value of cash flows for root finding 
    pvcashflows <- function(y) {
       t(cashflows[,i])%*%exp(-m[,i]*y)
    }

    # calculate bond yields
    
    bondyields[i,2] <- uniroot(pvcashflows, searchint, tol = tol,maxiter=3000)$root 
  }

  # return calculated bond yields matrix
  rownames(bondyields) <- colnames(cashflows)
  colnames(bondyields) <- c("Maturity","Yield")
  bondyields
}

###################################################################
#                      maturity range                             #
###################################################################

maturity_range <- 
  function(bonddata,lower,upper) {

m <- lapply(bonddata,create_maturities_matrix)
colmax <- function(m) apply(m,2,max)
m_max <- lapply(m,colmax)

bonds_in_range <- function(group) names(group[which(
                group>lower & group<upper)])

# list with ISINs of bonds in range
isins_range <- lapply(m_max,bonds_in_range)
index_set <- which(unlist(lapply(isins_range,length))>0)

# list with positions of bonds in isins_range
isins_range_pos <- list()

for (i in seq_along(bonddata)) {
        isins_range_pos[[i]] <- which(bonddata[[i]]$ISIN %in% 
            isins_range[[i]])
    }
names(isins_range_pos) <- names(bonddata)

# first part of bonddata for filtering

N <- which(names(bonddata[[1]]) %in% c("ISIN","MATURITYDATE","STARTDATE","COUPONRATE","PRICE","ACCRUED"))

first <- function(lst) lst[N]
filtered <- lapply(bonddata,first)

bonddata_range <- list()
for (i in seq_along(bonddata)) {
bonddata_range[[i]] <- as.list(as.data.frame(filtered[[i]])
                       [isins_range_pos[[i]],])
# convert to character                       
bonddata_range[[i]][["ISIN"]] <- as.character(bonddata_range[[i]][["ISIN"]])
bonddata_range[[i]][["MATURITYDATE"]] <- as.character(bonddata_range[[i]][["MATURITYDATE"]])
bonddata_range[[i]][["STARTDATE"]] <- as.character(bonddata_range[[i]][["STARTDATE"]])
}
names(bonddata_range) <- names(bonddata)

# list with positions of cashflows in isins_range

isins_range_pos <- list()
for (i in seq_along(bonddata)) {
isins_range_pos[[i]] <- which(bonddata[[i]][["CASHFLOWS"]]
                        [["ISIN"]]%in%isins_range[[i]])
}
names(isins_range_pos) <- names(bonddata)

for (i in seq_along(bonddata)) {
CASHFLOWS <- as.list(as.data.frame(bonddata[[i]]
             [["CASHFLOWS"]])[isins_range_pos[[i]],])
CASHFLOWS$ISIN <- as.character(CASHFLOWS$ISIN)
CASHFLOWS$DATE <- as.character(CASHFLOWS$DATE)             
bonddata_range[[i]][["CASHFLOWS"]] <- list()
bonddata_range[[i]][["CASHFLOWS"]] <- CASHFLOWS
}
names(bonddata_range) <- names(bonddata)
bonddata_range

# add TODAY from bonddata
for (i in seq_along(bonddata)) {
bonddata_range[[i]][["TODAY"]] <- bonddata[[i]][["TODAY"]]
}

# delete countries where no bonds are available
bonddata_range <- bonddata_range[index_set]
bonddata_range
}
   
###################################################################
#                 Duration                                        #
###################################################################

# duration, modified duration and weights for optimization


duration <-
function (cf_p,m_p,y) {
       y <- matrix(rep(y,nrow(m_p)),ncol=ncol(m_p),byrow=TRUE)
       # mac cauly duration
       d <- apply(cf_p*m_p*exp(-y*m_p),2,sum)/-cf_p[1,]
       # modified duration
       md <- d/(1+y[1,])
       omega <- (1/d)/sum(1/d)
       dur <- cbind(d,md,omega)
       colnames(dur) <- c("Duration","Modified duration","Weights")
       dur
    }


###################################################################
#                 Spotrate calculation                            #
###################################################################

 spotrates <- function(method,beta,m){ 
  switch(method,
 	"Nelson/Siegel" = nelson_siegel(beta,m),
 	"Svensson" = svensson(beta,m))
  }
	
###################################################################
#                 Forwardrate calculation                         #
###################################################################

 forwardrates <- function(method,beta,m){ 
  switch(method,
 	"Nelson/Siegel" = fwr_ns(beta,m),
 	"Svensson" = fwr_sv(beta,m))
  }
  
###################################################################
#              Implied forward rate calculation                   #
###################################################################

impl_fwr <- function(m,s) {
	
impl_fwr <- c(s[1],(s[-1]*m[-1] - s[-length(s)]*m[-length(m)])/(diff(m)))
impl_fwr[1] <- impl_fwr[2]
impl_fwr	
	
	}
  	
	
###################################################################
#                   Bond pricing function                         #
###################################################################

bond_prices <-
  function(method="Nelson/Siegel", beta, m, cf) {
     
  # calculate spot rates
  spot_rates <- spotrates(method,beta,m)

  # replace NaNs by zeros
  spot_rates[is.nan(spot_rates)] <- 0        
     
  # calculate discount factors
  discount_factors <- exp(-m*spot_rates)

  # calculate bond prices
  bond_prices <- apply(cf*discount_factors, 2, sum)  
  
  # return spot rates, discount factors and bond prices
  return (list(spot_rates=spot_rates,
               discount_factors=discount_factors,
               bond_prices=bond_prices))
  
}  

###################################################################
#                      Root mean squared error                    #
###################################################################

rmse <-
function (actual,estimated) {
	e <- actual - estimated
	sqrt(1/length(e)*sum((e-mean(e))^2))			
      }

###################################################################
#                 Average absolute error                          #
###################################################################
      
aabse <-
function (actual,estimated){
     e <- actual - estimated	
     1/length(e)*sum(abs(e-mean(e)))
     }   							

###################################################################
#                    Cubic functions                              #
###################################################################

gi <-
function(t,T,i,s){
  g <- rep(0,length(t))
  for(j in seq_along(t)){
    if(i==1){
    if(T[i]<=t[j]&t[j]<T[i+1]){
     g[j] <- (T[i])^2/6 + ((T[i])*(t[j]-T[i]))/2 + (t[j]-T[i])^2/2 - (t[j]-T[i])^3/(6*(T[i+1]-T[i]))
    }
    if(t[j]>=T[i+1]){
     g[j] <- (T[i+1])*((2*T[i+1]-T[i])/6 + (t[j]-T[i+1])/2)
    }   
  }
  if(i>1&i<length(T)){
    if(t[j]<T[i-1]){
     g[j] <- 0
    }
    if(T[i-1]<=t[j]&t[j]<T[i]){
     g[j] <- (t[j]-T[i-1])^3/(6*(T[i]-T[i-1]))
    }
    if(T[i]<=t[j]&t[j]<T[i+1]){
     g[j] <- (T[i]-T[i-1])^2/6 + ((T[i]-T[i-1])*(t[j]-T[i]))/2 + (t[j]-T[i])^2/2 - (t[j]-T[i])^3/(6*(T[i+1]-T[i]))
    }
    if(t[j]>=T[i+1]){
     g[j] <- (T[i+1]-T[i-1])*((2*T[i+1]-T[i]-T[i-1])/6 + (t[j]-T[i+1])/2)
    }
  }
   if(i==length(T)){
    if(t[j]<T[i-1]){
     g[j] <- 0
    }
    if(T[i-1]<=t[j]&t[j]<=T[i]){
     g[j] <- (t[j]-T[i-1])^3/(6*(T[i]-T[i-1]))
    }
  } 
  if(i==s){
    g[j] <- t[j]
  }
}
  g
}

###################################################################
#                    Bond removal function                        #
###################################################################

rm_bond <- function(bdata,ISIN,gr){
    cf_isin_index <- which(bdata[[gr]]$CASHFLOWS$ISIN %in% ISIN)
 	isin_index <- which(bdata[[gr]]$ISIN %in% ISIN)	

    	bdata[[gr]]$ISIN <-  bdata[[gr]]$ISIN[-isin_index]
    	bdata[[gr]]$MATURITYDATE <- bdata[[gr]]$MATURITYDATE[-isin_index]
    	bdata[[gr]]$STARTDATE <- bdata[[gr]]$STARTDATE[-isin_index]
    	bdata[[gr]]$COUPONRATE <- bdata[[gr]]$COUPONRATE[-isin_index]
    	bdata[[gr]]$PRICE <- bdata[[gr]]$PRICE[-isin_index]
    	bdata[[gr]]$ACCRUED <- bdata[[gr]]$ACCRUED[-isin_index]

		bdata[[gr]]$CASHFLOWS$ISIN <- bdata[[gr]]$CASHFLOWS$ISIN[-cf_isin_index]
		bdata[[gr]]$CASHFLOWS$CF <- bdata[[gr]]$CASHFLOWS$CF[-cf_isin_index]
		bdata[[gr]]$CASHFLOWS$DATE <- bdata[[gr]]$CASHFLOWS$DATE[-cf_isin_index]
	
	bdata
}


