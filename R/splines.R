###################################################################
#                        Cubic splines estimation                 #
###################################################################

splines_estim <-
  function(group,
           bonddata,
           matrange="all"
           ) {

    
  # select given group from bonddata
  bonddata <- bonddata[group]

  # select data according to chosen maturity range
  if (length(matrange)==1) {bonddata <- bonddata }else
   {bonddata <- maturity_range(bonddata,matrange[1],matrange[2]) }

  # number of groups 
  n_group <- length(bonddata) 
  
  # group sequence
  sgroup <- seq(n_group)
    
  # create cashflows matrix
  cf <- lapply(bonddata,create_cashflows_matrix)

  # create cashflows matrix including dirty price (needed for bond yield calculation)
  cf_p <- mapply(function(k) create_cashflows_matrix(bonddata[[k]],include_price=TRUE),
                 sgroup,SIMPLIFY=FALSE)
    
  # create maturities matrix
  m <- lapply(bonddata,create_maturities_matrix)

  # create maturities matrix including zeros (needed for bond yield calculation)
  m_p <- mapply(function(k) create_maturities_matrix(bonddata[[k]],include_price=TRUE),
                sgroup,SIMPLIFY=FALSE)
    
  # calculate dirty prices
  p <- mapply(function(k) bonddata[[k]]$PRICE + bonddata[[k]]$ACCRUED,sgroup,SIMPLIFY=FALSE)
  
  # index for ordering
  positions <- mapply(function(k) order(apply(m[[k]],2,max)),sgroup,SIMPLIFY=FALSE)
  
  # order matrices 
  cf <- mapply(function(k) cf[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  cf_p <- mapply(function(k) cf_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m <- mapply(function(k) m[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m_p <- mapply(function(k) m_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
   
  # calculate bond yields	
  y <- mapply(function(k) bond_yields(cf_p[[k]],m_p[[k]]),
                   sgroup,SIMPLIFY=FALSE)
    
  # Choosing knot points (McCulloch)
  # number of bonds in each group 
  K <- mapply(function(k) ncol(m[[k]]),sgroup,SIMPLIFY=FALSE)  
  # number of basis functions
  s <-  mapply(function(k) round(sqrt(K[[k]])),sgroup,SIMPLIFY=FALSE)
  
  # only perfom spline estimation if number of bonds per group >= 9
  if(sum(s>=3) != length(s))  stop(cat("Estimation aborted:
    For cubic splines estimation more than 9 observations per group are required","\n",
    "Check group(s):", group[which((s>=3)==FALSE)]),
    "\n" )
  
  # only used for knot point finding
  i <- mapply(function(k) 2:(max(2,(s[[k]]-2))),sgroup,SIMPLIFY=FALSE)  
  
  h <-  mapply(function(k) trunc(((i[[k]]-1)*K[[k]])/(s[[k]]-2)),sgroup,SIMPLIFY=FALSE)
             
  theta <- mapply(function(k)((i[[k]]-1)*K[[k]])/(s[[k]]-2)-h[[k]],sgroup,SIMPLIFY=FALSE)

  # knot points
  T <- mapply(function(k) if(s[[k]]>3) c(0,
       apply(as.matrix(m[[k]][,h[[k]]]),2,max)
       + theta[[k]]*(apply(as.matrix(m[[k]][,h[[k]]+1]),2,max)-apply(as.matrix(m[[k]][,h[[k]]]),2,max)),
       max(m[[k]][,ncol(m[[k]])])) else c(0,max(m[[k]][,ncol(m[[k]])])),sgroup,SIMPLIFY=FALSE)
 
 
  # parameter estimation with OLS
  Y <- mapply(function(k) apply(cf_p[[k]],2,sum),sgroup,SIMPLIFY=FALSE)
   
  X <- list()

  # k ... group index
  # j ... column index (bond)
  # sidx ... index for spline function  
 
  for (k in sgroup){  
  X[[k]] <- matrix(NA,ncol(m[[k]]),s[[k]])
                 
   for(sidx in 1:s[[k]]){
   X[[k]][,sidx] <- apply(cf[[k]]*mapply(function(j) gi(m[[k]][,j],T[[k]],sidx,s[[k]]),1:ncol(m[[k]])),2,sum)
   }
  }
  
  alpha <- mapply(function(k) coef(lm(-Y[[k]]~X[[k]]-1)),sgroup,SIMPLIFY=FALSE) # parameter vector
 
  # calculate discount factor matrix 
  dt <- list()
   for (k in sgroup){
   dt[[k]] <- matrix(1,nrow(m[[k]]),ncol(m[[k]]))
   for(sidx in 1:s[[k]]){
    dt[[k]] <- dt[[k]] + alpha[[k]][sidx]* mapply(function(j) gi(m[[k]][,j],T[[k]],sidx,s[[k]]),1:ncol(m[[k]]))
   }
  }  

   
  # calculate estimated prices 
  phat <- mapply(function(k) apply(cf[[k]]*dt[[k]],2,sum),sgroup,SIMPLIFY=FALSE)
  
  # calculate estimated yields 
  yhat <- mapply(function(k) bond_yields(rbind(-phat[[k]],cf[[k]]),m_p[[k]]),sgroup,SIMPLIFY=FALSE)
  
 
  # calculate estimated zero coupon yield curves
  t <- mapply(function(k) seq(0.01, max(T[[k]]),0.01), sgroup,SIMPLIFY=FALSE) 
  
  zcy_curves <- mapply(function(k) cbind(t[[k]],matrix(NA,nrow=length(t[[k]]),ncol=1)), sgroup,SIMPLIFY=FALSE)
        
  dt_zcy <- list()
  for( k in sgroup) {
   dt_zcy[[k]] <- rep(1,length(t[[k]]))
    for(sidx in 1:s[[k]]){  
    dt_zcy[[k]] <- dt_zcy[[k]] + alpha[[k]][sidx]*gi(t[[k]],T[[k]],sidx,s[[k]])
   }
  zcy_curves[[k]][,2] <- -log(dt_zcy[[k]])/t[[k]]          
  }
  
  
  # calculate spread curves              	    
 	if(n_group != 1) {  
   scurves <- as.matrix( mapply(function(k) (zcy_curves[[k]][1:nrow(zcy_curves[[which.min(mapply(function(k) min(length(zcy_curves[[k]][,1])), sgroup))]]),2] -
   zcy_curves[[1]][1:nrow(zcy_curves[[which.min(mapply(function(k) min(length(zcy_curves[[k]][,1])), sgroup))]]),2]), 2:n_group))
   
    } else scurves = "none" 
 

 # return list of results
 result <- list(  group=group,          # e.g. countries, rating classes
                  matrange=matrange,    # maturity range of bonds
                  n_group=n_group,      # number of groups,
                  T=T,                  #knot points 
                  zcy_curves=zcy_curves, # zero coupon yield curves
                  scurves=scurves,      # spread curves
                  cf=cf,                # cashflow matrix
                  m=m,                  # maturity matrix
                  p=p,                  # dirty prices
                  phat=phat,            # estimated prices
                  y=y,                  # maturities and yields
                  yhat=yhat,            # estimated yields
                  alpha=alpha           # cubic splines parameters                             
                 )
                 
  # assign names to results list 
  for ( i in 5:length(result)) names(result[[i]]) <- names(bonddata)
    
  class(result) <- "cubicsplines"
  result
 
}



