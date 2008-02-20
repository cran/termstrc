###################################################################
#               Nelson-Siegel, Svensson  estimation               #
###################################################################
   
nelson_estim <-
  function(group,
           bonddata,
           matrange="all",
           method="Nelson/Siegel",
           fit = "prices",
           weights="duration",
           startparam,
           control=list(eval.max=1000,iter.max= 500)
           
           ) {

  # check inputs  
  if(fit=="yields"&weights!="none"){
  warning("For minimization of yield errors no weights are needed")
  weights <- "none"}
    
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
  # assign ISIN 
  for(k in sgroup) names(p[[k]]) <- bonddata[[k]]$ISIN
  
  # index for ordering
  positions <- mapply(function(k) order(apply(m[[k]],2,max)),sgroup,SIMPLIFY=FALSE)
  
  
  # order matrices 
  cf <- mapply(function(k) cf[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  cf_p <- mapply(function(k) cf_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m <- mapply(function(k) m[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  m_p <- mapply(function(k) m_p[[k]][,positions[[k]]],sgroup,SIMPLIFY=FALSE)
  p <- mapply(function(k) p[[k]][positions[[k]]],sgroup,SIMPLIFY=FALSE)
  
  # calculate bond yields	
  y <- mapply(function(k) bond_yields(cf_p[[k]],m_p[[k]]),
                   sgroup,SIMPLIFY=FALSE)
 
  # calculate duration   
  duration <- mapply(function(k) duration(cf_p[[k]],m_p[[k]],y[[k]][,2]),
                   sgroup,SIMPLIFY=FALSE)
                   
  
  # objective function 
  obj_fct_prices <- function(b) {    # price error minimization
     loss_function(p[[k]],
     	bond_prices(method,b,m[[k]],cf[[k]])$bond_prices,duration[[k]][,3],weights)}
  
  obj_fct_yields <- function(b) {  
    loss_function(y[[k]][,2],bond_yields(rbind(
    -bond_prices(method,b,m[[k]],cf[[k]])$bond_prices,cf[[k]]),m_p[[k]])[,2],duration[[k]][,3],weights)} 
    
  obj_fct <- switch(fit,
                "prices" = obj_fct_prices,
                "yields" = obj_fct_yields)
  
  # lower and upper bounds for estimation   
  lower_bounds <- switch(method,
                       "Nelson/Siegel" = c(0, -Inf, -Inf, 0),
                       "Svensson" = c(0, -Inf, -Inf, 0, -Inf, 0))
 
  upper_bounds <- switch(method,
                       "Nelson/Siegel" = rep(Inf, 4),
                       "Svensson" = rep(Inf, 6))
 
  # calculate optimal parameter vector
  opt_result <- list()

  for (k in sgroup){
    opt_result[[k]] <- 
                  nlminb(startparam[k,],obj_fct, lower = lower_bounds,
                  upper = upper_bounds,control=control)
  }   
  
  # theoretical bond prices with estimated parameters
  phat <- mapply(function(k) bond_prices(method,opt_result[[k]]$par,
       m[[k]],cf[[k]])$bond_prices,sgroup,SIMPLIFY=FALSE)
       
  # price errors
  perrors <- mapply(function(k) cbind(y[[k]][,1],phat[[k]] - p[[k]]),sgroup,SIMPLIFY=FALSE)     
  for (k in sgroup) class(perrors[[k]]) <- "error"
  
  # calculate estimated yields 
  yhat <- mapply(function(k) bond_yields(rbind(-phat[[k]],cf[[k]]),m_p[[k]]),sgroup,SIMPLIFY=FALSE)
  
  # yield errors
  yerrors <- mapply(function(k) cbind(y[[k]][,1], yhat[[k]][,2] - y[[k]][,2]),sgroup,SIMPLIFY=FALSE)
  for (k in sgroup) class(yerrors[[k]]) <- "error"

  
  # maturity interval
  t <- seq(round(min(mapply(function(i) min(y[[i]][,1]), sgroup)),2),
                           ceiling(max(mapply(function(i) max(y[[i]][,1]), sgroup))),0.01)
  
  
  # calculate zero coupon yield curves  
  zcy_curves <- switch(method,
              "Nelson/Siegel" = mapply(function(k)
		            cbind(t,nelson_siegel(opt_result[[k]]$par,t)),sgroup, SIMPLIFY=FALSE),
              "Svensson" = mapply(function(k) 
              		cbind(t,svensson(opt_result[[k]]$par,t)),sgroup,SIMPLIFY=FALSE)
              		   )

  for (k in sgroup) class(zcy_curves[[k]]) <- "ir_curve"
  class(zcy_curves) <- "spot_curves"
                                      
  # calculate spread curves              	 
   if(n_group != 1) {  
   s_curves <- mapply(function(k) cbind(t,(zcy_curves[[k]][,2] - zcy_curves[[1]][,2])),sgroup,
   					SIMPLIFY=FALSE)
    } else s_curves = "none"
   
   for (k in sgroup) class(s_curves[[k]]) <- "ir_curve" 
   class(s_curves) <- "s_curves"
    
  # calculate extrapolation point                        
  expoints <- mapply(function(k) which(zcy_curves[[k]][,1] > 
                 mapply(function(i) max(y[[i]][,1]), seq(n_group))[k])[1],sgroup, SIMPLIFY=FALSE )  
        
  # calculate forward rate curves 
  fwr_curves <- switch(method,
              "Nelson/Siegel" = mapply(function(k)
		            cbind(t,fwr_ns(opt_result[[k]]$par,t)),sgroup, SIMPLIFY=FALSE),
              "Svensson" = mapply(function(k) 
              		cbind(t,fwr_sv(opt_result[[k]]$par,t)),sgroup,SIMPLIFY=FALSE)
                      )
  for (k in sgroup) class(fwr_curves[[k]]) <- "ir_curve"
  class(fwr_curves) <- "fwr_curves"

  
  # calculate discount factor curves 
  df_curves <- mapply(function(k) cbind(zcy_curves[[k]][,1],exp(-zcy_curves[[k]][,1]*
                                        zcy_curves[[k]][,2])),sgroup,SIMPLIFY=FALSE)
  
   for (k in sgroup) class(df_curves[[k]]) <- "ir_curve"
   class(df_curves) <- "df_curves"
  
   
 
 # return list of results 
 result <- list(group=group,           # e.g. countries, rating classes
                 matrange=matrange,    # maturity range of bonds
                 method=method,        # method (Nelson/Siegel or Svensson)
                 fit=fit,              # fitting method (prices or yields)
                 weights=weights,      # weighting type for estimation
                 n_group=n_group,      # number of groups,
                 spot=zcy_curves,      # zero coupon yield curves
                 spread=s_curves,      # spread curves
                 forward=fwr_curves,   # forward rate curves
                 discount=df_curves,   # discount factor curves
                 expoints=expoints,    # extrapolation points
       		 cf=cf,                # cashflow matrix
                 m=m,                  # maturity matrix
                 duration=duration,    # duration, modified duration, weights
                 p=p,                  # dirty prices
                 phat=phat,            # estimated prices
                 perrors=perrors,      # price errors
                 y=y,                  # maturities and yields
                 yhat=yhat,            # estimated yields
                 yerrors=yerrors,      # yield errors
                 opt_result=opt_result                             
                 )
              
  # assign names to results list 
  for ( i in 7:length(result)) names(result[[i]]) <- names(bonddata)
    
  class(result) <- "nelson"
  result
 }


###################################################################
#                        Spot rates Nelson/Siegel                 #
###################################################################

nelson_siegel <-
  function(beta, m) {
    (beta[1] + beta[2]*((1-exp(-m/beta[4]))/(m/beta[4]))
    + beta[3]*(((1-exp(-m/beta[4]))/(m/beta[4]))-exp(-m/beta[4])))}

###################################################################
#                        Spot rates Svensson                      #
###################################################################

svensson <-
  function(beta, m) {
  (beta[1] + beta[2] * ((1 - exp(-m/beta[4]))/(m/beta[4])) +
  beta[3] * (((1 - exp(-m/beta[4]))/(m/beta[4])) - exp(-m/beta[4])) +
  beta[5] * (((1 - exp(-m/beta[6]))/(m/beta[6])) - exp(-m/beta[6])))}
  
###################################################################
#                     Forward rates Nelson/Siegel                 #
###################################################################

fwr_ns <-
  function(beta, m) {
    (beta[1] + beta[2]*exp(-m/beta[4])
    + beta[3]*(m/beta[4]*exp(-m/beta[4])))}

###################################################################
#                     Forward rates Svensson                      #
###################################################################

fwr_sv <-
  function(beta, m) {
  (beta[1] + beta[2]*exp(-m/beta[4]) +
  beta[3] *m/beta[4]*exp(-m/beta[4]) +
  beta[5] *m/beta[6]*exp(-m/beta[6]))}
  
###################################################################
#                        loss function                            #
###################################################################

loss_function <-
  function(p,phat,omega,weights="none") {
  if (weights=="none") omega <- rep(1,length(p))
  sum(omega*((p-phat)^2))}

	



