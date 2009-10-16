
###################################################################
#                print-method for nelson                          #
###################################################################

print.nelson <- 
  function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Parameters for Nelson/Siegel, Svensson estimation:\n")
  cat("\n")
  cat("Method:",x$method,"\n")
  cat("Fitted:",x$fit,"\n")
  cat("Weights:",x$weights,"\n")
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  parameters <- mapply(function(i) x$opt_result[[i]]$par,seq_along(x$opt_result))
  colnames(parameters) <- names(x$opt_result)
  n_par <- as.character(nrow(parameters))
  rownames(parameters) <- switch(n_par,
          "4"=c("beta_0","beta_1","beta_2","tau_1"),
          "6"=c("beta_0","beta_1","beta_2","tau_1","beta_3","tau_2")) 
  print.default(parameters)
  cat("\n")
  #x
  }

###################################################################
#                    plot-method for nelson                       #
###################################################################

plot.nelson <-
  function(x,matrange=c(min(mapply(function(i) min(x$y[[i]][,1]),seq(x$n_group))),
                        max(mapply(function(i) max(x$y[[i]][,1]),seq(x$n_group))))
                        ,multiple=FALSE, expoints=unlist(x$expoints), ctype="spot",
                         errors="price",
                        lwd=2,lty=1,type="l",inset=c(0.8,0.1),ask=TRUE,
                        ...) {
     
     # min and max maturity of all bonds in the sample 
     samplemat <- c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                    max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group)))) 
  
     # check plot maturity conformity
    if(x$matrange[1] != "all") {
    if(matrange[2]>  x$matrange[2]) { matrange[2] <-  x$matrange[2]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
   
    if(matrange[1] <  x$matrange[1]) { matrange[1] <-  x$matrange[1]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
    }
   
    if( matrange[2] > samplemat[2]) {matrange[2] <-  samplemat[2]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
     
    if( matrange[1] < samplemat[1]) {matrange[1] <- samplemat[1]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
    
               				
    
    
    cdata <- switch(ctype, "spot" = x$spot,
    					   "forward" = x$forward,
    					   "discount" = x$discount
    					    )
    					   
    cname <- switch(ctype, "spot" = "Zero-coupon yield curve",
    					   "forward" = "Forward rate curve",
    					   "discount" = "Discount factor curve" )
    
    
    # plot all interst rate curves together
    if (multiple) {
    
    plot(x=cdata,multiple=multiple, expoints=expoints,lwd=lwd,type=type,...) }
  
    if (!multiple && ctype %in% c("spot", "forward", "discount")){
        old.par <- par(no.readonly = TRUE)

        if(x$n_group != 1) par(ask=ask)
        
    # plot each interest rate curve seperately
    for (k in seq(x$n_group)  ) 
    	{
    	
    	plot.ir_curve(cdata[[k]], ylim=c(0, max(cdata[[k]][,2]) + 0.01 )*100,
    	xlim=c(max(floor(min(x$y[[k]][,1])),matrange[1]),
             min(ceiling(max(x$y[[k]][,1])),matrange[2])), lwd=lwd,type=type, ...
    	)
    	 
    	title(names(x$opt_result)[k])
    	 
    	if(ctype=="spot") {points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red") 
    		legend("bottom",legend=c("Zero-coupon yield curve","Yield-to-maturity"),
                col=c("steelblue","red"), lty = c(1, -1), pch=c(-1,21))}
        else 	legend("bottom",legend=cname	,col=c("steelblue"), lty = lty , pch=(-1))

    	
    	}     
 	    on.exit(par(old.par))
     }
    
     # plot spread curves 
    if(ctype == "spread") {plot(x$spread,expoints=expoints,
    	xlim= c(max(floor(samplemat[1]),matrange[1]),
  	min(ceiling(samplemat[2]),matrange[2])),lwd=lwd,
  						    ...)
  	}
    # plot errors 
    if(errors %in% c("price", "yield")){
    	
    	edata <- switch(errors,"price" = x$perrors, "yield"= x$yerrors )
    	if(x$n_group == 1) ask= FALSE
       	for(k in seq(x$n_group)){
    		
     		plot.error(edata[[k]],ask=ask,main=x$group[k],
                           ylab=paste("Error ",paste(errors,"s)",sep=""),sep=" ("),...)
    		
    		legend("bottomright", legend=c(paste("  RMSE",
    		        switch(errors,"price" = round(rmse(x$p[[k]],x$phat[[k]]),4),
                       "yield" = round(rmse(x$y[[k]][,2],x$yhat[[k]][,2]),4)) ,sep=": "),
                        paste("AABSE",switch(errors,"price" = round(aabse(x$p[[k]],
                        x$phat[[k]]),4),
                        "yield" = round(aabse(x$y[[k]][,2],x$yhat[[k]][,2]),4)),
                        sep=": ")),bty="n", inset=inset) 
    		
    	}
    	
      }
    				
   
   
}  

###################################################################
#                 summary-method for nelson                       #
###################################################################

summary.nelson <-
    function(object,...) {
    x <- object
    RMSE_p <- mapply(function(i) rmse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    AABSE_p <- mapply(function(i) aabse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    RMSE_y <- mapply(function(i) rmse(x$y[[i]][,2],x$yhat[[i]][,2]),seq(x$n_group))
    AABSE_y <- mapply(function(i) aabse(x$y[[i]][,2],x$yhat[[i]][,2]),seq(x$n_group))
    
    gof <- rbind(RMSE_p,AABSE_p,RMSE_y,AABSE_y)
    colnames(gof) <- names(x$p)
    rownames(gof) <- c("RMSE-Prices","AABSE-Prices","RMSE-Yields","AABSE-Yields")
    convergencegroup <- as.matrix(apply(as.matrix(mapply(function(i) 
                              x$opt_result[[i]]$convergence,
                              seq_along(x$opt_result))),1,
                              function(x) if(x==1) "no convergence" else "converged"))
    colnames(convergencegroup) <- "Convergence ()"
    rownames(convergencegroup) <- x$group
    convergence <- as.matrix(mapply(function(i) x$opt_result[[i]]$message,
                   seq_along(x$opt_result)))
    colnames(convergence) <- "Solver message"
    rownames(convergence) <- x$group
    sumry <- list(gof,convergencegroup,convergence)
    names(sumry) <- c("gof","convergencegroup","convergence")
    class(sumry) <- "summary.nelson"
    sumry
}

###################################################################
#                 print-method for summary.nelson                 #
###################################################################

print.summary.nelson <-
    function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$gof)
    cat("\n")
    cat("\n")
    cat("---------------------------------------------------\n")
    cat("Convergence information:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$convergencegroup)
    cat("\n")
    print.default(x$convergence)
    cat("\n")
    cat("\n")
    x
}
 
###################################################################
#                print-method for cubic splines                   #
###################################################################

print.cubicsplines <- 
  function(x,...) {
  cat("---------------------------------------------------\n")
  cat("Parameters for Cubic splines estimation:\n")
  cat("\n")
  for(i in seq(x$n_group)) {
  print.default(paste(names(x$alpha)[[i]],":",sep=""))
  names(x$alpha[[i]]) <- paste("alpha",c(seq_along(x$alpha[[i]])))
  print.default(x$alpha[[i]])
  cat("\n")
  x
  }
 }
 
###################################################################
#            summary-method for cubic splines                      #
###################################################################

summary.cubicsplines <-
    function(object,...) {
    x <- object
    RMSE_p <- mapply(function(i) rmse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    AABSE_p <- mapply(function(i) aabse(x$p[[i]],x$phat[[i]]),seq(x$n_group))
    RMSE_y <- mapply(function(i) rmse(x$y[[i]][,2],x$yhat[[i]][,2]),seq(x$n_group))
    AABSE_y <- mapply(function(i) aabse(x$y[[i]][,2],x$yhat[[i]][,2]),seq(x$n_group))
    gof <- rbind(RMSE_p,AABSE_p,RMSE_y,AABSE_y)
    colnames(gof) <- names(x$p)
    rownames(gof) <- c("RMSE-Prices","AABSE-Prices","RMSE-Yields","AABSE-Yields")
    regsumry <- lapply(x$regout,summary)
    for (i in seq(x$n_group)) rownames(regsumry[[i]]$coefficients) <- 
	paste("alpha",c(seq_along(x$alpha[[i]])))
    sumry <- list(gof,regsumry)
    names(sumry) <- c("gof", "regsumry")
    class(sumry) <- "summary.cubicsplines"
    sumry
} 

###################################################################
#            print-method for summary.cubicsplines                #
###################################################################

print.summary.cubicsplines <-
    function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$gof)
    cat("\n")
    x$gof
    
    cat("---------------------------------------------------\n")
    cat("Summary statistics for the fitted models:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x$regsumry)
    cat("\n")
    x$regsumry

    
}

###################################################################
#                plot-method for cubic splines                    #
###################################################################

plot.cubicsplines <-
  function(x,matrange =c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                        max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group)))),
                        multiple=FALSE, ctype="spot",
                        lwd=2,lty=1,type="l",errors="price",inset=c(0.8,0.1),ask=TRUE, ...) {
       
     # min and max maturity of all bonds in the sample 
     samplemat <- c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                    max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group)))) 
  
     # check plot maturity conformity
    if(x$matrange[1] != "all") {
    if(matrange[2]>  x$matrange[2]) { matrange[2] <-  x$matrange[2]
       warning("The plot range for the maturity violates the estimation maturity range") 
    }
   
    if(matrange[1] <  x$matrange[1]) { matrange[1] <-  x$matrange[1]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
    }
   
    if( matrange[2] > samplemat[2]) {matrange[2] <-  samplemat[2]
       warning("The plot range for the maturity violates the estimation maturity range") 
     }
     
    cdata <- switch(ctype, "spot" = x$spot,
    					   "forward" = x$forward,
    					   "discount" = x$discount
    					    )
    					   
    cname <- switch(ctype, "spot" = "Zero-coupon yield curve",
    					   "forward" = "Forward rate curve",
    					   "discount" = "Discount factor curve" )
    
    
    # plot all interst rate curves together
    if (multiple) {
    
    plot(x=cdata,multiple=multiple, expoints=NULL,lwd=lwd,type=type,...) }
  
	 if (!multiple && ctype %in% c("spot", "forward", "discount")){
        old.par <- par(no.readonly = TRUE)
        if(x$n_group !=1) par(ask=ask)
        
    	# plot each interest rate curve seperately
    	for (k in seq(x$n_group)  ) 
    	{
    	
    	plot.ir_curve(cdata[[k]], ylim=c(0, max(cdata[[k]][,2]) + 0.01 )*100,
    	xlim=c(max(floor(min(x$y[[k]][,1])),matrange[1]),
             min(ceiling(max(x$y[[k]][,1])),matrange[2])), lwd=lwd,type=type,...)
    	
    	
    	 
    	title(x$group[k])
    	 
    	if(ctype=="spot") {points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red") 
    	  # lower ci         
          lines(cdata[[k]][,1],cdata[[k]][,3]*100, type="l", lty=3, col="steelblue" )   
          # upper ci 
          lines(cdata[[k]][,1],cdata[[k]][,4]*100, type="l", lty=3, col="steelblue")
    	  # knot points 
    	  abline(v=c(x$knotpoints[[k]]),lty=2, col="darkgrey")
    	  legend("bottom",legend=c("Zero-coupon yield curve",
    	  "95 % Confidence interval" ,"Yield-to-maturity", "Knot points"),
    	  col=c("steelblue","steelblue","red", "darkgrey"),
    	  lty = c(1,3,-1,2), pch=c(-1,-1,21,-1))
	
    	     } else  legend("bottom",legend=cname,col=c("steelblue"), lty = lty , pch=(-1))

    	
    	 }
        on.exit(par(old.par))
 	}
    	
    # plot spread curves 
    if(ctype == "spread") {plot(x$spread,expoints=NULL,
    	xlim= c(max(floor(samplemat[1]),matrange[1]),
  	     min(ceiling(samplemat[2]),matrange[2],max(mapply(function(i) 
	     max(x$spread[[i]][,1]),seq(x$spread))))),lwd=lwd ,...) 
       }
    						
    						
     # plot errors 
    if(errors %in% c("price", "yield")){
    	
    	edata <- switch(errors,"price" = x$perrors, "yield"= x$yerrors )

        if(x$n_group == 1) ask= FALSE

        for(k in seq(x$n_group)){
     		plot.error(edata[[k]],ask=ask
                ,main=x$group[k],ylab=paste("Error ",paste(errors,"s)",sep=""),sep=" ("),...)
    		
    		legend("bottomright", legend=c(paste("  RMSE",
    		switch(errors,"price" = round(rmse(x$p[[k]],x$phat[[k]]),4),
                       "yield" = round(rmse(x$y[[k]][,2],x$yhat[[k]][,2]),4)) ,sep=": "),
                        paste("AABSE",switch(errors,"price" = round(aabse(x$p[[k]],
                        x$phat[[k]]),4),
                        "yield" = round(aabse(x$y[[k]][,2],x$yhat[[k]][,2]),4)),
                        sep=": ")),bty="n", inset=inset) 
    		
    	  }
    	
     }						
    					     					        							
}  

###################################################################
#                    plot-method for ir_curve                     #
###################################################################
plot.ir_curve <- function(x,ylim=c(),xlim=c(),lwd=2, type="l",
        xlab="Maturity (years)",ylab="Percent", 
  	col="steelblue",lty=1, ...) 
				{
	plot(x[,1] ,x[,2]*100, type=type, ylim=ylim, xlim=xlim, xlab=xlab,
     ylab=ylab,lwd=lwd,lty=lty,col=col, ... )
      
}

###################################################################
#                    plot-method for spot_curves                  #
###################################################################

plot.spot_curves <- function(x,multiple= FALSE,
		   ylim= c(range(mapply(function(i) 
		   range(x[[i]][,2]),seq(x))))*100,xlim=c(),
                   type="l", lty=1, lwd=2, expoints=NULL, 
                   ylab= "Zero-coupon yields (percent)",
                   xlab= "Maturity (years)",main="Zero-coupon yield curves",
                            ...) {

	if(multiple) 
	{ plot(x[[which.max(mapply(function(i) max(x[[i]][,1]),
		seq(x)))]][,1], x[[which.max(mapply(function(i) 
		max(x[[i]][,1]), seq(x)))]][,2]*100, 
        type=type,col=which.max(mapply(function(i) max(x[[i]][,1]),seq(x))),
        lty=lty,lwd=lwd,xlab=xlab,ylab=ylab,ylim=ylim, ... )
   
	  for(k in c((seq(x))[-which.max(mapply(function(i) max(x[[i]][,1]), seq(x)))]))
	  { lines(x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) 
	  	else seq(nrow(x[[k]]))),1],x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) 
	  	else seq(nrow(x[[k]]))),2]*100,col=k,lwd=lwd,lty=lty, ... )
      
        if(is.numeric(expoints))
        {lines(x[[k]][((expoints[k]+1):nrow(x[[k]])) ,1],
        	x[[k]][((expoints[k]+1):nrow(x[[k]])),2]*100,col=k,lwd=lwd,lty=5, ... )
         }
    title(main)
    legend("bottom",legend=names(x),col=seq(x),lty=lty,lwd=lwd)
   }
  }
	else
	{   old.par <- par(no.readonly = TRUE)
		par(ask=TRUE)
		for(k in seq(x)) 
		{ 
		plot.ir_curve(x[[k]],...)
		title(names(x)[k])
      	legend("bottom",legend=main,col=c("steelblue"), lty = 1 , pch=c(-1))
        }	
        on.exit(par(old.par))	
	}
	

}

###################################################################
#                    plot-method for fwr_curves                   #
###################################################################

plot.fwr_curves <- function(x,multiple= FALSE,
   ylim= c(range(mapply(function(i) range(x[[i]][,2]),
   seq(x))))*100,xlim=c(),type="l", lty=1, 
   lwd=2, expoints=NULL, ylab= "Forward rate (percent)",
   xlab= "Maturity (years)",main="Forward rate curves",...) 
					
{ plot.spot_curves(x,ylab=ylab, xlab=xlab, main=main,
	multiple=multiple,expoints=expoints,lty=lty,lwd=lwd,type=type, ... )

}

###################################################################
#                    plot-method for df_curves                   #
###################################################################

plot.df_curves <- function(x,multiple= FALSE,
	ylim= c(range(mapply(function(i) range(x[[i]][,2]),
	seq(x))))*100,xlim=c(),type="l", lty=1,
	lwd=2, expoints=NULL, ylab="Discount factor (percent)",
	xlab= "Maturity (years)",main="Discount factor curves",...) 
{ plot.spot_curves(x,ylab=ylab, xlab=xlab, main=main,
	multiple=multiple,expoints=expoints,lty=lty,lwd=lwd,type=type, ... )

		}
	
###################################################################
#                    plot-method for s_curves                     #
###################################################################
plot.s_curves <- function(x,xlim=c(range(mapply(function(i) 
	   range(x[[i]][,1]),seq(x)))),
	   ylim=c(range(mapply(function(i) range(x[[i]][,2]),
	   seq(x))))*10000,expoints=NULL, xlab="Maturity (years)", 
	   ylab="Spread (basis points)", lwd=2,lty=1, main="Spread curves", ...)
					
{  if(!is.character(x))
   {
	
   plot(0,0, type="n",xlab=xlab,ylab=ylab,xlim=xlim, ylim=ylim,...)
 
   for(k in c(2:length(x)))
    { lines(x[[k]][(if(is.numeric(expoints))
    	seq(expoints[k]) else seq(nrow(x[[k]]))),1],
    	x[[k]][(if(is.numeric(expoints)) seq(expoints[k]) 
    	else seq(nrow(x[[k]]))),2]*10000,col=k,lwd=lwd,lty=lty, ... )
      
      if(is.numeric(expoints))
      {lines(x[[k]][((expoints[k]+1):nrow(x[[k]])) ,1],
     	x[[k]][((expoints[k]+1):nrow(x[[k]])),2]*10000,
     	col=k,lwd=lwd,lty=5, ... )
   
       }
    } 
   title(main)
   legend("topleft",legend=names(x)[-1],col=seq(x)[-1],lty=1,lwd=lwd)
   } else warning("No spread curves available")
} 

	
###################################################################
#              		 plot-method for error                         #
###################################################################    

plot.error <- function(x,type="b",main="", mar= c(7,6,6,2) + 0.1, oma=c(4,2,2,2) +0.1,
                       ylab="Error", ...) {
	old.par <- par(no.readonly = TRUE)
    par(mar=mar, oma=oma, ... )
    
   plot(x[,1],x[,2],axes=FALSE,pch=19,lwd=c(1,2),xlab="", ylab=ylab,type=type, ...)
   axis(1,x[,1],rownames(x),las=3,...)
   axis(2,...)
   axis(3,x[,1],round(x[,1],2),...)
   mtext("Maturity (years)",3,line=2.5)
   lines(x[,1],rep(0,nrow(x)),lty=2,lwd=1,... )
   title(xlab="ISIN", outer=TRUE,main=main,...) 
	 
   on.exit(par(old.par))
}               
	
