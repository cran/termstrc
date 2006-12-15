 
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
  x
  }

###################################################################
#                    plot-method for nelson                       #
###################################################################

plot.nelson <-
  function(x,matrange=c(min(mapply(function(i) min(x$y[[i]][,1]),seq(x$n_group))),
                        max(mapply(function(i) max(x$y[[i]][,1]),seq(x$n_group))))
                        ,pdf=FALSE, ...) {
   
    
     if(pdf) pdf( file="termstrc_results.pdf",... ) else par(ask=TRUE)  
     
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
                       				
    # plot each yield curve seperately
    for (k in seq(x$n_group)  ) {
         
      plot(x$zcy_curves[,1] ,x$zcy_curves[,k+1]*100,
      type="l",
      ylim=c(0, max(x$y[[k]][,2]) + 0.01 )*100,
      xlim=c(max(floor(min(x$y[[k]][,1])),matrange[1]),
             min(ceiling(max(x$y[[k]][,1])),matrange[2])),
      xlab="Maturity (years) ",
      ylab="Percent",
      lwd=2,
      col="steelblue")
      title(names(x$opt_result)[k])
      legend("bottomright",legend=c("Zero-coupon yield curve","Yield to maturity"),
              col=c("steelblue","red"), lty = c(1, -1), pch=c(-1,21))
      grid()
      points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red") 
      
    }
    
    # plot all zero coupon yield curves together
    if (is.numeric(x$scurves)) {
     plot(x$zcy_curves[,1], x$zcy_curves[,
      which.max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group))) +1 ]*100,
      type="l",col=which.max(mapply(function(i) max(x$y[[i]][,1]),seq(x$n_group))),
      lty=1,lwd=2,xlab="Maturity (years)",
      ylab="Zero-Coupon yields (%)",
      xlim=c(max(floor(samplemat[1]),matrange[1]), 
             min(ceiling(samplemat[2]),matrange[2])),
      ylim= c(0,max(x$zcy_curves[,2:(x$n_group+1)] ))*100)
   
	  for(k in c((seq(x$n_group))[-which.max(mapply(function(i) max(x$y[[i]][,1]),
                                         seq(x$n_group)))]))
	  {spoint <- which(x$zcy_curves[,1] > 
                      mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group))[k])[1] 
     lines(x$zcy_curves[1:spoint ,1],x$zcy_curves[1:spoint,k+1]*100,col=k,lwd=2)
     lines(x$zcy_curves[((spoint+1) : nrow(x$zcy_curves) ) ,1],
     x$zcy_curves[((spoint+1) : nrow(x$zcy_curves)),k+1]*100,col=k,lty=5,lwd=2)
 	  } 
    title("Zero-coupon yield curves")
    legend("bottomright",legend=names(x$opt_result),col=seq(x$n_group),lty=1,lwd=2)
    grid()
   }
                                        
    # plot spread curves    
   if (is.numeric(x$scurves)) {
    plot(0,0, type="n",
    col=(which.max(mapply(function(i) max(x$y[[i]][,1]),
                                         seq(x$n_group))[-1]) + 1),lty=1,lwd=2,
    xlab="Maturity (years)",
    ylab="Spread (basis points)",
    xlim= c(max(floor(samplemat[1]),matrange[1]),
            min(ceiling(samplemat[2]),matrange[2])),
    ylim=c(min(x$scurves[,seq(x$n_group)-1]),max(x$scurves[,seq(x$n_group)-1]))*10000)
    
    for(k in c(2:x$n_group))
    {spoint <- which(x$zcy_curves[,1] > mapply(function(i) max(x$y[[i]][,1]),seq(x$n_group))[k])[1] 
     lines(x$zcy_curves[1:spoint ,1],x$scurves[1:spoint ,k-1]*10000,col=k,lwd=2)
     lines(x$zcy_curves[((spoint+1) : nrow(x$zcy_curves) ) ,1],
     x$scurves[((spoint+1) : nrow(x$zcy_curves) ) ,k-1]*10000,col=k,lty=5,lwd=2)
 	  } 
    title("Spread curves")
    legend("topleft",legend=names(x$opt_result[-1]),col=2:x$n_group,lty=1,lwd=2)
    grid()
   
   }                    
  
   
   if(pdf) dev.off() else par(ask=FALSE) 
   
   
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
    convergencegroup <- as.matrix(apply(as.matrix(mapply(function(i) x$opt_result[[i]]$convergence,
                              seq_along(x$opt_result))),1,
                              function(x) if(x==1) "no convergence" else "converged"))
    colnames(convergencegroup) <- "Convergence ()"
    rownames(convergencegroup) <- x$group
    convergence <- as.matrix(mapply(function(i) x$opt_result[[i]]$message,seq_along(x$opt_result)))
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
    cat("Goodness of fit tests:\n")
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
    class(gof) <- "summary.cubicsplines"
    gof
} 

###################################################################
#            print-method for summary.cubicsplines                #
###################################################################

print.summary.cubicsplines <-
    function(x,...) {
    cat("---------------------------------------------------\n")
    cat("Goodness of fit tests:\n")
    cat("---------------------------------------------------\n")
    cat("\n")
    print.default(x)
    cat("\n")
    x
}

###################################################################
#                plot-method for cubic splines                    #
###################################################################

plot.cubicsplines <-
  function(x,matrange =c(min(mapply(function(i) min(x$y[[i]][,1]), seq(x$n_group))),
                        max(mapply(function(i) max(x$y[[i]][,1]), seq(x$n_group))))
                        ,pdf=FALSE, ...) {
  
     if(pdf) pdf( file="termstrc_results.pdf",... )  else par(ask=TRUE) 
     
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
                   				
    # plot each zero cupon yield curve seperately
    for (k in seq(x$n_group)  ) {  
      plot(x$zcy_curves[[k]][,1] ,x$zcy_curves[[k]][,2]*100,
      type="l",
      ylim=c(0,max(x$zcy_curves[[k]][,2]) + 0.01 )*100,
      xlim=c(max(0,matrange[1]),min(max(x$zcy_curves[[k]][,1]),matrange[2])),
      xlab="Maturity (years)",
      ylab="Percent",
      lwd=2,
      col="steelblue")
      title(x$group[k])
      legend("bottomright",legend=c("Zero-coupon yield curve","Yield to maturity"),
              col=c("steelblue","red"), lty = c(1, -1), pch=c(-1,21))
      grid()
      points(x$y[[k]][,1],x$y[[k]][,2]*100,col="red")
      abline(v=c(x$T[[k]]),lty=2, col="darkgrey") 
      if(pdf == FALSE) par(ask=TRUE) 
    }
    
    
    # plot all zero cupon yield curves together
    if (is.numeric(x$scurves)){
    plot(x$zcy_curves[[which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]),
                                 seq(x$n_group)))]][,1],
      		x$zcy_curves[[which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]),
                                  seq(x$n_group)))]][,2]*100,
     type="l",
     ylim=c(0,
                max(x$zcy_curves[[which.max(mapply(function(k) 
                max(x$zcy_curves[[k]][,1]), seq(x$n_group)))]][,2])+ 0.01 )*100,
     xlim=c(max(0,matrange[1]),min(matrange[2],
              max(x$zcy_curves[[which.max(mapply(function(k) 
              max(x$zcy_curves[[k]][,1]), seq(x$n_group)))]][,1]))),
     xlab="Maturity (years)",
     ylab="Percent",
     lwd=2,
     col=which.max(mapply(function(k) max(x$zcy_curves[[k]][,1]), seq(x$n_group))))
     grid()
     title("Zero coupon yield curves") 
      
	  for(k in c( (seq(x$n_group))[- which.max(mapply(function(k) 
              max(x$zcy_curves[[k]][,1]), seq(x$n_group))) ]))
	  {lines(x$zcy_curves[[k]][,1] ,
      		x$zcy_curves[[k]][,2]*100, lwd=2,col=k )
		}
	  legend("bottomright",legend=x$group,col=seq(x$n_group), lty=1, lwd=2)	
    }
    
    # plot spread curves    
    if (is.numeric(x$scurves)) {
    matplot(x$zcy_curves[[which.min(mapply(function(k)
             max(x$zcy_curves[[k]][,1]), seq(x$n_group)))]][,1],
             x$scurves[,1:(x$n_group-1)]*10000, type="l",
    col=2:x$n_group,lty=1,lwd=2,
    xlab="Maturity (years)",
    ylab="Spread (basis points)",
    xlim= c(max(0,matrange[1]),
            min(max(x$zcy_curves[[which.min(mapply(function(k)
            max(x$zcy_curves[[k]][,1]), seq(x$n_group)))]][,1]),matrange[2])),
    ylim=c(min(x$scurves),max(x$scurves ))*10000)
    title("Spread curves")
    legend("topleft",legend=x$group[-1],col=2:x$n_group,lty=1,lwd=2)
    grid()
    }                    
  
   if(pdf) dev.off() else par(ask=FALSE) 
}  
