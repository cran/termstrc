data(eurobonds)


group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- c(2,10)
method <- "Nelson/Siegel"
fit <- "prices"
weights <- "none"
control <- list(eval.max=100000)

b <- matrix(c(0.02547394, -0.012162592, -0.02547394,    1,
 			0.02611532, -0.011367422, -0.02611532,    1,
			0.02578871, -0.015207250, -0.02578871,    1),
			nrow=3,ncol=4,byrow=TRUE)
			
rownames(b) <- group

colnames(b) <- c("beta0","beta1","beta2","tau1")

x <- nelson_estim(group, bonddata, matrange, 
                  method, fit, weights, startparam=b,control)

print(x)
summary(x)
plot(x,pdf=TRUE)

