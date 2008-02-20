data(govbonds)

group <- c("GERMANY","FRANCE", "BELGIUM", "SPAIN")
bonddata <- govbonds
matrange <- c(0,30)
method <- "Nelson/Siegel"
fit <- "prices"
weights <- "duration"
b <- matrix(rep(c(0,0,0, 1),4),nrow=4,byrow=TRUE)
rownames(b) <- group
colnames(b) <- c("beta0","beta1","beta2","tau1")

x <- nelson_estim(group, bonddata, matrange, method, fit, weights, b)

print(x)
summary(x)
plot(x,errors="none")
plot(x,multiple=TRUE,errors="none")
