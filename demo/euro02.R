data(eurobonds)

group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- "all"  

x <- splines_estim(group, bonddata, matrange)

print(x)
summary(x)
plot(x)
