data(eurobonds)

group <- c("GERMANY", "AUSTRIA", "ITALY")
bonddata <- eurobonds
matrange <- c(0, 19)  

x <- splines_estim(group, bonddata, matrange)

print(x)
summary(x)
plot(x,errors="none")
plot(x,multiple=TRUE,errors="none")

