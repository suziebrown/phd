### PLOTS ILLUSTRATING KINGMAN COALESCENT PROPERTIES
#library(ggplot2)

nmax<-100
nvals <- seq(2,nmax)

ET1 <- function(n){2*(1-(1/n))}
VT1 <- function(n){4* cumsum(n^(-2)*(n-1)^(-2))}
EL <- function(n){2* cumsum((n-1)^(-1))}
VL <- function(n){4* cumsum((n-1)^(-2))}

plot(nvals, ET1(nvals), type='l')
lines(nvals, VT1(nvals), lty=2)
plot(nvals, EL(nvals), type='l')
lines(nvals, VL(nvals), lty=2)


#plotdata <- data.frame(n=nvals, ET1=ET1(nvals), VT1=VT1(nvals), lower=ET1(nvals)-sqrt(VT1(nvals)), upper=ET1(nvals)+(VT1(nvals)))
#ggplot(plotdata) + geom_line(aes(x=n, y=ET1)) + geom_ribbon(aes(x=n, ymin=lower, ymax=upper), alpha=0.3)

