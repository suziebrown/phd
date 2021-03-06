### Importance sampling illustration

# shades to match beaver colour theme in beamer
beaverred <- rgb(0.8,0,0) 
beavergrey1 <- rgb(242,242,242)
beavergrey2 <- rgb(217,217,217)

set.seed(127)

x <- seq(-4,4,0.001)
q <- function(x) 0.8*exp(-(0.6*x)^2) # proposal
p <- function(x) 1.2*exp(-(2*(x-1))^2) + 3*exp(-(5*(x))^2) + 0.8*exp(-(x+0.5)^2) # target

par(mar = rep(0.1,4))
plot(0,type='n',axes=FALSE,ann=FALSE, xlim=c(-3,3), ylim=c(0,3.6)) # blank plot
abline(h=0) # zero line
polygon(c(-4,x,4), c(0,p(x),0), col=rgb(0.85,0.85,0.85), lty='blank') # shade target density
lines(x,q(x), lwd=2, lty='dotted') # proposal density

sam <- rnorm(10,0,2) # sample from proposal
points(sam, rep(0,10), pch=16, col=beaverred) # plot sampled points
weights <- p(sam)/q(sam) # calculate importance weights
wsum <- sum(weights)
weights <- weights/wsum # normalise weights
segments(sam,rep(0,10),sam,weights*10, col=beaverred, lwd=2) # add vertical bars for weights

#---
## Resampling illustration

plot(0,type='n',axes=FALSE,ann=FALSE, xlim=c(-2,2), ylim=c(0,0.3)) # blank plot
abline(h=0) # zero line
points(sam, rep(0,10), pch=16, col=beaverred) # plot sampled points
segments(sam,rep(0,10),sam,weights, col=beaverred, lwd=2) # add vertical bars for weights

resam <- rmultinom(1,10,weights)
resam

plot(0,type='n',axes=FALSE,ann=FALSE, xlim=c(-2,2), ylim=c(0,0.3)) # blank plot
abline(h=0) # zero line
# plot crosses on the positions that are not resampled
points(sam[resam==0], rep(0, sum(resam==0)), pch=4, col=beaverred)
# plot a stack of particles above each resampled position
layer <- 0
while (sum(resam)>0) {
  points(sam[resam>0], rep(0.02*layer, sum(resam>0)), pch=16, col=beaverred)
  resam <- resam -1 * (resam>0) # subtract one from each count, or keep equal to zero
  layer <- layer +1 # go up one level on the y axis
}
