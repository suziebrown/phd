w1 <- seq(from=0, to=1, by=0.01)
nreps=100

v1mn <- matrix(nrow=length(w1), ncol=nreps)
v2mn <- matrix(nrow=length(w1), ncol=nreps)
v2csmc <- matrix(nrow=length(w1), ncol=nreps)
v1csmc <- matrix(nrow=length(w1), ncol=nreps)

for (j in 1:nreps){
  for (i in 1:length(w1)) {
    res <- rbinom(1, 1, w1[i])
    v1csmc[i,j] <- 1 + res
    v2csmc[i,j] <- 1 - res
  }
}
v1csmcsum = rowSums(v1csmc)/nreps
v2csmcsum = rowSums(v2csmc)/nreps

for (j in 1:nreps){
  for (i in 1:length(w1)) {
    v1mn[i,j] <- rbinom(1, 2, w1[i])
    v2mn[i,j] <- 2 - v1mn[i,j]
  }
}
v1mnsum = rowSums(v1mn)/nreps
v2mnsum = rowSums(v2mn)/nreps

plot(w1, v1csmcsum, ylim=c(0,2), ylab='average of v1')
points(w1, v1mnsum, col='red')

plot(w1, v2csmcsum, ylim=c(0,2), ylab='average of v2')
points(w1, v2mnsum, col='red')

coalcsmc <- v1csmc * (v1csmc-1) + v2csmc * (v2csmc-1)
coalmn <- v1mn *(v1mn-1) + v2mn *(v2mn-1)
coalcsmcsum <- rowSums(coalcsmc)
coalmnsum <- rowSums(coalmn)

plot(w1, coalcsmcsum, ylab='average of cN')
points(w1, coalmnsum, col='red')

