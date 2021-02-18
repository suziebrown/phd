### Plot a realisation of n-coalescent
## This actually works now!

n <- 50

rates <- (2:n)*(2:n -1)/2 # calculate {k choose 2} for k=1,...,n-1 (the exponential rates of the waiting times)
times <- rev(rexp(n-1,rates)) # holding times (going up the tree)
sumtimes <- c(0,cumsum(times)) # event times (going up the tree)

partner <- rev(round(runif(n-1)*(1:(n-1)) + 0.5 )) # index of lefthand partner in coalescing pair (going up the tree)
partner <- sample(n=1:(n-1))

plot(0,type='n',axes=FALSE,ann=FALSE, xlim=c(0,n), ylim=c(0,sumtimes[n])) # blank plot

pos <- 1:n # reset pos values (pos = where each vertical line is located in x direction)
for (i in 1:(n-1)){
  segments(x0=pos, y0=rep(sumtimes[i],n), x1=pos, y1=rep(sumtimes[i+1],n)) # add vertical segments up to next event time
  segments(x0=pos[partner[i]], y0=sumtimes[i+1], x1=pos[partner[i]+1], y1=sumtimes[i+1]) # add horizontal segment for coalescing pair
  pos[partner[i]] <- mean(c(pos[partner[i]], pos[partner[i] +1])) # new pos is in middle of two coalescing positions
  pos <- pos[-(partner[i]+1)] # remove the extra value from pos; no. of lineages goes down by one.
}
