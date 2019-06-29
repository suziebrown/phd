## One iteration of PaRIS, going from (ksi(t), omega(t), tau(t)) to (ksi(t+1), omega(t+1), tau(t+1)):

paris.iter <- function(ksi, omega, tau, yt, Ntilde, par.trans, par.em){
  ## ksi, omega, tau passed in from previous iteration
  ## yt the most recent observation (scalar)
  ## param to be passed to PF.iter
  ## Ntilde number of indices to sample in Paris step (precision parameter)
  
  # check inputs okay:
  check.dim<-c(length(ksi),length(omega),length(tau))
  if(all(check.dim==check.dim[1])==FALSE) stop("The length of ksi, omega and tau must be the same!")
  # get number of particles:
  N <- length(ksi)
  tau_new <- numeric(N)
  
  # standard particle update:
  PFout <- PF.iter(ksi, omega, yt, par.trans, par.em)
  ksi_new <- PFout$ksi
  omega_new <- PFout$omega
  
  for (i in 1:N){
    # choose indices:
    J <- sample(1:N, Ntilde, replace=T, prob = omega*d.trans(ksi, ksi_new[i], par.trans)) 
    # update particle estimate of h:
    tau_new[i] <- (sum(tau[J]) + sum(htilde(ksi[J], ksi_new[i]))) / Ntilde 
  }
  
  # return:
  list(ksi=ksi_new, omega=omega_new, tau=tau_new)
}



## Wrapper for full PaRIS algorithm when implemented offline:

paris <- function(y, N, Ntilde, par.init, par.trans, par.em){
  ## y=(y_1, ..., y_Nobs) vector of observations from a HMM
  ## param to be passed to PF.iter
  ## N number of particles
  ## Ntilde number of indices to sample in Paris step (precision parameter)
  
  # number of time steps:
  Nobs <- length(y)
  ksi <- matrix(NA, Nobs, N)
  omega <- matrix(NA, Nobs, N)
  tau <- matrix(NA, Nobs, N)
  
  # intialisation:
  ksi[1,] <- r.init(N, par.init)
  omega[1,] <- d.em(y[1], ksi[1,], par.em)
  omega[1,] <- omega[1,]/sum(omega[1,]) 
  tau[1,] <- rep(0,N)
  
  # iterate the ol' Paris:
  for (i in 1:(Nobs-1)){
    parisout <- paris.iter(ksi[i,], omega[i,], tau[i,], y[i], Ntilde, par.trans, par.em)
    ksi[i+1,] <- parisout$ksi
    omega[i+1,] <- parisout$omega
    tau[i+1,] <- parisout$tau
  }
  
  # return:
  list(ksi=ksi, omega=omega, tau=tau)
}

obs.sim <- function(Nobs, par.init, par.trans, par.em){ 
  X <- numeric(Nobs) 
  Y <- numeric(Nobs)
  X[1] <- r.init(1, par.init)
  for (t in 1:Nobs){
    Y[t] <- r.em(1, X[t], par.em)
    X[t+1] <- r.trans(1, X[t], par.trans)
  }
  Y
}


##~~~~~~~~~~
## functions for Gaussian DLM

r.trans <- function(N, xt, par.trans){ # transition distribution for x_t|x_{t-1}
  rnorm(N, par.trans[1]*xt, par.trans[2])
}
d.trans <- function(x1, x2, par.trans){ # forward transition density for X
  dnorm(x1, par.trans[1]*x2, par.trans[2])
}
d.em <- function(yt, xt, par.em){ # emission density for y_t|x_t
  dnorm(yt, par.em[1]*xt, par.em[2])
}
r.em <- function(N, xt, par.em){ # emission distribution for y_t|x_t
  rnorm(N, par.em[1]*xt, par.em[2])
}
r.init <- function(N, par.init){ # initial distribution for x_1
  rnorm(N, par.init[1], par.init[2])
}

PF.smooth <- function(N, y, par.init, par.trans, par.em){
  # N: number of particles
  # y=(y_1,...,y_T): observations
  # param: parameters for initial/trasition/emission distributions (in "super Gaussian" case c(mean.x, sigma.x, sigma.y))
  Nobs <- length(y)
  xtilde <- r.init(N, par.init) 
  W1 <- d.em(y[1], xtilde, par.em) 
  W1 <- W1/sum(W1)
  xtilde <- sample(xtilde, N, prob=W1, replace=TRUE)
  Xtilde <- matrix(0, Nobs, N)
  Xtilde[1,] <- xtilde
  W <- matrix(0, Nobs, N)
  W[1,] <- W1
  for(t in 2:Nobs){
    Xtilde[t,] <- r.trans(N, Xtilde[t-1,], par.trans) 
    W[t,] <- d.em(y[t], Xtilde[t,], par.em)
    W[t,] <- W[t,]/sum(W[t,])
    index <- sample(1:N, N, prob=W[t,], replace=TRUE)
    Xtilde <- Xtilde[,index] 
  }
  return(list(ksi=Xtilde, omega=W))
}
# out <- PF.smooth(10, c(1,3,2,3,3,7,10,5,6,11), c(1,10,4))
# out

PF.filter <- function(N, y, par.init, par.trans, par.em){
  # N: number of particles
  # y=(y_1,...,y_Nobs): observations
  # param: parameters for initial/trasition/emission distributions (in "super Gaussian" case c(mean.x, sigma.x, sigma.y))
  Nobs <- length(y)
  xtilde <- r.init(N, par.init)
  W1 <- d.em(y[1], xtilde, par.em)
  W1 <- W1/sum(W1)
  Xtilde <- matrix(0, Nobs, N)
  Xtilde[1,] <- xtilde
  W <- matrix(0, Nobs, N)
  W[1,] <- W1
  for(t in 1:(Nobs-1)){
    index <- sample(1:N, N, replace=TRUE, prob=W[t,])
    Xtilde[t+1,] <- r.trans(N, Xtilde[t,index], par.trans)
    W[t+1,] <- d.em(y[t+1], Xtilde[t+1,], par.em)
    W[t+1,] <- W[t+1,]/sum(W[t+1,])
  }
  return(list(ksi=Xtilde[Nobs,], omega=W[Nobs,]))
}
# out1<-PF.filter(10,c(1,3,2,3,3,7,10,5,6,11),c(0,1),c(0,1),c(0,1))
# out1

PF.iter <- function(ksi, omega, yt, par.trans, par.em){
  # ksi=(ksi_t^1,...,ksi_t^N): states at time t of each particle
  # omega=(omega_t^1,...,omega_t^N): weights at time t for each particle
  # yt=y_t: observation at time t
  # param: parameters for initial/trasition/emission distributions (in "super Gaussian" case c(mean.x, sigma.x, sigma.y))
  N <- length(ksi)
  index <- sample(1:N, N, replace=TRUE, prob=omega)
  Xtilde <- r.trans(N, ksi[index], par.trans)
  W <- d.em(yt, Xtilde, par.em)
  W <- W/sum(W)
  return(list(ksi=Xtilde, omega=W))
}

N.paris<-10
Ntilde<-2
#N.FFBSm<-30
par.init<-c(0,1)
par.trans<-c(0.7,0.08)
par.em<-c(1,1)
Nobs <- 30

y<-obs.sim(Nobs, par.init, par.trans, par.em)

htilde <- function(x, x_new){
  rep(x_new, length(x))
}

out <- paris(y, N.paris, Ntilde, par.init, par.trans, par.em)
tau <- out$tau
out2 <- PF.smooth(N.paris, y, par.init, par.trans, par.em)
ksi <- out2$ksi

par(mfrow=c(1,1))

xlim1 <- c(1,Nobs)

ylim2 <- range(ksi)
plot(NA,ylim=ylim2, xlim=xlim1, xlab="time", ylab="", axes=FALSE)
title(ylab="position", line=0)
axis(side=1, at=c(0,5,10,15,20,25,30))
for (i in 1:Nobs){
  points(rep(i, N.paris), ksi[i,], pch=16, col="purple")
}
for (i in 1:N.paris){
  lines(ksi[,i], col="purple")
}


