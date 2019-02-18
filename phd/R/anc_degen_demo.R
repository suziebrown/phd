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

N.paris<-20
Ntilde<-2
par.init<-c(0,1)
par.trans<-c(0.7,0.2)
par.em<-c(1,1)
Nobs <- 40

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

y<-obs.sim(Nobs, par.init, par.trans, par.em)

htilde <- function(x, x_new){
  rep(x_new, length(x))
}

out2 <- PF.smooth(N.paris, y, par.init, par.trans, par.em)
ksi <- out2$ksi


xlim1 <- c(1,Nobs)
ylim2 <- range(ksi)

plot(NA,ylim=ylim2, xlim=xlim1, xlab="time", ylab="position")
#for (i in 1:Nobs){
#  points(rep(i, N.paris), ksi[i,], pch=16, col="purple")
#}
for (i in 1:N.paris){
  lines(ksi[,i], col="purple")
}
points(rep(Nobs, N.paris), ksi[Nobs,], pch=16, col="purple")
