## distance-attraction model of coalescence

N <- 10  ## total number of particles
maxmove <- 0.05 ## max distance a particle can move in x/y direction at each time step
T = 10  ## no. time steps to simulate

## initial positions
t <- 0
posx <- runif(N)
posy <- runif(N)
plot(posx, posy, xlim=c(0,1), ylim = c(0,1), main=paste("t = ",t), pch=16, asp=1)

for (t in 1:T){
  posx <- posx + runif(N, -maxmove, maxmove)
  posy <- posy + runif(N, -maxmove, maxmove)
  
  ## periodic boundaries
  posx[posx<0] <- posx[posx<0] +1
  posx[posx>1] <- posx[posx>1] -1
  posy[posy<0] <- posy[posy<0] +1
  posy[posy>1] <- posy[posy>1] -1
  
  plot(posx, posy, xlim=c(0,1), ylim = c(0,1), main=paste("t = ",t), pch=16)
  #points(posx, posy, cex=5)  ## add "interaction radius" around points
  Sys.sleep(0.5)  ## "animate" the plot
}
  
  