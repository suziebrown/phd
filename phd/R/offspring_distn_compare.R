#-----
# Probability density over offspring counts nu, in multinomial and residual-mn resampling

probnu_mn <- function(nu, w){
  # add error catching for suitable w, nu values (same length, w in 0--1 sum to 1, nu integer <=N)
  N <- 3 #length(w)
  factorial(N)*prod(w^nu)/prod(factorial(nu))
}

probnu_res <- function(nu, w){
  # add error catching for suitable w, nu values (same length, w in 0--1 sum to 1, nu integer <=N)
  N <- 3 #length(w)
  flnw <- floor(N*w)
  nubar <- nu - flnw
  wbar <- N*w - flnw
  R <- sum(wbar)
  if (R>0) {wbar <- wbar/R}
  if (all(nubar>=0)){
    out <- factorial(R)*prod(wbar^nubar)/prod(factorial(nubar))
  }
  else{
    out <- 0
  }
  out
}


#-----
# Helper functions for ternary plots

prob111_mn <- function(a, b, c){6*a*b*c}

prob111_res <- function(a, b, c){
  N <- 3 #length(w)
  w <- c(a,b,c)
  flnw <- floor(N*w)
  nubar <- c(1,1,1) - flnw
  wbar <- N*w - flnw
  R <- sum(wbar)
  if (R>0) {wbar <- wbar/R}
  
  all(nubar>=0)*factorial(R)*prod(wbar^nubar)
}

#-----
# Plotting attempt 1: square heatmap (ignore values where x+y>1)

prob2d_mn <- function(x, y) {
  6*x*y*(1-x-y)*(x+y<=1)
}

w1 <- seq(from=0, to=1, by=0.01)
w2 <- seq(from=0, to=1, by=0.01)
mymatrix <- outer(w1,w2, prob2d_mn)
heatmap(mymatrix, Rowv = NA, Colv = NA)


#-----
# Plotting attempt 2: using ternary plot package

library('Ternary')
par(mar=rep(0.2, 4), mfrow=c(1,1))

TernaryPlot(alab = 'w_1', blab = 'w_2', clab = 'w_3', axis.labels = seq(0, 1, by = 0.1), grid.minor.lines=0)
values <- TernaryPointValues(prob111_mn, resolution=48L)
ColourTernary(values)
#TernaryContour(prob111_mn, resolution=10L)

TernaryPlot(alab = 'w_1', blab = 'w_2', clab = 'w_3', axis.labels = seq(0, 1, by = 0.1), grid.minor.lines=0)
values <- TernaryPointValues(prob111_res, resolution=4L)
ColourTernary(values)
#TernaryContour(prob111_res, resolution=100L)
