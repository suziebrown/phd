### I used the following code to fiddle around with values of w and u
#   until I found values that illustrate things well, for use in the
#   illustrative diagrams of Background>Resampling section in thesis.

using Plots

N = 6

w = fill(0.0, N)
u = fill(0.0, N)

# set values in 6*simplex for weights:
w[1] = 1.5
w[2] = 0.3
w[3] = 0.6
w[4] = 2.1
w[5] = 1.2
w[6] = 0.3 #N - sum(w)
println("w = ", round.(w/N, digits=2))
grid = [0.0; cumsum(w)] # ends of weight intervals
gridsorted = [0.0; cumsum(sort(w))]

u[3] = 0.27
u[2] = 0.29
u[6] = 0.36
u[5] = 0.54
u[1] = 0.78
u[4] = 0.92
println("u = ", round.(u, digits=2))
usyst = u[1] .+ collect(0:(N-1)) # systematic sampling points
ustrat = u + collect(0:(N-1)) # stratified sampling points
u = u * N # multinomial sampling points


# plot with weights not sorted
vline(0:N, leg=false, color=:lightgray, yaxis=[0,0.4])
vline!(grid, color=:black)
scatter!(u, fill(0.3, N), color=:pink)
scatter!(usyst, fill(0.1, N), color=:red)
scatter!(ustrat, fill(0.2, N), color=:orange)

# plot with weights sorted
scatter(0:N, fill(0.0,N+1), leg=false, color=:white)
vline!(gridsorted, color=:gray)
scatter!(u, fill(0.3, N), color=:pink)
scatter!(usyst, fill(0.1, N), color=:red)
scatter!(ustrat, fill(0.2, N), color=:orange)
