### I wrote this script to try to see whether Hol et al's suggestion
#   of a way to generate order stats from Unif[0,1] is actually correct.
#   TBH this experiment was totally inconclusive :-)

using Plots
N=10
nrep=10000
X = Array{Float64,2}(undef,N,nrep)

for t in 1:nrep
    U = rand(N)
    X[N,t] = U[N]^(1/N)
    for i in (N-1):-1:1
        X[i,t] = X[i+1,t] * U[i]^(1/i)
    end
end

histogram(vec(X), leg=false, bins=10)
