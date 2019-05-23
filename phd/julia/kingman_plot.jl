# ---- Kingman n-coalescent example plot ----

using Random
using Distributions
using Plots
using PhyloTrees


N=10

# initialise
exprv = Array{Float64, 1}(undef, N-1)
pair = Array{Int64, 2}(undef, N-1, 2)

for t in 0:(N-2)
    # number of distinct lineages remaining
    nblocks = N-t

    # generate exponential waiting times between mergers
    exprv[t+1] = rand(Exponential((binomial(nblocks,2))^(-1)))

    # choose a pair of indices to merge at random
    pair[t+1,:] = sample(activeblocks, 2, replace=false)
end
expsum = cumsum(exprv)
